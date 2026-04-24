#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

TEXT_SUFFIXES = {".md", ".html", ".txt"}
TABLE_SUFFIXES = {".tsv", ".csv"}
JSON_SUFFIXES = {".json"}


def parse_args():
    parser = argparse.ArgumentParser(description="Compare baseline and refactored workflow outputs.")
    parser.add_argument("--baseline-root", required=True)
    parser.add_argument("--candidate-root", required=True)
    parser.add_argument("--output-json", required=True)
    parser.add_argument("--output-md", required=True)
    parser.add_argument("--relpath", action="append", default=[], help="Explicit relative paths to compare")
    parser.add_argument("--float-tol", type=float, default=1e-8)
    parser.add_argument("--key-column", action="append", default=[], help="Override key columns as relpath:col1,col2")
    parser.add_argument("--ignore-column", action="append", default=[], help="Ignore columns as relpath:col1,col2")
    return parser.parse_args()


def ensure_parent(path_value: str | Path) -> Path:
    path = Path(path_value)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def discover_paths(root: Path):
    return sorted(
        str(path.relative_to(root)).replace("\\", "/")
        for path in root.rglob("*")
        if path.is_file()
    )


def parse_key_overrides(raw_values):
    mapping = {}
    for raw in raw_values:
        relpath, _, cols = str(raw).partition(":")
        if not relpath or not cols:
            raise ValueError(f"Invalid --key-column: {raw}")
        mapping[relpath.replace("\\", "/")] = [item.strip() for item in cols.split(",") if item.strip()]
    return mapping


def parse_column_mapping(raw_values):
    mapping = {}
    for raw in raw_values:
        relpath, _, cols = str(raw).partition(":")
        if not relpath or not cols:
            raise ValueError(f"Invalid column mapping: {raw}")
        mapping[relpath.replace("\\", "/")] = [item.strip() for item in cols.split(",") if item.strip()]
    return mapping


def detect_key_columns(df: pd.DataFrame, override):
    if override:
        return [column for column in override if column in df.columns]
    for candidate in (
        ["sample_id", "chrom", "start", "end", "state"],
        ["sample_id", "chrom", "start", "end", "expected_state"],
        ["sample_id"],
    ):
        if all(column in df.columns for column in candidate):
            return candidate
    return []


def compare_tables(path_a: Path, path_b: Path, float_tol: float, key_override, ignored_columns):
    import pandas as pd

    sep = "\t" if path_a.suffix.lower() == ".tsv" else ","
    df_a = pd.read_csv(path_a, sep=sep)
    df_b = pd.read_csv(path_b, sep=sep)
    result = {
        "kind": "table",
        "columns_equal": list(df_a.columns) == list(df_b.columns),
        "row_count_equal": int(len(df_a)) == int(len(df_b)),
        "baseline_rows": int(len(df_a)),
        "candidate_rows": int(len(df_b)),
    }
    if list(df_a.columns) != list(df_b.columns):
        result["status"] = "schema_changed"
        result["baseline_columns"] = list(df_a.columns)
        result["candidate_columns"] = list(df_b.columns)
        return result

    key_columns = detect_key_columns(df_a, key_override)
    if key_columns:
        df_a = df_a.sort_values(key_columns).reset_index(drop=True)
        df_b = df_b.sort_values(key_columns).reset_index(drop=True)
        result["key_columns"] = key_columns

    if len(df_a) != len(df_b):
        result["status"] = "row_count_changed"
        return result

    ignored = [column for column in ignored_columns if column in df_a.columns]
    if ignored:
        result["ignored_columns"] = ignored

    float_diffs = []
    nonfloat_changes = []
    for column in df_a.columns:
        if column in ignored:
            continue
        series_a = df_a[column]
        series_b = df_b[column]
        numeric_a = pd.to_numeric(series_a, errors="coerce")
        numeric_b = pd.to_numeric(series_b, errors="coerce")
        if numeric_a.notna().all() and numeric_b.notna().all():
            # Normalize bool-like numerics to float before subtraction.
            numeric_a = numeric_a.astype(float)
            numeric_b = numeric_b.astype(float)
            max_diff = float((numeric_a - numeric_b).abs().max()) if len(df_a) else 0.0
            if max_diff > float_tol:
                float_diffs.append({"column": column, "max_abs_diff": max_diff})
            continue
        equal = series_a.fillna("").astype(str).equals(series_b.fillna("").astype(str))
        if not equal:
            nonfloat_changes.append(column)

    if nonfloat_changes:
        result["status"] = "value_changed"
        result["changed_columns"] = nonfloat_changes
    elif float_diffs:
        result["status"] = "float_tolerance_only"
        result["float_differences"] = float_diffs
    else:
        result["status"] = "identical"
    return result


def compare_json(path_a: Path, path_b: Path, float_tol: float):
    payload_a = json.loads(path_a.read_text(encoding="utf-8"))
    payload_b = json.loads(path_b.read_text(encoding="utf-8"))
    result = {"kind": "json"}

    def normalize(payload):
        if isinstance(payload, dict):
            return {key: normalize(value) for key, value in sorted(payload.items())}
        if isinstance(payload, list):
            return [normalize(value) for value in payload]
        if isinstance(payload, float) and math.isfinite(payload):
            return round(payload, 12)
        return payload

    norm_a = normalize(payload_a)
    norm_b = normalize(payload_b)
    if norm_a == norm_b:
        result["status"] = "identical"
        return result
    keys_a = sorted(payload_a.keys()) if isinstance(payload_a, dict) else []
    keys_b = sorted(payload_b.keys()) if isinstance(payload_b, dict) else []
    result["status"] = "value_changed"
    result["baseline_keys"] = keys_a
    result["candidate_keys"] = keys_b
    result["float_tolerance"] = float_tol
    return result


def compare_text(path_a: Path, path_b: Path):
    text_a = path_a.read_text(encoding="utf-8")
    text_b = path_b.read_text(encoding="utf-8")
    return {
        "kind": "text",
        "status": "identical" if text_a == text_b else "value_changed",
        "baseline_length": len(text_a),
        "candidate_length": len(text_b),
    }


def compare_file(path_a: Path, path_b: Path, float_tol: float, key_override, ignored_columns):
    suffix = path_a.suffix.lower()
    if suffix in TABLE_SUFFIXES:
        return compare_tables(path_a, path_b, float_tol, key_override, ignored_columns)
    if suffix in JSON_SUFFIXES:
        return compare_json(path_a, path_b, float_tol)
    if suffix in TEXT_SUFFIXES:
        return compare_text(path_a, path_b)
    return {
        "kind": "binary_or_other",
        "status": "identical" if path_a.read_bytes() == path_b.read_bytes() else "value_changed",
    }


def build_report(results):
    lines = [
        "# Output Comparison Report",
        "",
        "| Relative path | Status | Kind |",
        "| --- | --- | --- |",
    ]
    for relpath, item in sorted(results.items()):
        lines.append(f"| `{relpath}` | `{item['status']}` | `{item['kind']}` |")
    return "\n".join(lines) + "\n"


def main():
    args = parse_args()
    baseline_root = Path(args.baseline_root)
    candidate_root = Path(args.candidate_root)
    overrides = parse_key_overrides(args.key_column)
    ignored_column_map = parse_column_mapping(args.ignore_column)

    relpaths = sorted(set(args.relpath)) if args.relpath else sorted(set(discover_paths(baseline_root)) | set(discover_paths(candidate_root)))
    results = {}
    for relpath in relpaths:
        relpath_norm = relpath.replace("\\", "/")
        path_a = baseline_root / relpath_norm
        path_b = candidate_root / relpath_norm
        if not path_a.exists() or not path_b.exists():
            results[relpath_norm] = {
                "kind": "missing",
                "status": "missing",
                "baseline_exists": path_a.exists(),
                "candidate_exists": path_b.exists(),
            }
            continue
        results[relpath_norm] = compare_file(
            path_a,
            path_b,
            args.float_tol,
            overrides.get(relpath_norm, []),
            ignored_column_map.get(relpath_norm, []),
        )

    summary = {
        "baseline_root": str(baseline_root),
        "candidate_root": str(candidate_root),
        "float_tol": args.float_tol,
        "results": results,
    }
    ensure_parent(args.output_json).write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    ensure_parent(args.output_md).write_text(build_report(results), encoding="utf-8")


if __name__ == "__main__":
    main()
