#!/usr/bin/env python
import argparse
import csv
import json
import sys
from pathlib import Path


REQUIRED_COLUMNS = [
    "sample_id",
    "chrom",
    "expected_state",
    "abnormal_cell_fraction",
]

OPTIONAL_KEY_COLUMNS = ["start", "end"]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Validate fraction truth TSV before using it for CNV benchmark/evaluation."
    )
    parser.add_argument("--input-tsv", required=True, help="Truth TSV with abnormal_cell_fraction column.")
    parser.add_argument("--output-json", default="", help="Optional validation summary JSON path.")
    parser.add_argument(
        "--require-complete",
        action="store_true",
        help="Exit non-zero when any abnormal_cell_fraction is missing or invalid.",
    )
    return parser.parse_args()


def normalize_text(value):
    return str(value or "").strip()


def build_truth_key(row):
    key_fields = [
        normalize_text(row.get("sample_id")),
        normalize_text(row.get("chrom")),
        normalize_text(row.get("expected_state")).lower(),
        normalize_text(row.get("start")),
        normalize_text(row.get("end")),
    ]
    return tuple(key_fields)


def parse_fraction(raw_value):
    text = normalize_text(raw_value)
    if not text:
        return None, "missing"
    try:
        value = float(text)
    except ValueError:
        return None, "not_numeric"
    if not (0.0 <= value <= 1.0):
        return None, "out_of_range"
    return value, ""


def validate_truth_table(input_tsv):
    path = Path(input_tsv)
    if not path.exists():
        raise FileNotFoundError(f"Truth TSV not found: {path}")

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames or []
        rows = list(reader)

    missing_columns = [column for column in REQUIRED_COLUMNS if column not in header]
    summary = {
        "status": "ok",
        "input_tsv": str(path),
        "row_count": len(rows),
        "required_columns": REQUIRED_COLUMNS,
        "missing_required_columns": missing_columns,
        "fraction_column": "abnormal_cell_fraction",
        "complete_fraction_row_count": 0,
        "missing_fraction_row_count": 0,
        "invalid_fraction_row_count": 0,
        "duplicate_truth_key_count": 0,
        "missing_fraction_rows": [],
        "invalid_fraction_rows": [],
        "duplicate_truth_keys": [],
    }
    if missing_columns:
        summary["status"] = "invalid_schema"
        return summary

    seen_keys = {}
    for index, row in enumerate(rows, start=2):
        truth_key = build_truth_key(row)
        seen_keys.setdefault(truth_key, []).append(index)

        _, error = parse_fraction(row.get("abnormal_cell_fraction"))
        row_label = {
            "line": index,
            "sample_id": normalize_text(row.get("sample_id")),
            "chrom": normalize_text(row.get("chrom")),
            "expected_state": normalize_text(row.get("expected_state")),
            "start": normalize_text(row.get("start")),
            "end": normalize_text(row.get("end")),
        }
        if error == "missing":
            summary["missing_fraction_row_count"] += 1
            summary["missing_fraction_rows"].append(row_label)
        elif error:
            summary["invalid_fraction_row_count"] += 1
            row_label["error"] = error
            row_label["raw_fraction"] = normalize_text(row.get("abnormal_cell_fraction"))
            summary["invalid_fraction_rows"].append(row_label)
        else:
            summary["complete_fraction_row_count"] += 1

    duplicate_keys = []
    for key_fields, lines in seen_keys.items():
        if len(lines) <= 1:
            continue
        duplicate_keys.append(
            {
                "truth_key": {
                    "sample_id": key_fields[0],
                    "chrom": key_fields[1],
                    "expected_state": key_fields[2],
                    "start": key_fields[3],
                    "end": key_fields[4],
                },
                "lines": lines,
            }
        )
    summary["duplicate_truth_key_count"] = len(duplicate_keys)
    summary["duplicate_truth_keys"] = duplicate_keys

    if summary["missing_fraction_row_count"] or summary["invalid_fraction_row_count"] or missing_columns:
        summary["status"] = "incomplete"
    if summary["duplicate_truth_key_count"]:
        summary["status"] = "needs_review" if summary["status"] == "ok" else summary["status"]
    return summary


def main():
    args = parse_args()
    summary = validate_truth_table(args.input_tsv)

    if args.output_json:
        output_path = Path(args.output_json)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

    print(json.dumps(summary, indent=2, ensure_ascii=False))

    if args.require_complete and summary["status"] != "ok":
        sys.exit(1)


if __name__ == "__main__":
    main()
