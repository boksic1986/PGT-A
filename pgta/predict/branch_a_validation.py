#!/biosoftware/miniconda/envs/snakemake_env/bin/python
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import pandas as pd

from pgta.predict.truth import (
    filter_truth_to_samples,
    normalize_chrom,
    normalize_state,
    overlap_fraction,
)


OUTPUT_COLUMNS = [
    "reference_id",
    "binsize",
    "preprocess_mask_version",
    "wisecondorx_predict_command",
    "blacklist_status",
    "minrefbins",
    "zscore",
    "alpha",
    "maskrepeats",
    "sample",
    "qc_status",
    "sex_call",
    "truth_event_count",
    "truth_detected_count",
    "FN_count",
    "branch_a_candidate_count",
    "branch_a_strong_count",
    "branch_a_sensitive_count",
    "top_branch_a_signal",
    "top_branch_a_abs_z",
    "H6_chr21_status",
]


TRUTH_COLUMNS = [
    "sample",
    "truth_index",
    "chrom",
    "start",
    "end",
    "state",
    "detected",
    "matched_candidate_id",
    "matched_a_abs_zscore",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarize Branch A no-FN validation under a fixed reference contract."
    )
    parser.add_argument("--candidate-tsv", action="append", default=[])
    parser.add_argument("--truth-tsv", default="")
    parser.add_argument("--qc-tsv", action="append", default=[])
    parser.add_argument("--gender-tsv", action="append", default=[])
    parser.add_argument("--sample-id", action="append", default=[])
    parser.add_argument("--reference-id", required=True)
    parser.add_argument("--binsize", required=True)
    parser.add_argument("--preprocess-mask-version", default="")
    parser.add_argument("--wisecondorx-predict-command", default="")
    parser.add_argument("--blacklist-bed", default="")
    parser.add_argument("--minrefbins", default="")
    parser.add_argument("--zscore", default="")
    parser.add_argument("--alpha", default="")
    parser.add_argument("--maskrepeats", default="")
    parser.add_argument("--output-sample-summary", required=True)
    parser.add_argument("--output-truth-metrics", required=True)
    parser.add_argument("--output-summary", required=True)
    return parser.parse_args()


def _read_table(path_value: str, columns: list[str] | None = None) -> pd.DataFrame:
    path = Path(path_value)
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame(columns=columns or [])
    return pd.read_csv(path, sep="\t")


def _load_candidates(candidate_paths: list[str]) -> pd.DataFrame:
    frames = []
    for path in candidate_paths:
        df = _read_table(path)
        if df.empty:
            sample_id = _sample_from_candidate_path(path)
            frames.append(pd.DataFrame([{"sample_id": sample_id}]).iloc[0:0])
            continue
        frames.append(df)
    if not frames:
        return pd.DataFrame(columns=["sample_id", "chrom", "start", "end", "state"])
    out = pd.concat(frames, ignore_index=True, sort=False)
    if "sample" in out.columns and "sample_id" not in out.columns:
        out = out.rename(columns={"sample": "sample_id"})
    if "sample_id" not in out.columns:
        out["sample_id"] = ""
    for column in ["chrom", "state"]:
        if column not in out.columns:
            out[column] = ""
    for column in ["start", "end", "a_abs_zscore", "a_zscore"]:
        if column not in out.columns:
            out[column] = math.nan
    out["sample_id"] = out["sample_id"].astype(str)
    out["chrom"] = out["chrom"].map(normalize_chrom)
    out["state"] = out["state"].map(normalize_state)
    out["start"] = pd.to_numeric(out["start"], errors="coerce")
    out["end"] = pd.to_numeric(out["end"], errors="coerce")
    out["a_abs_zscore"] = pd.to_numeric(out["a_abs_zscore"], errors="coerce")
    out["a_zscore"] = pd.to_numeric(out["a_zscore"], errors="coerce")
    return out


def _sample_from_candidate_path(path_value: str) -> str:
    name = Path(path_value).name
    suffix = ".candidate_events.tsv"
    return name[: -len(suffix)] if name.endswith(suffix) else Path(path_value).stem


def _load_truth(truth_tsv: str, sample_ids: list[str]) -> pd.DataFrame:
    if not truth_tsv:
        return pd.DataFrame(columns=["sample_id", "chrom", "start", "end", "state"])
    df = _read_table(truth_tsv)
    if df.empty:
        return pd.DataFrame(columns=["sample_id", "chrom", "start", "end", "state"])
    if "sample" in df.columns and "sample_id" not in df.columns:
        df = df.rename(columns={"sample": "sample_id"})
    if "expected_state" in df.columns and "state" not in df.columns:
        df = df.rename(columns={"expected_state": "state"})
    if "svtype" in df.columns and "state" not in df.columns:
        df = df.rename(columns={"svtype": "state"})
    required = {"sample_id", "chrom", "state"}
    missing = sorted(required.difference(df.columns))
    if missing:
        raise ValueError(f"truth table missing columns: {','.join(missing)}")
    for column in ["start", "end"]:
        if column not in df.columns:
            df[column] = math.nan
    df = filter_truth_to_samples(df, sample_ids)
    df["sample_id"] = df["sample_id"].astype(str)
    df["chrom"] = df["chrom"].map(normalize_chrom)
    df["state"] = df["state"].map(normalize_state)
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    return df.reset_index(drop=True)


def _load_status(paths: list[str], sample_ids: list[str], kind: str) -> dict[str, str]:
    status = {sample_id: "" for sample_id in sample_ids}
    for path_value in paths:
        df = _read_table(path_value)
        if df.empty:
            continue
        if "sample" in df.columns and "sample_id" not in df.columns:
            df = df.rename(columns={"sample": "sample_id"})
        sample_column = "sample_id" if "sample_id" in df.columns else ""
        if not sample_column:
            sample_id = _sample_from_generic_path(path_value)
            rows = [df.iloc[0]]
        else:
            rows = [row for _, row in df.iterrows()]
        for row in rows:
            sample_id = str(row.get(sample_column, _sample_from_generic_path(path_value)))
            if kind == "qc":
                value = _first_present(row, ["qc_status", "status", "decision", "pass_fail"], "")
            else:
                value = _first_present(row, ["sex_call", "gender", "sex", "prediction"], "")
            status[sample_id] = str(value)
    return status


def _sample_from_generic_path(path_value: str) -> str:
    name = Path(path_value).name
    for suffix in [".qc.tsv", ".gender.tsv", ".tsv"]:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(path_value).stem


def _first_present(row, columns: list[str], default: str):
    for column in columns:
        if column in row and not pd.isna(row[column]):
            return row[column]
    return default


def _candidate_matches_truth(candidates: pd.DataFrame, truth_row) -> pd.DataFrame:
    if candidates.empty:
        return candidates
    sample_id = str(truth_row.sample_id)
    chrom = normalize_chrom(truth_row.chrom)
    state = normalize_state(truth_row.state)
    subset = candidates[
        candidates["sample_id"].astype(str).eq(sample_id)
        & candidates["chrom"].astype(str).eq(chrom)
        & candidates["state"].astype(str).eq(state)
    ].copy()
    if subset.empty:
        return subset
    if pd.isna(truth_row.start) or pd.isna(truth_row.end):
        return subset
    overlaps = subset.apply(
        lambda row: overlap_fraction(row["start"], row["end"], truth_row.start, truth_row.end),
        axis=1,
    )
    return subset[overlaps.fillna(0.0) > 0.0].copy()


def _top_candidate_signal(candidates: pd.DataFrame) -> tuple[str, float]:
    if candidates.empty:
        return "", math.nan
    ranked = candidates.copy()
    ranked["_rank_z"] = pd.to_numeric(ranked["a_abs_zscore"], errors="coerce").fillna(-1.0)
    row = ranked.sort_values("_rank_z", ascending=False).iloc[0]
    abs_z = float(row["_rank_z"]) if row["_rank_z"] >= 0 else math.nan
    signal = (
        f"{row.get('chrom', '')}:{int(row.get('start', 0))}-{int(row.get('end', 0))} "
        f"{row.get('state', '')} z={row.get('a_zscore', '')}"
    )
    return signal, abs_z


def _blacklist_status(blacklist_bed: str) -> str:
    if not str(blacklist_bed).strip():
        return "not_configured"
    return "configured"


def build_branch_a_validation(
    candidate_paths: list[str],
    truth_tsv: str,
    sample_ids: list[str],
    reference_id: str,
    binsize,
    preprocess_mask_version: str,
    wisecondorx_predict_command: str,
    blacklist_bed: str,
    minrefbins,
    zscore,
    alpha,
    maskrepeats,
    output_sample_summary: str,
    output_truth_metrics: str,
    output_summary: str,
    qc_paths: list[str] | None = None,
    gender_paths: list[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    inferred_samples = [_sample_from_candidate_path(path) for path in candidate_paths]
    ordered_samples = list(dict.fromkeys([str(s) for s in sample_ids if str(s)] + inferred_samples))
    candidates = _load_candidates(candidate_paths)
    truth = _load_truth(truth_tsv, ordered_samples)
    qc_status = _load_status(qc_paths or [], ordered_samples, "qc")
    sex_call = _load_status(gender_paths or [], ordered_samples, "gender")

    truth_rows = []
    for truth_index, row in truth.iterrows():
        matches = _candidate_matches_truth(candidates, row)
        if matches.empty:
            truth_rows.append(
                {
                    "sample": str(row.sample_id),
                    "truth_index": int(truth_index),
                    "chrom": row.chrom,
                    "start": row.start,
                    "end": row.end,
                    "state": row.state,
                    "detected": 0,
                    "matched_candidate_id": "",
                    "matched_a_abs_zscore": math.nan,
                }
            )
            continue
        match = matches.sort_values("a_abs_zscore", ascending=False).iloc[0]
        truth_rows.append(
            {
                "sample": str(row.sample_id),
                "truth_index": int(truth_index),
                "chrom": row.chrom,
                "start": row.start,
                "end": row.end,
                "state": row.state,
                "detected": 1,
                "matched_candidate_id": str(match.get("candidate_id", "")),
                "matched_a_abs_zscore": float(match.get("a_abs_zscore", math.nan)),
            }
        )
    truth_metrics = pd.DataFrame(truth_rows, columns=TRUTH_COLUMNS)

    sample_rows = []
    for sample_id in ordered_samples:
        sample_candidates = candidates[candidates["sample_id"].astype(str).eq(sample_id)].copy()
        sample_truth = truth[truth["sample_id"].astype(str).eq(sample_id)].copy()
        sample_truth_metrics = truth_metrics[truth_metrics["sample"].astype(str).eq(sample_id)].copy()
        truth_event_count = int(len(sample_truth))
        truth_detected_count = int(sample_truth_metrics["detected"].sum()) if not sample_truth_metrics.empty else 0
        fn_count = max(truth_event_count - truth_detected_count, 0)
        top_signal, top_abs_z = _top_candidate_signal(sample_candidates)
        strong_count = (
            int(sample_candidates["a_support_level"].astype(str).eq("strong").sum())
            if "a_support_level" in sample_candidates.columns
            else int(sample_candidates["a_abs_zscore"].ge(float(zscore or 0)).sum())
        )
        sensitive_count = int(len(sample_candidates) - strong_count)
        h6_chr21_status = "not_applicable"
        if sample_id == "H6":
            h6_truth = sample_truth[
                sample_truth["chrom"].astype(str).eq("chr21")
                & sample_truth["state"].astype(str).eq("gain")
            ]
            if h6_truth.empty:
                h6_chr21_status = "truth_not_available"
            else:
                h6_detected = sample_truth_metrics[
                    sample_truth_metrics["chrom"].astype(str).eq("chr21")
                    & sample_truth_metrics["state"].astype(str).eq("gain")
                    & sample_truth_metrics["detected"].astype(int).eq(1)
                ]
                h6_chr21_status = "detected" if not h6_detected.empty else "not_detected"
        sample_rows.append(
            {
                "reference_id": reference_id,
                "binsize": binsize,
                "preprocess_mask_version": preprocess_mask_version,
                "wisecondorx_predict_command": wisecondorx_predict_command,
                "blacklist_status": _blacklist_status(blacklist_bed),
                "minrefbins": minrefbins,
                "zscore": zscore,
                "alpha": alpha,
                "maskrepeats": maskrepeats,
                "sample": sample_id,
                "qc_status": qc_status.get(sample_id, ""),
                "sex_call": sex_call.get(sample_id, ""),
                "truth_event_count": truth_event_count,
                "truth_detected_count": truth_detected_count,
                "FN_count": fn_count,
                "branch_a_candidate_count": int(len(sample_candidates)),
                "branch_a_strong_count": strong_count,
                "branch_a_sensitive_count": sensitive_count,
                "top_branch_a_signal": top_signal,
                "top_branch_a_abs_z": top_abs_z,
                "H6_chr21_status": h6_chr21_status,
            }
        )
    sample_summary = pd.DataFrame(sample_rows, columns=OUTPUT_COLUMNS)
    h6_rows = sample_summary[sample_summary["sample"].astype(str).eq("H6")]
    h6_status = str(h6_rows["H6_chr21_status"].iloc[0]) if not h6_rows.empty else "not_applicable"
    summary = {
        "status": "ready",
        "reference_id": reference_id,
        "sample_count": int(len(sample_summary)),
        "truth_event_count": int(sample_summary["truth_event_count"].sum()) if not sample_summary.empty else 0,
        "truth_detected_count": int(sample_summary["truth_detected_count"].sum()) if not sample_summary.empty else 0,
        "FN_count": int(sample_summary["FN_count"].sum()) if not sample_summary.empty else 0,
        "branch_a_candidate_count": int(sample_summary["branch_a_candidate_count"].sum()) if not sample_summary.empty else 0,
        "H6_chr21_status": h6_status,
        "branch_b_used_for_pass_fail": False,
        "branch_s_used_for_pass_fail": False,
        "report_outputs_used_for_pass_fail": False,
    }

    for path_value, df in [
        (output_sample_summary, sample_summary),
        (output_truth_metrics, truth_metrics),
    ]:
        path = Path(path_value)
        path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(path, sep="\t", index=False)
    summary_path = Path(output_summary)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return sample_summary, truth_metrics, summary


def main():
    args = parse_args()
    build_branch_a_validation(
        candidate_paths=args.candidate_tsv,
        truth_tsv=args.truth_tsv,
        qc_paths=args.qc_tsv,
        gender_paths=args.gender_tsv,
        sample_ids=args.sample_id,
        reference_id=args.reference_id,
        binsize=args.binsize,
        preprocess_mask_version=args.preprocess_mask_version,
        wisecondorx_predict_command=args.wisecondorx_predict_command,
        blacklist_bed=args.blacklist_bed,
        minrefbins=args.minrefbins,
        zscore=args.zscore,
        alpha=args.alpha,
        maskrepeats=args.maskrepeats,
        output_sample_summary=args.output_sample_summary,
        output_truth_metrics=args.output_truth_metrics,
        output_summary=args.output_summary,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
