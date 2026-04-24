#!/biosoftware/miniconda/envs/snakemake_env/bin/python
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from pgta.core.logging import setup_logger


def parse_args():
    parser = argparse.ArgumentParser(description="Evaluate CNV candidate events and sample-level technical metrics.")
    parser.add_argument("--event-tsv", action="append", default=[])
    parser.add_argument("--gender-tsv", action="append", default=[])
    parser.add_argument("--qc-tsv", action="append", default=[])
    parser.add_argument("--truth-tsv", default="")
    parser.add_argument("--branch-b-z-threshold", type=float, default=2.5)
    parser.add_argument("--output-sample-metrics", required=True)
    parser.add_argument("--output-event-metrics", required=True)
    parser.add_argument("--output-calibration", required=True)
    parser.add_argument("--output-summary", required=True)
    parser.add_argument("--log", default="")
    return parser.parse_args()


def ensure_parent(path_value):
    path = Path(path_value)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def load_tables(paths):
    frames = []
    for path_value in paths:
        path = Path(path_value)
        if not path.exists():
            continue
        df = pd.read_csv(path, sep="\t")
        if not df.empty:
            frames.append(df)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def sample_id_from_path(path_value):
    return Path(path_value).stem.split(".")[0]


def normalize_chrom(value):
    chrom = str(value).strip()
    if not chrom:
        return chrom
    return chrom if chrom.startswith("chr") else f"chr{chrom}"


def normalize_state(value):
    text = str(value).strip().lower()
    mapping = {
        "gain": "gain",
        "dup": "gain",
        "duplication": "gain",
        "amplification": "gain",
        "+": "gain",
        "loss": "loss",
        "del": "loss",
        "deletion": "loss",
        "-": "loss",
    }
    return mapping.get(text, text)


def load_gender_tables(paths):
    rows = []
    for path_value in paths:
        path = Path(path_value)
        if not path.exists():
            continue
        df = pd.read_csv(path, sep="\t")
        if df.empty:
            continue
        row = df.iloc[0].to_dict()
        row.setdefault("sample_id", sample_id_from_path(path_value))
        rows.append(row)
    return pd.DataFrame(rows)


def load_qc_tables(paths):
    rows = []
    for path_value in paths:
        path = Path(path_value)
        if not path.exists():
            continue
        df = pd.read_csv(path, sep="\t")
        if df.empty:
            continue
        row = df.iloc[0].to_dict()
        row.setdefault("sample_id", sample_id_from_path(path_value))
        rows.append(row)
    return pd.DataFrame(rows)


def overlap_fraction(left_start, left_end, right_start, right_end):
    if any(pd.isna(value) for value in [left_start, left_end, right_start, right_end]):
        return np.nan
    overlap = max(0.0, min(float(left_end), float(right_end)) - max(float(left_start), float(right_start)) + 1.0)
    length = max(float(left_end) - float(left_start) + 1.0, 1.0)
    return overlap / length


def event_support_z(frame):
    support_columns = [
        column for column in ["calibrated_mean_z", "calibrated_median_z", "event_corr_adjusted_z"] if column in frame.columns
    ]
    if not support_columns or frame.empty:
        return pd.Series(np.nan, index=frame.index, dtype=float)
    support_frame = frame[support_columns].apply(pd.to_numeric, errors="coerce").abs()
    return support_frame.max(axis=1, skipna=True)


def compute_truth_metrics(events_df, truth_tsv, branch_b_z_threshold):
    path = Path(truth_tsv)
    if not truth_tsv or not path.exists():
        return pd.DataFrame(), {}
    truth_df = pd.read_csv(path, sep="\t")
    if truth_df.empty or not {"sample_id", "chrom", "expected_state"}.issubset(truth_df.columns):
        return pd.DataFrame(), {}
    truth_df = truth_df.copy()
    truth_df["chrom"] = truth_df["chrom"].map(normalize_chrom)
    truth_df["expected_state"] = truth_df["expected_state"].map(normalize_state)
    for column in ["start", "end"]:
        if column not in truth_df.columns:
            truth_df[column] = np.nan
    kept = events_df[events_df["keep_event"] == 1].copy()
    if kept.empty:
        metrics_df = truth_df[["sample_id", "chrom", "expected_state"]].copy()
        metrics_df["matched"] = 0
        metrics_df["detected"] = 0
        metrics_df["top_support_z"] = np.nan
        return metrics_df, {
            "truth_row_count": int(len(metrics_df)),
            "truth_match_count": 0,
            "truth_detected_count": 0,
            "truth_match_rate": 0.0 if len(metrics_df) else None,
            "truth_recall": 0.0 if len(metrics_df) else None,
        }
    kept["chrom"] = kept["chrom"].map(normalize_chrom)
    kept["state"] = kept["state"].map(normalize_state)
    kept["support_z"] = event_support_z(kept)
    rows = []
    for truth in truth_df.itertuples(index=False):
        sample_events = kept[
            kept["sample_id"].astype(str).eq(str(truth.sample_id))
            & kept["chrom"].astype(str).eq(str(truth.chrom))
            & kept["state"].astype(str).eq(str(truth.expected_state))
        ].copy()
        if not sample_events.empty and pd.notna(truth.start) and pd.notna(truth.end) and {"start", "end"}.issubset(sample_events.columns):
            sample_events["truth_overlap_fraction"] = sample_events.apply(
                lambda row: overlap_fraction(truth.start, truth.end, row.start, row.end),
                axis=1,
            )
            sample_events = sample_events[sample_events["truth_overlap_fraction"].fillna(0.0) > 0.0].copy()
        top_support = float(sample_events["support_z"].max()) if not sample_events.empty else np.nan
        rows.append(
            {
                "sample_id": truth.sample_id,
                "chrom": truth.chrom,
                "expected_state": truth.expected_state,
                "matched": int(not sample_events.empty),
                "detected": int(np.isfinite(top_support) and top_support >= float(branch_b_z_threshold)),
                "top_support_z": top_support,
            }
        )
    metrics_df = pd.DataFrame(rows)
    summary = {
        "truth_row_count": int(len(metrics_df)),
        "truth_match_count": int(metrics_df["matched"].sum()) if not metrics_df.empty else 0,
        "truth_detected_count": int(metrics_df["detected"].sum()) if not metrics_df.empty else 0,
        "truth_match_rate": float(metrics_df["matched"].mean()) if not metrics_df.empty else None,
        "truth_recall": float(metrics_df["detected"].mean()) if not metrics_df.empty else None,
    }
    return metrics_df, summary


def main():
    args = parse_args()
    logger = setup_logger("cnv_evaluation", args.log or None)
    events_df = load_tables(args.event_tsv)
    genders_df = load_gender_tables(args.gender_tsv)
    qc_df = load_qc_tables(args.qc_tsv)

    if events_df.empty:
        empty = pd.DataFrame(columns=["sample_id"])
        ensure_parent(args.output_sample_metrics)
        empty.to_csv(args.output_sample_metrics, sep="\t", index=False)
        empty.to_csv(args.output_event_metrics, sep="\t", index=False)
        empty.to_csv(args.output_calibration, sep="\t", index=False)
        ensure_parent(args.output_summary).write_text(json.dumps({"status": "empty"}, indent=2), encoding="utf-8")
        logger.info("evaluation skipped: no events")
        return

    event_metrics = (
        events_df.groupby(["chrom", "state", "artifact_status"], dropna=False)
        .agg(
            event_count=("event_id", "size"),
            kept_event_count=("keep_event", "sum"),
            median_priority_score=("priority_score", "median"),
            median_abs_calibrated_z=("calibrated_mean_z", lambda values: float(np.nanmedian(np.abs(values)))),
        )
        .reset_index()
    )

    sample_metrics = (
        events_df.groupby("sample_id", dropna=False)
        .agg(
            total_events=("event_id", "size"),
            kept_events=("keep_event", "sum"),
            pass_events=("artifact_status", lambda values: int((values == "pass").sum())),
            review_events=("artifact_status", lambda values: int((values == "review").sum())),
            top_priority_score=("priority_score", "max"),
            max_abs_calibrated_z=("calibrated_mean_z", lambda values: float(np.nanmax(np.abs(values)))),
        )
        .reset_index()
    )
    if "high_risk_bin_fraction" in events_df.columns:
        clean_fraction = (
            pd.to_numeric(events_df["clean_bin_fraction"], errors="coerce").fillna(0.0)
            if "clean_bin_fraction" in events_df.columns
            else pd.Series(0.0, index=events_df.index)
        )
        risk_metrics = (
            events_df.assign(
                high_risk_event=(events_df["high_risk_bin_fraction"].fillna(0.0) >= 0.50).astype(int),
                clean_supported_event=(clean_fraction >= 0.50).astype(int),
            )
            .groupby("sample_id", dropna=False)
            .agg(
                high_risk_events=("high_risk_event", "sum"),
                clean_supported_events=("clean_supported_event", "sum"),
            )
            .reset_index()
        )
        sample_metrics = sample_metrics.merge(risk_metrics, on="sample_id", how="left")
    if not genders_df.empty:
        sample_metrics = sample_metrics.merge(
            genders_df[[column for column in ["sample_id", "sex_call", "sex_call_source", "predict_gender"] if column in genders_df.columns]],
            on="sample_id",
            how="left",
        )
    if not qc_df.empty:
        sample_metrics = sample_metrics.merge(
            qc_df[[column for column in ["sample_id", "status", "mad_log1p", "nonzero_fraction"] if column in qc_df.columns]].rename(
                columns={"status": "qc_status"}
            ),
            on="sample_id",
            how="left",
        )

    calibration_df = (
        events_df.assign(
            qvalue_bin=pd.cut(
                events_df["empirical_qvalue"].fillna(1.0),
                bins=[-0.01, 0.01, 0.05, 0.10, 0.25, 1.0],
                labels=["<=0.01", "0.01-0.05", "0.05-0.10", "0.10-0.25", ">0.25"],
            )
        )
        .groupby("qvalue_bin", dropna=False)
        .agg(
            event_count=("event_id", "size"),
            kept_event_count=("keep_event", "sum"),
            median_priority_score=("priority_score", "median"),
        )
        .reset_index()
    )

    truth_metrics_df, truth_summary = compute_truth_metrics(events_df, args.truth_tsv, args.branch_b_z_threshold)
    if not truth_metrics_df.empty:
        truth_metrics_df.to_csv(str(Path(args.output_sample_metrics).with_name("truth_metrics.tsv")), sep="\t", index=False)

    ensure_parent(args.output_sample_metrics)
    sample_metrics.to_csv(args.output_sample_metrics, sep="\t", index=False)
    event_metrics.to_csv(args.output_event_metrics, sep="\t", index=False)
    calibration_df.to_csv(args.output_calibration, sep="\t", index=False)
    summary = {
        "status": "completed",
        "sample_count": int(sample_metrics["sample_id"].nunique()),
        "event_count": int(len(events_df)),
        "kept_event_count": int(events_df["keep_event"].sum()),
        "pass_event_count": int((events_df["artifact_status"] == "pass").sum()),
        "review_event_count": int((events_df["artifact_status"] == "review").sum()),
        "artifact_event_count": int((events_df["artifact_status"] == "artifact").sum()),
    }
    if "high_risk_bin_fraction" in events_df.columns:
        summary["high_risk_event_count"] = int((events_df["high_risk_bin_fraction"].fillna(0.0) >= 0.50).sum())
        summary["clean_supported_event_count"] = int((events_df["clean_bin_fraction"].fillna(0.0) >= 0.50).sum())
    summary.update(truth_summary)
    ensure_parent(args.output_summary).write_text(json.dumps(summary, indent=2), encoding="utf-8")
    logger.info("evaluation completed: samples=%d events=%d", summary["sample_count"], summary["event_count"])


if __name__ == "__main__":
    main()
