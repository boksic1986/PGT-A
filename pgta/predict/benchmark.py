#!/biosoftware/miniconda/envs/snakemake_env/bin/python
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from pgta.core.logging import setup_logger


def parse_args():
    parser = argparse.ArgumentParser(description="Run truth-based simulation and admixture benchmarking for CNV candidates.")
    parser.add_argument("--event-tsv", action="append", default=[])
    parser.add_argument("--a-branch-bed", action="append", default=[])
    parser.add_argument("--truth-tsv", default="")
    parser.add_argument("--admixture-level", action="append", type=float, default=[])
    parser.add_argument("--low-fraction-threshold", action="append", type=float, default=[])
    parser.add_argument("--branch-b-z-threshold", type=float, default=2.5)
    parser.add_argument("--branch-a-z-threshold", type=float, default=5.0)
    parser.add_argument("--output-simulation", required=True)
    parser.add_argument("--output-admixture", required=True)
    parser.add_argument("--output-summary", required=True)
    parser.add_argument("--log", default="")
    return parser.parse_args()


def ensure_parent(path_value):
    path = Path(path_value)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


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


def safe_float(value, default=np.nan):
    try:
        return float(value)
    except (TypeError, ValueError):
        return float(default)


def load_truth_table(path_value):
    path = Path(path_value)
    if not path_value or not path.exists():
        return pd.DataFrame()
    truth_df = pd.read_csv(path, sep="\t")
    if truth_df.empty:
        return truth_df
    truth_df = truth_df.copy()
    if "chrom" in truth_df.columns:
        truth_df["chrom"] = truth_df["chrom"].map(normalize_chrom)
    if "expected_state" in truth_df.columns:
        truth_df["expected_state"] = truth_df["expected_state"].map(normalize_state)
    fraction_column_candidates = [
        "abnormal_cell_fraction",
        "expected_abnormal_cell_fraction",
        "expected_fraction",
        "truth_fraction",
        "mosaic_fraction",
        "abnormal_fraction",
    ]
    truth_fraction_column = next((column for column in fraction_column_candidates if column in truth_df.columns), "")
    truth_df["truth_fraction"] = (
        pd.to_numeric(truth_df[truth_fraction_column], errors="coerce")
        if truth_fraction_column
        else np.nan
    )
    truth_df["truth_fraction_source"] = truth_fraction_column
    for column in ["start", "end"]:
        if column not in truth_df.columns:
            truth_df[column] = np.nan
    return truth_df


def load_branch_b_events(paths):
    frames = []
    for path_value in paths:
        path = Path(path_value)
        if not path.exists():
            continue
        df = pd.read_csv(path, sep="\t")
        if df.empty:
            continue
        df = df.copy()
        df["chrom"] = df["chrom"].map(normalize_chrom)
        df["state"] = df["state"].map(normalize_state)
        if "sample_id" not in df.columns:
            df["sample_id"] = path.stem.split(".")[0]
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    merged = pd.concat(frames, ignore_index=True)
    if "keep_event" in merged.columns:
        merged = merged[merged["keep_event"].fillna(0).astype(int) == 1].copy()
    if "high_risk_bin_fraction" in merged.columns:
        merged["risk_bucket"] = np.where(merged["high_risk_bin_fraction"].fillna(0.0) >= 0.50, "high_risk", "clean_or_mixed")
    return merged


def top_branch_b_row(frame):
    if frame.empty:
        return None
    ordered = frame.assign(_abs_score=event_support_z(frame)).sort_values("_abs_score", ascending=False)
    return ordered.iloc[0]


def event_support_z(frame):
    support_columns = [
        column for column in ["calibrated_mean_z", "calibrated_median_z", "event_corr_adjusted_z"] if column in frame.columns
    ]
    if not support_columns or frame.empty:
        return pd.Series(np.nan, index=frame.index, dtype=float)
    support_frame = frame[support_columns].apply(pd.to_numeric, errors="coerce").abs()
    return support_frame.max(axis=1, skipna=True)


def extract_branch_b_fraction(row):
    if row is None:
        return {
            "point": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
            "status": "",
            "reliable": 0,
        }
    return {
        "point": safe_float(row.get("biopsy_abnormal_cell_fraction_point", np.nan)),
        "ci_low": safe_float(row.get("biopsy_abnormal_cell_fraction_ci_low", np.nan)),
        "ci_high": safe_float(row.get("biopsy_abnormal_cell_fraction_ci_high", np.nan)),
        "status": str(row.get("biopsy_abnormal_cell_fraction_status", "")),
        "reliable": int(safe_float(row.get("biopsy_abnormal_cell_fraction_reliable", 0), default=0.0)),
    }


def classify_fraction_bin(value):
    if not np.isfinite(value):
        return "unknown"
    if value <= 0.05:
        return "<=5%"
    if value <= 0.10:
        return "5-10%"
    if value <= 0.15:
        return "10-15%"
    if value <= 0.20:
        return "15-20%"
    if value <= 0.30:
        return "20-30%"
    return ">30%"


def load_a_branch_events(paths):
    rows = []
    for path_value in paths:
        path = Path(path_value)
        if not path.exists():
            continue
        df = pd.read_csv(path, sep="\t")
        if df.empty:
            continue
        sample_id = path.name.split("_")[0]
        df = df.copy()
        df["sample_id"] = sample_id
        df["chrom"] = df["chr"].map(normalize_chrom)
        df["state"] = df["type"].map(normalize_state)
        df["abs_zscore"] = df["zscore"].astype(float).abs()
        rows.append(df[["sample_id", "chrom", "state", "start", "end", "zscore", "abs_zscore"]])
    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True)


def overlap_fraction(left_start, left_end, right_start, right_end):
    if any(pd.isna(value) for value in [left_start, left_end, right_start, right_end]):
        return np.nan
    overlap = max(0.0, min(float(left_end), float(right_end)) - max(float(left_start), float(right_start)) + 1.0)
    length = max(float(left_end) - float(left_start) + 1.0, 1.0)
    return overlap / length


def select_matches(df, sample_id, chrom, state, truth_start, truth_end, start_col, end_col):
    if df.empty:
        return df
    subset = df[
        df["sample_id"].astype(str).eq(str(sample_id))
        & df["chrom"].astype(str).eq(str(chrom))
        & df["state"].astype(str).eq(str(state))
    ].copy()
    if subset.empty:
        return subset
    if pd.notna(truth_start) and pd.notna(truth_end) and start_col in subset.columns and end_col in subset.columns:
        subset["truth_overlap_fraction"] = subset.apply(
            lambda row: overlap_fraction(truth_start, truth_end, row[start_col], row[end_col]),
            axis=1,
        )
        subset = subset[subset["truth_overlap_fraction"].fillna(0.0) > 0.0].copy()
    return subset


def build_fraction_summary(simulation_df):
    if simulation_df.empty or "truth_fraction" not in simulation_df.columns:
        return {}
    truth_count = int(simulation_df["truth_fraction"].notna().sum())
    fraction_eval = simulation_df[
        simulation_df["truth_fraction"].notna()
        & simulation_df["branch_b_top_fraction"].notna()
    ].copy()
    if fraction_eval.empty:
        return (
            {
                "truth_fraction_event_count": truth_count,
                "evaluable_event_count": 0,
                "reliable_event_count": 0,
                "mae": None,
                "rmse": None,
                "median_abs_error": None,
                "ci_coverage": None,
            }
            if truth_count > 0
            else {}
        )
    squared_error = np.square(fraction_eval["branch_b_fraction_error"].astype(float))
    coverage_values = fraction_eval["branch_b_fraction_ci_contains_truth"].dropna().astype(float)
    return {
        "truth_fraction_event_count": truth_count,
        "evaluable_event_count": int(len(fraction_eval)),
        "reliable_event_count": int(fraction_eval["branch_b_top_fraction_reliable"].fillna(0).astype(int).sum()),
        "mae": float(fraction_eval["branch_b_fraction_abs_error"].mean()),
        "rmse": float(np.sqrt(np.mean(squared_error))),
        "median_abs_error": float(fraction_eval["branch_b_fraction_abs_error"].median()),
        "ci_coverage": float(coverage_values.mean()) if not coverage_values.empty else None,
    }


def build_low_fraction_detection(admixture_df, thresholds):
    rows = []
    if admixture_df.empty or "simulated_target_fraction" not in admixture_df.columns:
        return rows
    for threshold in thresholds:
        subset = admixture_df[
            admixture_df["simulated_target_fraction"].notna()
            & (admixture_df["simulated_target_fraction"].astype(float) <= float(threshold))
        ].copy()
        rows.append(
            {
                "fraction_threshold": float(threshold),
                "branch_b_detection_rate": float(subset["branch_b_detected"].mean()) if not subset.empty else None,
                "a_branch_detection_rate": float(subset["a_branch_detected"].mean()) if not subset.empty else None,
                "truth_event_count": int(len(subset)),
            }
        )
    return rows


def main():
    args = parse_args()
    logger = setup_logger("cnv_benchmark", args.log or None)
    truth_df = load_truth_table(args.truth_tsv)
    branch_b_df = load_branch_b_events(args.event_tsv)
    branch_a_df = load_a_branch_events(args.a_branch_bed)
    admixture_levels = args.admixture_level or [1.0, 0.75, 0.5, 0.3, 0.2, 0.1]
    low_fraction_thresholds = sorted(
        {float(item) for item in (args.low_fraction_threshold or [0.05, 0.10, 0.15, 0.20, 0.30])}
    )

    if truth_df.empty:
        empty = pd.DataFrame(columns=["sample_id"])
        ensure_parent(args.output_simulation)
        empty.to_csv(args.output_simulation, sep="\t", index=False)
        empty.to_csv(args.output_admixture, sep="\t", index=False)
        ensure_parent(args.output_summary).write_text(
            json.dumps({"status": "skipped", "reason": "truth_unavailable"}, indent=2),
            encoding="utf-8",
        )
        logger.info("benchmark skipped: truth table unavailable")
        return

    simulation_rows = []
    admixture_rows = []
    for index, truth in enumerate(truth_df.itertuples(index=False), start=1):
        sample_id = getattr(truth, "sample_id")
        chrom = normalize_chrom(getattr(truth, "chrom"))
        state = normalize_state(getattr(truth, "expected_state"))
        truth_start = getattr(truth, "start", np.nan)
        truth_end = getattr(truth, "end", np.nan)
        truth_fraction = safe_float(getattr(truth, "truth_fraction", np.nan))
        truth_fraction_source = str(getattr(truth, "truth_fraction_source", "") or "")

        branch_b_matches = select_matches(branch_b_df, sample_id, chrom, state, truth_start, truth_end, "start", "end")
        branch_a_matches = select_matches(branch_a_df, sample_id, chrom, state, truth_start, truth_end, "start", "end")
        branch_b_top_row = top_branch_b_row(branch_b_matches)
        branch_b_fraction = extract_branch_b_fraction(branch_b_top_row)

        branch_b_top_z = float(event_support_z(branch_b_matches).max()) if not branch_b_matches.empty else np.nan
        branch_a_top_z = float(branch_a_matches["abs_zscore"].max()) if not branch_a_matches.empty else np.nan
        fraction_error = (
            branch_b_fraction["point"] - truth_fraction
            if np.isfinite(branch_b_fraction["point"]) and np.isfinite(truth_fraction)
            else np.nan
        )
        ci_contains_truth = (
            int(branch_b_fraction["ci_low"] <= truth_fraction <= branch_b_fraction["ci_high"])
            if np.isfinite(truth_fraction) and np.isfinite(branch_b_fraction["ci_low"]) and np.isfinite(branch_b_fraction["ci_high"])
            else np.nan
        )
        simulation_rows.append(
            {
                "truth_id": f"truth_{index}",
                "sample_id": sample_id,
                "chrom": chrom,
                "expected_state": state,
                "truth_start": truth_start,
                "truth_end": truth_end,
                "truth_fraction": truth_fraction,
                "truth_fraction_source": truth_fraction_source,
                "branch_b_match_count": int(len(branch_b_matches)),
                "branch_b_top_abs_z": branch_b_top_z,
                "branch_b_detected": int(np.isfinite(branch_b_top_z) and branch_b_top_z >= args.branch_b_z_threshold),
                "branch_b_top_fraction": branch_b_fraction["point"],
                "branch_b_top_fraction_ci_low": branch_b_fraction["ci_low"],
                "branch_b_top_fraction_ci_high": branch_b_fraction["ci_high"],
                "branch_b_top_fraction_status": branch_b_fraction["status"],
                "branch_b_top_fraction_reliable": branch_b_fraction["reliable"],
                "branch_b_fraction_error": fraction_error,
                "branch_b_fraction_abs_error": abs(fraction_error) if np.isfinite(fraction_error) else np.nan,
                "branch_b_fraction_ci_contains_truth": ci_contains_truth,
                "a_branch_match_count": int(len(branch_a_matches)),
                "a_branch_top_abs_z": branch_a_top_z,
                "a_branch_detected": int(np.isfinite(branch_a_top_z) and branch_a_top_z >= args.branch_a_z_threshold),
            }
        )

        for admixture_level in admixture_levels:
            simulated_branch_b = branch_b_top_z * admixture_level if np.isfinite(branch_b_top_z) else np.nan
            simulated_branch_a = branch_a_top_z * admixture_level if np.isfinite(branch_a_top_z) else np.nan
            simulated_truth_fraction = truth_fraction * admixture_level if np.isfinite(truth_fraction) else np.nan
            simulated_fraction_proxy = (
                branch_b_fraction["point"] * admixture_level
                if np.isfinite(branch_b_fraction["point"])
                else np.nan
            )
            simulated_target_fraction = (
                simulated_truth_fraction
                if np.isfinite(simulated_truth_fraction)
                else simulated_fraction_proxy
            )
            simulated_fraction_source = (
                "truth_fraction"
                if np.isfinite(simulated_truth_fraction)
                else ("branch_b_estimate" if np.isfinite(simulated_fraction_proxy) else "")
            )
            admixture_rows.append(
                {
                    "truth_id": f"truth_{index}",
                    "sample_id": sample_id,
                    "chrom": chrom,
                    "expected_state": state,
                    "admixture_level": float(admixture_level),
                    "simulated_branch_b_abs_z": simulated_branch_b,
                    "branch_b_detected": int(np.isfinite(simulated_branch_b) and simulated_branch_b >= args.branch_b_z_threshold),
                    "simulated_a_branch_abs_z": simulated_branch_a,
                    "a_branch_detected": int(np.isfinite(simulated_branch_a) and simulated_branch_a >= args.branch_a_z_threshold),
                    "simulated_truth_fraction": simulated_truth_fraction,
                    "simulated_branch_b_fraction_proxy": simulated_fraction_proxy,
                    "simulated_target_fraction": simulated_target_fraction,
                    "simulated_target_fraction_source": simulated_fraction_source,
                    "simulated_fraction_bin": classify_fraction_bin(simulated_target_fraction),
                }
            )

    simulation_df = pd.DataFrame(simulation_rows)
    admixture_df = pd.DataFrame(admixture_rows)
    ensure_parent(args.output_simulation)
    simulation_df.to_csv(args.output_simulation, sep="\t", index=False)
    admixture_df.to_csv(args.output_admixture, sep="\t", index=False)

    detection_by_level = []
    if not admixture_df.empty:
        for admixture_level, group in admixture_df.groupby("admixture_level", dropna=False):
            detection_by_level.append(
                {
                    "admixture_level": float(admixture_level),
                    "branch_b_detection_rate": float(group["branch_b_detected"].mean()),
                    "a_branch_detection_rate": float(group["a_branch_detected"].mean()),
                    "truth_event_count": int(len(group)),
                }
            )

    summary = {
        "status": "completed",
        "truth_event_count": int(len(simulation_df)),
        "branch_b_detected_truth_events": int(simulation_df["branch_b_detected"].sum()) if not simulation_df.empty else 0,
        "a_branch_detected_truth_events": int(simulation_df["a_branch_detected"].sum()) if not simulation_df.empty else 0,
        "branch_b_detection_rate": float(simulation_df["branch_b_detected"].mean()) if not simulation_df.empty else None,
        "a_branch_detection_rate": float(simulation_df["a_branch_detected"].mean()) if not simulation_df.empty else None,
        "admixture_detection": detection_by_level,
        "low_fraction_detection": build_low_fraction_detection(admixture_df, low_fraction_thresholds),
        "branch_b_z_threshold": float(args.branch_b_z_threshold),
        "a_branch_z_threshold": float(args.branch_a_z_threshold),
    }
    fraction_summary = build_fraction_summary(simulation_df)
    if fraction_summary:
        summary["fraction_estimation"] = fraction_summary
    if not branch_b_df.empty and "risk_bucket" in branch_b_df.columns:
        summary["branch_b_event_risk_buckets"] = branch_b_df["risk_bucket"].value_counts(dropna=False).to_dict()
    ensure_parent(args.output_summary).write_text(json.dumps(summary, indent=2), encoding="utf-8")
    logger.info("benchmark completed: truth_events=%d", len(simulation_df))


if __name__ == "__main__":
    main()
