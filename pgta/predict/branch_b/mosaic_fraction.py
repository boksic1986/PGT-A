#!/biosoftware/miniconda/envs/snakemake_env/bin/python
import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd

from pgta.predict.branch_b.common import effective_sample_size, read_bins_and_candidates, weighted_mean, write_json, write_table
from pgta.core.logging import setup_logger


def parse_args():
    parser = argparse.ArgumentParser(description="Estimate biopsy abnormal-cell fraction for Branch B CNV events.")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--input-bins", required=True)
    parser.add_argument("--input-candidates", required=True)
    parser.add_argument("--output-candidates", required=True)
    parser.add_argument("--output-summary", required=True)
    parser.add_argument("--min-effective-bins", type=float, default=5.0)
    parser.add_argument("--min-clean-fraction", type=float, default=0.50)
    parser.add_argument("--max-high-risk-fraction", type=float, default=0.25)
    parser.add_argument("--min-abs-log2-ratio", type=float, default=0.03)
    parser.add_argument("--low-fraction-threshold", type=float, default=0.15)
    parser.add_argument("--baseline-min-bins", type=int, default=200)
    parser.add_argument("--ci-zscore", type=float, default=1.96)
    parser.add_argument("--log", default="")
    return parser.parse_args()


def is_autosome(chrom):
    text = str(chrom).strip()
    return bool(text.startswith("chr") and text[3:].isdigit() and 1 <= int(text[3:]) <= 22)


def safe_float(value, default=np.nan):
    try:
        return float(value)
    except (TypeError, ValueError):
        return float(default)


def choose_signal_column(bins_df):
    for column in ["signal_for_calling", "normalized_signal"]:
        if column in bins_df.columns:
            return column
    raise ValueError("Neither signal_for_calling nor normalized_signal is available in bins table.")


def weighted_std(values, weights, mean_value):
    values = np.asarray(values, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    valid = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    if not np.any(valid):
        return np.nan
    variance = np.average(np.square(values[valid] - mean_value), weights=weights[valid])
    return float(math.sqrt(max(variance, 0.0)))


def build_baseline(bins_df, signal_column, args):
    baseline_mask = (
        bins_df["chrom"].map(is_autosome)
        & bins_df[signal_column].notna()
    )
    if "hmm_state" in bins_df.columns:
        baseline_mask &= bins_df["hmm_state"].astype(str).eq("neutral")
    if "calibration_null_eligible" in bins_df.columns:
        baseline_mask &= bins_df["calibration_null_eligible"].fillna(0).astype(int).eq(1)
    clean_mask = baseline_mask.copy()
    if "region_risk_class" in bins_df.columns:
        clean_mask &= bins_df["region_risk_class"].astype(str).eq("clean")

    baseline_df = bins_df.loc[clean_mask].copy()
    if len(baseline_df) < args.baseline_min_bins:
        baseline_df = bins_df.loc[baseline_mask].copy()
    if len(baseline_df) < args.baseline_min_bins:
        baseline_df = bins_df.loc[bins_df["chrom"].map(is_autosome) & bins_df[signal_column].notna()].copy()
    if baseline_df.empty:
        raise ValueError("Unable to build autosomal baseline for mosaic fraction estimation.")

    return float(np.median(baseline_df[signal_column].to_numpy(dtype=np.float64))), int(len(baseline_df))


def estimate_fraction_for_event(row, bins_df, signal_column, baseline_signal, args):
    default_result = {
        "event_copy_ratio": np.nan,
        "event_log2_ratio": np.nan,
        "neutral_baseline_signal": baseline_signal,
        "biopsy_abnormal_cell_fraction_point": np.nan,
        "biopsy_abnormal_cell_fraction_ci_low": np.nan,
        "biopsy_abnormal_cell_fraction_ci_high": np.nan,
        "biopsy_abnormal_cell_fraction_status": "not_estimated",
        "biopsy_abnormal_cell_fraction_reliable": 0,
        "biopsy_abnormal_cell_fraction_method": "",
        "biopsy_abnormal_cell_fraction_reason": "",
    }

    chrom = str(row.chrom)
    state = str(row.state).strip().lower()
    if state not in {"gain", "loss"}:
        default_result["biopsy_abnormal_cell_fraction_reason"] = "unsupported_state"
        return default_result
    if not is_autosome(chrom):
        default_result["biopsy_abnormal_cell_fraction_reason"] = "non_autosomal_event"
        return default_result

    subset = bins_df[
        bins_df["chrom"].astype(str).eq(chrom)
        & bins_df["bin_index"].between(int(row.start_bin), int(row.end_bin))
    ].copy()
    if subset.empty:
        default_result["biopsy_abnormal_cell_fraction_reason"] = "event_bins_unavailable"
        return default_result

    weights = np.clip(
        subset.get("calling_weight", subset.get("bin_weight", pd.Series(1.0, index=subset.index))).to_numpy(dtype=np.float64),
        0.10,
        None,
    )
    signal_values = subset[signal_column].to_numpy(dtype=np.float64)
    event_mean_signal = weighted_mean(signal_values, weights)
    effective_bins = effective_sample_size(weights)
    clean_fraction = safe_float(getattr(row, "clean_bin_fraction", np.nan), default=np.nan)
    high_risk_fraction = safe_float(getattr(row, "high_risk_bin_fraction", np.nan), default=np.nan)

    if not np.isfinite(event_mean_signal) or not np.isfinite(baseline_signal) or baseline_signal <= 0.0:
        default_result["biopsy_abnormal_cell_fraction_reason"] = "invalid_signal_baseline"
        return default_result

    ratio = event_mean_signal / baseline_signal
    ratio = max(ratio, 1e-6)
    log2_ratio = math.log(ratio, 2)
    default_result["event_copy_ratio"] = ratio
    default_result["event_log2_ratio"] = log2_ratio

    if effective_bins < args.min_effective_bins:
        default_result["biopsy_abnormal_cell_fraction_reason"] = "effective_bins_below_minimum"
        return default_result

    if state == "gain" and ratio <= 1.0:
        default_result["biopsy_abnormal_cell_fraction_reason"] = "gain_without_ratio_shift"
        return default_result
    if state == "loss" and ratio >= 1.0:
        default_result["biopsy_abnormal_cell_fraction_reason"] = "loss_without_ratio_shift"
        return default_result

    if abs(log2_ratio) < args.min_abs_log2_ratio:
        default_result["biopsy_abnormal_cell_fraction_reason"] = "log2_ratio_below_minimum"
        return default_result

    if state == "gain":
        fraction_point = 2.0 * (ratio - 1.0)
        method = "single_copy_autosome_gain_ratio_linear"
    else:
        fraction_point = 2.0 * (1.0 - ratio)
        method = "single_copy_autosome_loss_ratio_linear"

    signal_std = weighted_std(signal_values, weights, event_mean_signal)
    se_signal = signal_std / math.sqrt(max(effective_bins, 1.0)) if np.isfinite(signal_std) else np.nan
    se_ratio = se_signal / baseline_signal if np.isfinite(se_signal) and baseline_signal > 0.0 else np.nan
    se_fraction = 2.0 * se_ratio if np.isfinite(se_ratio) else np.nan
    ci_low = fraction_point - args.ci_zscore * se_fraction if np.isfinite(se_fraction) else np.nan
    ci_high = fraction_point + args.ci_zscore * se_fraction if np.isfinite(se_fraction) else np.nan

    fraction_point = float(np.clip(fraction_point, 0.0, 1.0))
    ci_low = float(np.clip(ci_low, 0.0, 1.0)) if np.isfinite(ci_low) else np.nan
    ci_high = float(np.clip(ci_high, 0.0, 1.0)) if np.isfinite(ci_high) else np.nan

    status = "estimated"
    reliable = 1
    reasons = []
    if np.isfinite(clean_fraction) and clean_fraction < args.min_clean_fraction:
        status = "estimated_low_confidence"
        reliable = 0
        reasons.append("clean_support_below_threshold")
    if np.isfinite(high_risk_fraction) and high_risk_fraction > args.max_high_risk_fraction:
        status = "estimated_low_confidence"
        reliable = 0
        reasons.append("high_risk_fraction_above_threshold")
    if fraction_point < args.low_fraction_threshold:
        status = "low_level_signal"
        reliable = 0
        reasons.append("fraction_below_low_level_threshold")

    return {
        "event_copy_ratio": ratio,
        "event_log2_ratio": log2_ratio,
        "neutral_baseline_signal": baseline_signal,
        "biopsy_abnormal_cell_fraction_point": fraction_point,
        "biopsy_abnormal_cell_fraction_ci_low": ci_low,
        "biopsy_abnormal_cell_fraction_ci_high": ci_high,
        "biopsy_abnormal_cell_fraction_status": status,
        "biopsy_abnormal_cell_fraction_reliable": reliable,
        "biopsy_abnormal_cell_fraction_method": method,
        "biopsy_abnormal_cell_fraction_reason": ";".join(reasons) if reasons else "supported_autosomal_single_copy_model",
    }


def main():
    args = parse_args()
    logger = setup_logger("cnv_mosaic_fraction", args.log or None)
    bins_df, events_df = read_bins_and_candidates(
        args.input_bins,
        args.input_candidates,
        bins_required_columns=["chrom", "bin_index"],
        empty_candidates_ok=True,
    )
    signal_column = choose_signal_column(bins_df)
    baseline_signal, baseline_bins = build_baseline(bins_df, signal_column, args)

    if events_df.empty:
        write_table(args.output_candidates, events_df)
        write_json(
            args.output_summary,
            {
                "sample_id": args.sample_id,
                "status": "empty",
                "signal_column": signal_column,
                "neutral_baseline_signal": baseline_signal,
                "baseline_bin_count": baseline_bins,
            },
        )
        logger.info("no candidate events to estimate fraction")
        return

    estimates = [estimate_fraction_for_event(row, bins_df, signal_column, baseline_signal, args) for row in events_df.itertuples(index=False)]
    estimate_df = pd.DataFrame(estimates)
    output_df = pd.concat([events_df.reset_index(drop=True), estimate_df], axis=1)

    summary = {
        "sample_id": args.sample_id,
        "status": "completed",
        "signal_column": signal_column,
        "neutral_baseline_signal": baseline_signal,
        "baseline_bin_count": baseline_bins,
        "event_count": int(len(output_df)),
        "estimated_event_count": int(output_df["biopsy_abnormal_cell_fraction_status"].eq("estimated").sum()),
        "low_level_event_count": int(output_df["biopsy_abnormal_cell_fraction_status"].eq("low_level_signal").sum()),
        "low_confidence_event_count": int(output_df["biopsy_abnormal_cell_fraction_status"].eq("estimated_low_confidence").sum()),
        "reliable_event_count": int(output_df["biopsy_abnormal_cell_fraction_reliable"].fillna(0).astype(int).sum()),
    }
    write_table(args.output_candidates, output_df)
    write_json(args.output_summary, summary)
    logger.info(
        "mosaic fraction finished: sample=%s events=%d reliable=%d",
        args.sample_id,
        len(output_df),
        summary["reliable_event_count"],
    )


if __name__ == "__main__":
    main()
