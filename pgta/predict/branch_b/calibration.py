#!/biosoftware/miniconda/envs/snakemake_env/bin/python
import argparse
import math

import numpy as np
import pandas as pd

from pgta.predict.branch_b.common import (
    annotate_region_risk,
    autocorrelation_inflation,
    benjamini_hochberg,
    effective_sample_size,
    estimate_autocorrelation_by_chrom,
    read_bins_and_candidates,
    stouffer_weighted_z,
    weighted_mean,
    write_json,
    write_table,
)
from pgta.core.logging import setup_logger


def parse_args():
    parser = argparse.ArgumentParser(description="Calibrate candidate CNV events with robust z and empirical null.")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--input-bins", required=True)
    parser.add_argument("--input-candidates", required=True)
    parser.add_argument("--output-bins", required=True)
    parser.add_argument("--output-candidates", required=True)
    parser.add_argument("--output-summary", required=True)
    parser.add_argument("--null-quantile-low", type=float, default=0.1)
    parser.add_argument("--null-quantile-high", type=float, default=0.9)
    parser.add_argument("--min-null-bins", type=int, default=200)
    parser.add_argument("--event-z-threshold", type=float, default=2.5)
    parser.add_argument("--log", default="")
    return parser.parse_args()


EVENT_REGION_COLUMNS = [
    "par_overlap_fraction",
    "xtr_overlap_fraction",
    "sex_homology_overlap_fraction",
    "segmental_duplication_overlap_fraction",
    "low_mappability_overlap_fraction",
    "gap_centromere_telomere_overlap_fraction",
    "repeat_rich_overlap_fraction",
    "blacklist_overlap_fraction",
    "ambiguous_alignment_overlap_fraction",
]


def summarize_event_bins(event_row, bins_df, lag_correlations):
    subset = bins_df[
        bins_df["chrom"].astype(str).eq(str(event_row.chrom))
        & bins_df["bin_index"].between(int(event_row.start_bin), int(event_row.end_bin))
    ].copy()
    if subset.empty:
        return {
            "event_weighted_z": np.nan,
            "event_corr_adjusted_z": np.nan,
            "effective_bin_count": 0.0,
            "autocorr_inflation": 1.0,
            "clean_bin_fraction": np.nan,
            "moderate_risk_bin_fraction": np.nan,
            "high_risk_bin_fraction": np.nan,
            "region_risk_score_mean": np.nan,
            "region_risk_score_max": np.nan,
            "high_risk_boundary_crossing": 0,
            "region_risk_weight_mean": np.nan,
            **{column: np.nan for column in EVENT_REGION_COLUMNS},
        }

    weights = np.clip(
        subset.get("calling_weight", subset["bin_weight"] * subset["region_risk_weight"]).to_numpy(dtype=np.float64),
        0.10,
        None,
    )
    effective_size_weights = np.clip(subset.get("effective_size", pd.Series(1.0, index=subset.index)).to_numpy(dtype=np.float64), 1.0, None)
    event_weighted_z = stouffer_weighted_z(subset["calibrated_z"].to_numpy(dtype=np.float64), weights)
    effective_bin_count = effective_sample_size(weights)
    inflation = autocorrelation_inflation(max(len(subset), int(round(effective_bin_count))), lag_correlations)
    event_corr_adjusted_z = (
        float(event_weighted_z) / math.sqrt(inflation) if np.isfinite(event_weighted_z) and inflation > 0.0 else np.nan
    )

    risk_classes = subset["region_risk_class"].astype(str)
    clean_fraction = float(np.average((risk_classes == "clean").astype(float), weights=weights))
    moderate_fraction = float(np.average((risk_classes == "moderate").astype(float), weights=weights))
    high_fraction = float(np.average((risk_classes == "high").astype(float), weights=weights))
    first_risk = str(risk_classes.iloc[0])
    last_risk = str(risk_classes.iloc[-1])
    high_risk_boundary_crossing = int(
        first_risk != last_risk and {"clean", "high"}.issubset(set(risk_classes.unique()))
    )

    result = {
        "event_weighted_z": event_weighted_z,
        "event_corr_adjusted_z": event_corr_adjusted_z,
        "effective_bin_count": effective_bin_count,
        "autocorr_inflation": inflation,
        "clean_bin_fraction": clean_fraction,
        "moderate_risk_bin_fraction": moderate_fraction,
        "high_risk_bin_fraction": high_fraction,
        "region_risk_score_mean": weighted_mean(subset["region_risk_score"], weights),
        "region_risk_score_max": float(np.nanmax(subset["region_risk_score"].to_numpy(dtype=np.float64))),
        "high_risk_boundary_crossing": high_risk_boundary_crossing,
        "region_risk_weight_mean": weighted_mean(subset["region_risk_weight"], weights),
    }
    for column in EVENT_REGION_COLUMNS:
        result[column] = weighted_mean(subset[column], effective_size_weights)
    return result


def main():
    args = parse_args()
    logger = setup_logger("cnv_calibration", args.log or None)
    bins_df, events_df = read_bins_and_candidates(
        args.input_bins,
        args.input_candidates,
        bins_required_columns=["chrom", "bin_index", "robust_z"],
        empty_candidates_ok=True,
    )
    bins_df = annotate_region_risk(bins_df)

    neutral_bins = bins_df[
        bins_df["hmm_state"].astype(str).eq("neutral")
        & bins_df["calibration_null_eligible"].fillna(0).astype(int).eq(1)
    ].copy()
    if len(neutral_bins) < args.min_null_bins:
        null_pool = bins_df[bins_df["calibration_null_eligible"].fillna(0).astype(int).eq(1)].copy()
        if null_pool.empty:
            null_pool = bins_df.copy()
        low = null_pool["robust_z"].quantile(args.null_quantile_low)
        high = null_pool["robust_z"].quantile(args.null_quantile_high)
        neutral_bins = null_pool[(null_pool["robust_z"] >= low) & (null_pool["robust_z"] <= high)].copy()
    if len(neutral_bins) < args.min_null_bins:
        neutral_bins = bins_df.copy()

    null_values = neutral_bins["robust_z"].to_numpy(dtype=np.float64)
    null_median = float(np.median(null_values))
    null_mad = float(np.median(np.abs(null_values - null_median)))
    null_scale = max(1.4826 * null_mad, float(np.std(null_values)), 1e-6)
    bins_df["calibrated_z"] = (bins_df["robust_z"] - null_median) / null_scale
    neutral_bins = bins_df.loc[neutral_bins.index].copy()
    lag_correlations = estimate_autocorrelation_by_chrom(neutral_bins, "calibrated_z")

    if not events_df.empty:
        calibrated_mean = []
        calibrated_median = []
        weighted_z = []
        corr_adjusted_z = []
        effective_bin_count = []
        autocorr_inflation_values = []
        clean_bin_fraction = []
        moderate_risk_bin_fraction = []
        high_risk_bin_fraction = []
        region_risk_score_mean = []
        region_risk_score_max = []
        region_risk_weight_mean = []
        high_risk_boundary_crossing = []
        overlap_metrics = {column: [] for column in EVENT_REGION_COLUMNS}
        pvalues = []
        for row in events_df.itertuples(index=False):
            mean_z = (float(row.segment_mean_robust_z) - null_median) / null_scale
            median_z = (float(row.segment_median_robust_z) - null_median) / null_scale
            event_metrics = summarize_event_bins(row, bins_df, lag_correlations)
            calibrated_mean.append(mean_z)
            calibrated_median.append(median_z)
            weighted_z.append(event_metrics["event_weighted_z"])
            corr_adjusted_z.append(event_metrics["event_corr_adjusted_z"])
            effective_bin_count.append(event_metrics["effective_bin_count"])
            autocorr_inflation_values.append(event_metrics["autocorr_inflation"])
            clean_bin_fraction.append(event_metrics["clean_bin_fraction"])
            moderate_risk_bin_fraction.append(event_metrics["moderate_risk_bin_fraction"])
            high_risk_bin_fraction.append(event_metrics["high_risk_bin_fraction"])
            region_risk_score_mean.append(event_metrics["region_risk_score_mean"])
            region_risk_score_max.append(event_metrics["region_risk_score_max"])
            region_risk_weight_mean.append(event_metrics["region_risk_weight_mean"])
            high_risk_boundary_crossing.append(event_metrics["high_risk_boundary_crossing"])
            for column in EVENT_REGION_COLUMNS:
                overlap_metrics[column].append(event_metrics[column])
            observed = abs(event_metrics["event_corr_adjusted_z"]) if np.isfinite(event_metrics["event_corr_adjusted_z"]) else 0.0
            pvalues.append(math.erfc(observed / math.sqrt(2.0)))
        events_df["calibrated_mean_z"] = calibrated_mean
        events_df["calibrated_median_z"] = calibrated_median
        events_df["event_weighted_z"] = weighted_z
        events_df["event_corr_adjusted_z"] = corr_adjusted_z
        events_df["effective_bin_count"] = effective_bin_count
        events_df["autocorr_inflation"] = autocorr_inflation_values
        events_df["clean_bin_fraction"] = clean_bin_fraction
        events_df["moderate_risk_bin_fraction"] = moderate_risk_bin_fraction
        events_df["high_risk_bin_fraction"] = high_risk_bin_fraction
        events_df["region_risk_score_mean"] = region_risk_score_mean
        events_df["region_risk_score_max"] = region_risk_score_max
        events_df["region_risk_weight_mean"] = region_risk_weight_mean
        events_df["high_risk_boundary_crossing"] = high_risk_boundary_crossing
        for column, values in overlap_metrics.items():
            events_df[column] = values
        events_df["empirical_pvalue"] = pvalues
        events_df["empirical_qvalue"] = benjamini_hochberg(events_df["empirical_pvalue"].to_numpy(dtype=np.float64))
        events_df["keep_event"] = (
            (events_df["empirical_qvalue"] <= 0.25)
            | (events_df["calibrated_mean_z"].abs() >= args.event_z_threshold)
            | (events_df["calibrated_median_z"].abs() >= args.event_z_threshold)
            | (events_df["event_corr_adjusted_z"].abs() >= args.event_z_threshold)
        ).astype(int)

    summary = {
        "sample_id": args.sample_id,
        "null_bin_count": int(len(neutral_bins)),
        "null_median": null_median,
        "null_scale": null_scale,
        "lag_correlations": {str(key): value for key, value in lag_correlations.items()},
        "event_count": int(len(events_df)),
        "retained_event_count": int(events_df["keep_event"].sum()) if "keep_event" in events_df.columns else 0,
    }
    write_table(args.output_bins, bins_df)
    write_table(args.output_candidates, events_df)
    write_json(args.output_summary, summary)
    logger.info("wrote calibrated bins=%s events=%s", args.output_bins, args.output_candidates)


if __name__ == "__main__":
    main()
