#!/biosoftware/miniconda/envs/snakemake_env/bin/python
import argparse
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

from pgta.predict.branch_b.common import (
    OPTIONAL_REGION_BOOL_COLUMNS,
    OPTIONAL_REGION_FRACTION_COLUMNS,
    annotate_region_risk,
    interval_overlap,
    load_sample_bins,
    read_table,
    write_json,
    write_table,
)
from pgta.core.logging import setup_logger


def parse_args():
    parser = argparse.ArgumentParser(description="Explicit Branch B bin correction using 2D local regression on GC and mappability proxy.")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--npz", required=True)
    parser.add_argument("--annotations", required=True)
    parser.add_argument("--combined-mask", required=True)
    parser.add_argument("--output-bins", required=True)
    parser.add_argument("--output-summary", required=True)
    parser.add_argument("--correction-model", default="2d_loess_gc_mappability")
    parser.add_argument("--loess-frac", type=float, default=0.2)
    parser.add_argument("--min-valid-bins", type=int, default=200)
    parser.add_argument("--robust-iters", type=int, default=2)
    parser.add_argument("--include-mask-label", action="append", default=[])
    parser.add_argument("--log", default="")
    return parser.parse_args()


ANNOTATION_VALUE_COLUMNS = [
    "gc_fraction",
    "atgc_fraction",
    "n_fraction",
    "effective_size",
    "mappability_score",
    *OPTIONAL_REGION_FRACTION_COLUMNS.keys(),
]
ANNOTATION_BOOL_COLUMNS = ["is_autosome", *OPTIONAL_REGION_BOOL_COLUMNS.keys()]
MASK_LABEL_PRIORITY = {"pass": 0, "soft": 1, "dynamic": 2, "hard": 3}
ACROCENTRIC_CHROMS = {"chr13", "chr14", "chr15", "chr21", "chr22"}
REPEAT_HOTSPOT_CHROMS = {"chr9", "chr16"}


def aggregate_reference_to_sample_bins(sample_bins_df, source_df, value_columns, bool_columns=None):
    if source_df.empty:
        return pd.DataFrame(columns=["chrom", "start", "end", *value_columns, *(bool_columns or [])])

    bool_columns = list(bool_columns or [])
    target_bins = sample_bins_df[["chrom", "start", "end"]].drop_duplicates().sort_values(["chrom", "start", "end"]).copy()
    aggregated_rows = []
    for chrom, chrom_targets in target_bins.groupby("chrom", sort=False):
        chrom_source = source_df[source_df["chrom"].astype(str).eq(str(chrom))].sort_values(["start", "end"]).copy()
        if chrom_source.empty:
            for row in chrom_targets.itertuples(index=False):
                payload = {"chrom": row.chrom, "start": int(row.start), "end": int(row.end)}
                payload.update({column: np.nan for column in value_columns})
                payload.update({column: 0 for column in bool_columns})
                aggregated_rows.append(payload)
            continue

        source_starts = chrom_source["start"].to_numpy(dtype=np.int64)
        source_ends = chrom_source["end"].to_numpy(dtype=np.int64)
        source_values = {
            column: pd.to_numeric(chrom_source[column], errors="coerce").to_numpy(dtype=np.float64)
            for column in value_columns
            if column in chrom_source.columns
        }
        source_bools = {
            column: chrom_source[column].fillna(0).astype(int).to_numpy(dtype=np.int64)
            for column in bool_columns
            if column in chrom_source.columns
        }
        source_ptr = 0
        for row in chrom_targets.itertuples(index=False):
            start = int(row.start)
            end = int(row.end)
            span = max(end - start, 1)
            while source_ptr < len(chrom_source) and source_ends[source_ptr] <= start:
                source_ptr += 1
            cursor = source_ptr
            value_totals = {column: 0.0 for column in value_columns}
            bool_hits = {column: 0 for column in bool_columns}
            while cursor < len(chrom_source) and source_starts[cursor] < end:
                overlap_bp = interval_overlap(start, end, source_starts[cursor], source_ends[cursor])
                if overlap_bp > 0:
                    overlap_fraction = overlap_bp / span
                    for column, values in source_values.items():
                        value = values[cursor]
                        if np.isfinite(value):
                            value_totals[column] += float(value) * overlap_fraction
                    for column, values in source_bools.items():
                        if int(values[cursor]) > 0:
                            bool_hits[column] = 1
                cursor += 1
            payload = {"chrom": row.chrom, "start": start, "end": end}
            payload.update(value_totals)
            payload.update(bool_hits)
            aggregated_rows.append(payload)
    return pd.DataFrame(aggregated_rows)


def aggregate_mask_to_sample_bins(sample_bins_df, mask_df):
    if mask_df.empty:
        return pd.DataFrame(columns=["chrom", "start", "end", "mask_label", "mask_reason"])

    target_bins = sample_bins_df[["chrom", "start", "end"]].drop_duplicates().sort_values(["chrom", "start", "end"]).copy()
    aggregated_rows = []
    for chrom, chrom_targets in target_bins.groupby("chrom", sort=False):
        chrom_mask = mask_df[mask_df["chrom"].astype(str).eq(str(chrom))].sort_values(["start", "end"]).copy()
        if chrom_mask.empty:
            for row in chrom_targets.itertuples(index=False):
                aggregated_rows.append(
                    {"chrom": row.chrom, "start": int(row.start), "end": int(row.end), "mask_label": "pass", "mask_reason": ""}
                )
            continue

        mask_starts = chrom_mask["start"].to_numpy(dtype=np.int64)
        mask_ends = chrom_mask["end"].to_numpy(dtype=np.int64)
        mask_labels = chrom_mask["mask_label"].fillna("pass").astype(str).to_numpy(dtype=object)
        mask_reasons = chrom_mask["mask_reason"].fillna("").astype(str).to_numpy(dtype=object)
        mask_ptr = 0
        for row in chrom_targets.itertuples(index=False):
            start = int(row.start)
            end = int(row.end)
            span = max(end - start, 1)
            while mask_ptr < len(chrom_mask) and mask_ends[mask_ptr] <= start:
                mask_ptr += 1
            cursor = mask_ptr
            label_bp = Counter()
            reason_bp = Counter()
            while cursor < len(chrom_mask) and mask_starts[cursor] < end:
                overlap_bp = interval_overlap(start, end, mask_starts[cursor], mask_ends[cursor])
                if overlap_bp > 0:
                    label = str(mask_labels[cursor] or "pass")
                    label_bp[label] += overlap_bp
                    reason = str(mask_reasons[cursor] or "")
                    if reason:
                        reason_bp[reason] += overlap_bp
                cursor += 1

            label_fraction = {label: bp / span for label, bp in label_bp.items()}
            if label_fraction.get("hard", 0.0) >= 0.50:
                mask_label = "hard"
            elif label_fraction.get("dynamic", 0.0) >= 0.30:
                mask_label = "dynamic"
            elif label_fraction.get("soft", 0.0) >= 0.30:
                mask_label = "soft"
            elif label_fraction.get("hard", 0.0) > 0.0:
                mask_label = "dynamic"
            elif label_fraction.get("dynamic", 0.0) > 0.0:
                mask_label = "dynamic"
            elif label_fraction.get("soft", 0.0) > 0.0:
                mask_label = "soft"
            else:
                mask_label = "pass"

            dominant_reason = ""
            if reason_bp:
                dominant_reason = max(reason_bp.items(), key=lambda item: (item[1], len(str(item[0]))))[0]
            aggregated_rows.append(
                {"chrom": row.chrom, "start": start, "end": end, "mask_label": mask_label, "mask_reason": dominant_reason}
            )
    return pd.DataFrame(aggregated_rows)


def escalate_hotspot_masks(reference_df):
    frame = reference_df.copy()
    homology_fraction = np.maximum.reduce(
        [
            frame["par_overlap_fraction"].to_numpy(dtype=np.float64),
            frame["xtr_overlap_fraction"].to_numpy(dtype=np.float64),
            frame["sex_homology_overlap_fraction"].to_numpy(dtype=np.float64),
            frame["ambiguous_alignment_overlap_fraction"].to_numpy(dtype=np.float64),
        ]
    )
    acrocentric_gap = (
        frame["chrom"].astype(str).isin(ACROCENTRIC_CHROMS).to_numpy(dtype=bool)
        & (frame["gap_centromere_telomere_overlap_fraction"].to_numpy(dtype=np.float64) >= 0.50)
    )
    repeat_hotspot = (
        frame["chrom"].astype(str).isin(REPEAT_HOTSPOT_CHROMS).to_numpy(dtype=bool)
        & (frame["segmental_duplication_overlap_fraction"].to_numpy(dtype=np.float64) >= 0.35)
    )
    xy_homology = (
        frame["chrom"].astype(str).isin({"chrX", "chrY"}).to_numpy(dtype=bool)
        & (homology_fraction >= 0.10)
    )
    strong_xy_homology = (
        frame["chrom"].astype(str).isin({"chrX", "chrY"}).to_numpy(dtype=bool)
        & (homology_fraction >= 0.50)
    )

    reasons = frame["mask_reason"].fillna("").astype(str)
    labels = frame["mask_label"].fillna("pass").astype(str)

    dynamic_hotspot = repeat_hotspot | xy_homology
    hard_hotspot = acrocentric_gap | strong_xy_homology
    for index in np.where(dynamic_hotspot | hard_hotspot)[0]:
        desired = "hard" if hard_hotspot[index] else "dynamic"
        current = str(labels.iloc[index])
        if MASK_LABEL_PRIORITY.get(desired, 0) <= MASK_LABEL_PRIORITY.get(current, 0):
            continue
        hotspot_reason = "branch_b_hotspot"
        if acrocentric_gap[index]:
            hotspot_reason = "branch_b_hotspot:acrocentric_gap"
        elif strong_xy_homology[index] or xy_homology[index]:
            hotspot_reason = "branch_b_hotspot:xy_homology"
        elif repeat_hotspot[index]:
            hotspot_reason = "branch_b_hotspot:repeat_hotspot"
        labels.iloc[index] = desired
        reasons.iloc[index] = hotspot_reason if not reasons.iloc[index] else f"{reasons.iloc[index]};{hotspot_reason}"

    frame["mask_label"] = labels
    frame["mask_reason"] = reasons
    return frame


def read_reference_tables(annotation_path, combined_mask_path, sample_bins_df):
    annotations = read_table(annotation_path, required_columns=["chrom", "start", "end"], empty_ok=False)
    mask_df = read_table(combined_mask_path, required_columns=["chrom", "start", "end"], empty_ok=False)
    keep_cols = ["chrom", "start", "end", *ANNOTATION_VALUE_COLUMNS, *ANNOTATION_BOOL_COLUMNS]
    annotations = annotations[[column for column in keep_cols if column in annotations.columns]].copy()
    mask_cols = ["chrom", "start", "end", "mask_label", "mask_reason"]
    mask_df = mask_df[[column for column in mask_cols if column in mask_df.columns]].copy()

    projected_annotations = aggregate_reference_to_sample_bins(
        sample_bins_df=sample_bins_df,
        source_df=annotations,
        value_columns=ANNOTATION_VALUE_COLUMNS,
        bool_columns=ANNOTATION_BOOL_COLUMNS,
    )
    projected_mask = aggregate_mask_to_sample_bins(sample_bins_df=sample_bins_df, mask_df=mask_df)
    merged = projected_annotations.merge(projected_mask, on=["chrom", "start", "end"], how="left")
    merged["mask_label"] = merged["mask_label"].fillna("pass")
    merged["mask_reason"] = merged["mask_reason"].fillna("")
    return escalate_hotspot_masks(merged)


def tricube(distance, radius):
    if radius <= 0.0:
        return np.ones_like(distance, dtype=np.float64)
    scaled = np.clip(distance / radius, 0.0, None)
    weights = np.zeros_like(scaled, dtype=np.float64)
    inside = scaled < 1.0
    weights[inside] = np.power(1.0 - np.power(scaled[inside], 3.0), 3.0)
    return weights


def bisquare_robust_weights(residuals):
    residuals = np.asarray(residuals, dtype=np.float64)
    centered = residuals - np.median(residuals)
    mad = np.median(np.abs(centered))
    if not np.isfinite(mad) or mad <= 0.0:
        return np.ones_like(residuals, dtype=np.float64)
    scale = 6.0 * 1.4826 * mad
    if scale <= 0.0:
        return np.ones_like(residuals, dtype=np.float64)
    u = np.abs(centered) / scale
    weights = np.zeros_like(u, dtype=np.float64)
    inside = u < 1.0
    weights[inside] = np.power(1.0 - np.power(u[inside], 2.0), 2.0)
    return weights


def local_linear_predict(train_features, train_signal, query_features, base_weights, frac):
    n_train = train_features.shape[0]
    n_query = query_features.shape[0]
    if n_train == 0:
        return np.zeros(n_query, dtype=np.float64)
    neighbor_count = max(25, min(n_train, int(np.ceil(max(frac, 0.01) * n_train))))
    predictions = np.zeros(n_query, dtype=np.float64)
    for index in range(n_query):
        center = query_features[index]
        deltas = train_features - center
        distances = np.sqrt(np.sum(np.square(deltas), axis=1))
        radius = float(np.partition(distances, neighbor_count - 1)[neighbor_count - 1])
        if radius <= 0.0:
            positive = distances[distances > 0.0]
            radius = float(np.min(positive)) if positive.size else 1.0
        weights = tricube(distances, radius) * base_weights
        if not np.any(weights > 0.0):
            predictions[index] = float(np.average(train_signal, weights=np.clip(base_weights, 1e-6, None)))
            continue
        design = np.column_stack(
            [
                np.ones(n_train, dtype=np.float64),
                deltas[:, 0],
                deltas[:, 1],
            ]
        )
        weighted_design = design * np.sqrt(weights)[:, None]
        weighted_signal = train_signal * np.sqrt(weights)
        try:
            beta, _, _, _ = np.linalg.lstsq(weighted_design, weighted_signal, rcond=None)
            predictions[index] = float(beta[0])
        except np.linalg.LinAlgError:
            predictions[index] = float(np.average(train_signal, weights=np.clip(weights, 1e-6, None)))
    return predictions


def fit_surface(feature_frame, signal, base_weights, frac, robust_iters):
    robust = np.ones_like(signal, dtype=np.float64)
    train_features = feature_frame.to_numpy(dtype=np.float64)
    for _ in range(max(int(robust_iters), 1)):
        predictions = local_linear_predict(train_features, signal, train_features, base_weights * robust, frac)
        residuals = signal - predictions
        robust = bisquare_robust_weights(residuals)
    return lambda query_frame: local_linear_predict(
        train_features,
        signal,
        query_frame.to_numpy(dtype=np.float64),
        base_weights * robust,
        frac,
    )


def main():
    args = parse_args()
    logger = setup_logger("cnv_correction", args.log or None)
    bins_df, binsize, quality = load_sample_bins(args.npz)
    reference_df = read_reference_tables(args.annotations, args.combined_mask, bins_df)
    bins_df = bins_df.merge(reference_df, on=["chrom", "start", "end"], how="left")
    bins_df["gc_fraction"] = bins_df["gc_fraction"].astype(float)
    bins_df["atgc_fraction"] = bins_df["atgc_fraction"].astype(float)
    bins_df["mask_label"] = bins_df["mask_label"].fillna("pass")
    bins_df["mask_reason"] = bins_df["mask_reason"].fillna("")
    bins_df["mappability_proxy"] = bins_df.get("mappability_score", bins_df["atgc_fraction"]).fillna(1.0)
    bins_df = annotate_region_risk(bins_df)
    bins_df["signal_for_calling"] = bins_df["normalized_signal"]
    bins_df["correction_expected_signal"] = bins_df["normalized_signal"].median()
    bins_df["correction_residual"] = 0.0
    bins_df["correction_applied"] = 0

    include_labels = {label.strip() for label in args.include_mask_label if str(label).strip()}
    if not include_labels:
        include_labels = {"pass", "soft"}

    valid = bins_df[
        bins_df["chrom"].str.match(r"^chr([1-9]|1[0-9]|2[0-2])$")
        & bins_df["mask_label"].isin(include_labels)
        & bins_df["correction_fit_eligible"].eq(1)
        & bins_df["gc_fraction"].notna()
        & bins_df["mappability_proxy"].notna()
    ].copy()

    correction_status = "skipped"
    if len(valid) >= args.min_valid_bins:
        valid["gc_z"] = (valid["gc_fraction"] - valid["gc_fraction"].median()) / max(valid["gc_fraction"].std(ddof=0), 1e-6)
        valid["map_z"] = (valid["mappability_proxy"] - valid["mappability_proxy"].median()) / max(valid["mappability_proxy"].std(ddof=0), 1e-6)
        feature_frame = valid[["gc_z", "map_z"]]
        fit = fit_surface(
            feature_frame=feature_frame,
            signal=valid["normalized_signal"].to_numpy(dtype=np.float64),
            base_weights=np.clip(
                valid["bin_weight"].to_numpy(dtype=np.float64) * valid["region_risk_weight"].to_numpy(dtype=np.float64),
                0.25,
                None,
            ),
            frac=args.loess_frac,
            robust_iters=args.robust_iters,
        )

        all_features = bins_df.copy()
        all_features["gc_z"] = (all_features["gc_fraction"] - valid["gc_fraction"].median()) / max(valid["gc_fraction"].std(ddof=0), 1e-6)
        all_features["map_z"] = (all_features["mappability_proxy"] - valid["mappability_proxy"].median()) / max(valid["mappability_proxy"].std(ddof=0), 1e-6)
        all_features["gc_z"] = all_features["gc_z"].fillna(0.0)
        all_features["map_z"] = all_features["map_z"].fillna(0.0)
        expected = fit(all_features[["gc_z", "map_z"]])
        baseline_level = float(np.median(valid["normalized_signal"].to_numpy(dtype=np.float64)))
        corrected = bins_df["normalized_signal"].to_numpy(dtype=np.float64) - expected + baseline_level
        bins_df["correction_expected_signal"] = expected
        bins_df["correction_residual"] = bins_df["normalized_signal"] - bins_df["correction_expected_signal"]
        bins_df["signal_for_calling"] = corrected
        bins_df["correction_applied"] = 1
        correction_status = "applied"
    else:
        logger.warning(
            "insufficient valid bins for explicit correction: sample=%s valid_bins=%d min_valid_bins=%d",
            args.sample_id,
            len(valid),
            args.min_valid_bins,
        )

    summary = {
        "sample_id": args.sample_id,
        "binsize": int(binsize),
        "quality": None if not np.isfinite(quality) else float(quality),
        "correction_model": args.correction_model,
        "correction_status": correction_status,
        "valid_bin_count": int(len(valid)),
        "fit_eligible_bin_count": int(bins_df["correction_fit_eligible"].sum()),
        "include_mask_labels": sorted(include_labels),
        "risk_class_counts": bins_df["region_risk_class"].value_counts(dropna=False).to_dict(),
        "mappability_source": "annotation_mappability_score_or_atgc_fraction_proxy",
        "signal_median_before": float(np.median(bins_df["normalized_signal"].to_numpy(dtype=np.float64))),
        "signal_median_after": float(np.median(bins_df["signal_for_calling"].to_numpy(dtype=np.float64))),
    }
    write_table(args.output_bins, bins_df)
    write_json(args.output_summary, summary)
    logger.info(
        "correction finished: sample=%s status=%s valid_bins=%d output=%s",
        args.sample_id,
        correction_status,
        len(valid),
        args.output_bins,
    )


if __name__ == "__main__":
    main()
