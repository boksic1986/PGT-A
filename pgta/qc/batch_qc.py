#!/biosoftware/miniconda/envs/snakemake_env/bin/python
from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

@dataclass(frozen=True)
class GcCorrectionConfig:
    method: str
    frac: float
    poly_degree: int
    min_valid_bins: int
    robust_iters: int


def add_threshold_args(parser):
    parser.add_argument("--min-mapped-warn", type=int, default=2000000)
    parser.add_argument("--min-mapped-fail", type=int, default=1000000)
    parser.add_argument("--max-zero-frac-warn", type=float, default=0.15)
    parser.add_argument("--max-zero-frac-fail", type=float, default=0.30)
    parser.add_argument("--max-bin-cv-warn", type=float, default=1.50)
    parser.add_argument("--max-bin-cv-fail", type=float, default=2.50)
    parser.add_argument("--max-adj-mad-warn", type=float, default=0.35)
    parser.add_argument("--max-adj-mad-fail", type=float, default=0.55)
    parser.add_argument("--max-gini-warn", type=float, default=0.28)
    parser.add_argument("--max-gini-fail", type=float, default=0.35)
    parser.add_argument("--min-pearson-warn", type=float, default=0.92)
    parser.add_argument("--min-pearson-fail", type=float, default=0.85)
    parser.add_argument("--min-spearman-warn", type=float, default=0.90)
    parser.add_argument("--min-spearman-fail", type=float, default=0.83)
    parser.add_argument("--max-median-abs-z-warn", type=float, default=1.0)
    parser.add_argument("--max-median-abs-z-fail", type=float, default=1.5)
    parser.add_argument("--max-outlier3-warn", type=float, default=0.15)
    parser.add_argument("--max-outlier3-fail", type=float, default=0.30)


def add_gc_correction_args(parser):
    parser.add_argument("--gc-correction-method", choices=["none", "loess", "poly2"], default="loess")
    parser.add_argument("--gc-correction-frac", type=float, default=0.2)
    parser.add_argument("--gc-correction-poly-degree", type=int, default=2)
    parser.add_argument("--gc-correction-min-valid-bins", type=int, default=200)
    parser.add_argument("--gc-correction-robust-iters", type=int, default=1)


def build_thresholds_from_args(args):
    from pgta.qc.sample_qc import QcThresholds

    return QcThresholds(
        min_mapped_warn=args.min_mapped_warn,
        min_mapped_fail=args.min_mapped_fail,
        max_zero_frac_warn=args.max_zero_frac_warn,
        max_zero_frac_fail=args.max_zero_frac_fail,
        max_bin_cv_warn=args.max_bin_cv_warn,
        max_bin_cv_fail=args.max_bin_cv_fail,
        max_adj_mad_warn=args.max_adj_mad_warn,
        max_adj_mad_fail=args.max_adj_mad_fail,
        max_gini_warn=args.max_gini_warn,
        max_gini_fail=args.max_gini_fail,
        min_pearson_warn=args.min_pearson_warn,
        min_pearson_fail=args.min_pearson_fail,
        min_spearman_warn=args.min_spearman_warn,
        min_spearman_fail=args.min_spearman_fail,
        max_median_abs_z_warn=args.max_median_abs_z_warn,
        max_median_abs_z_fail=args.max_median_abs_z_fail,
        max_outlier3_warn=args.max_outlier3_warn,
        max_outlier3_fail=args.max_outlier3_fail,
    )


def build_gc_correction_config_from_args(args):
    return GcCorrectionConfig(
        method=args.gc_correction_method,
        frac=args.gc_correction_frac,
        poly_degree=args.gc_correction_poly_degree,
        min_valid_bins=args.gc_correction_min_valid_bins,
        robust_iters=args.gc_correction_robust_iters,
    )


def sample_id_from_bam(bam_path):
    return Path(bam_path).name.replace(".sorted.bam", "")


def binsize_label(bin_size_bp):
    value = int(bin_size_bp)
    if value % 1000000 == 0:
        return f"{value // 1000000}mb"
    if value % 1000 == 0:
        return f"{value // 1000}kb"
    return f"{value}bp"


def build_target_spec(base_layout, gc_fraction, target_bin_size):
    import numpy as np
    import pandas as pd

    assignments = np.zeros(base_layout.total_bins, dtype=np.int64)
    target_rows = []
    base_lengths = base_layout.bin_ends - base_layout.bin_starts
    target_idx = 0

    for chrom_index, chrom in enumerate(base_layout.chroms):
        start_idx = base_layout.chrom_start_bin[chrom]
        stop_idx = start_idx + base_layout.chrom_bin_counts[chrom]
        chrom_idx = np.arange(start_idx, stop_idx, dtype=np.int64)
        chrom_starts = base_layout.bin_starts[chrom_idx]
        chrom_length = int(base_layout.chrom_lengths[chrom_index])

        start = 0
        pos = 0
        while pos < len(chrom_idx):
            end = min(start + target_bin_size, chrom_length)
            target_rows.append((chrom, int(start), int(end)))
            while pos < len(chrom_idx) and int(chrom_starts[pos]) < end:
                assignments[int(chrom_idx[pos])] = target_idx
                pos += 1
            target_idx += 1
            start = end

    valid_gc = np.isfinite(gc_fraction)
    gc_weighted = np.bincount(
        assignments,
        weights=np.where(valid_gc, gc_fraction * base_lengths, 0.0),
        minlength=target_idx,
    )
    gc_length = np.bincount(
        assignments,
        weights=np.where(valid_gc, base_lengths, 0.0),
        minlength=target_idx,
    )
    target_gc = np.divide(
        gc_weighted,
        gc_length,
        out=np.full(target_idx, np.nan, dtype=np.float64),
        where=gc_length > 0,
    )
    return {
        "bins_df": pd.DataFrame(target_rows, columns=["chrom", "start", "end"]),
        "assignments": assignments,
        "n_bins": target_idx,
        "gc_fraction": target_gc,
    }


def aggregate_count_matrix(raw_count_matrix, assignments, n_bins):
    import numpy as np

    aggregated = []
    for counts in raw_count_matrix:
        agg = np.bincount(assignments, weights=counts.astype(np.float64), minlength=n_bins)
        aggregated.append(np.rint(agg).astype(np.int64))
    return np.vstack(aggregated)


def write_sample_outputs(
    sample_dir,
    bins_df,
    gc_fraction,
    raw_counts,
    raw_signal,
    gc_trend,
    signal,
    ref_median,
    ref_sigma,
    z_scores,
    qc_row,
    ref_rows,
):
    import pandas as pd

    from pgta.qc.sample_qc import plot_gc_bias, plot_profiles

    sample_dir.mkdir(parents=True, exist_ok=True)

    pd.DataFrame([qc_row]).to_csv(sample_dir / "qc_metrics.tsv", sep="\t", index=False)

    profile_df = bins_df.copy()
    profile_df["gc_fraction"] = gc_fraction
    profile_df["raw_count"] = raw_counts
    profile_df["raw_signal_log2cpm"] = raw_signal
    profile_df["gc_correction_trend"] = gc_trend
    profile_df["signal_log2cpm"] = signal
    profile_df["ref_median"] = ref_median
    profile_df["ref_sigma"] = ref_sigma
    profile_df["z_score"] = z_scores
    profile_df.to_csv(sample_dir / "target_bin_profile.tsv", sep="\t", index=False)

    pd.DataFrame(ref_rows).to_csv(sample_dir / "reference_summary.tsv", sep="\t", index=False)
    plot_profiles(raw_signal, signal, ref_median, z_scores, sample_dir / "target_vs_ref_profile.png")
    plot_gc_bias(gc_fraction, raw_signal, signal, sample_dir / "gc_bias_plot.png")


def run_batch_bam_qc(
    bam_paths,
    reference_fasta,
    target_bin_sizes,
    mapq,
    threads,
    outdir,
    thresholds,
    gc_config,
    nest_binsize_subdir,
):
    import numpy as np
    import pysam

    from pgta.qc.sample_qc import (
        apply_gc_correction,
        build_bin_layout,
        compute_gc_bias_metrics,
        compute_gc_fraction_for_layout,
        compute_reference_stats,
        compute_uniformity_metrics,
        count_multiple_bams,
        infer_autosome_names,
        layout_to_dataframe,
        log2_cpm,
        qc_decision,
        safe_corr,
    )

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    bam_paths = [Path(path) for path in bam_paths]
    if len(bam_paths) < 2:
        raise ValueError("Batch baseline QC requires at least two BAMs.")

    target_bin_sizes = sorted({int(item) for item in target_bin_sizes})
    base_bin_size = math.gcd(*target_bin_sizes)

    with pysam.AlignmentFile(str(bam_paths[0]), "rb") as first_bam:
        chroms = infer_autosome_names(first_bam)
        chrom_lengths = [first_bam.get_reference_length(chrom) for chrom in chroms]
    layout = build_bin_layout(chroms, chrom_lengths, base_bin_size)
    gc_fraction_base = compute_gc_fraction_for_layout(Path(reference_fasta), layout)

    count_results = count_multiple_bams(
        bam_paths=bam_paths,
        chroms=layout.chroms,
        chrom_start_bin=layout.chrom_start_bin,
        bin_size=base_bin_size,
        mapq=mapq,
        threads=threads,
        total_bins=layout.total_bins,
    )

    ordered_bams = [str(path) for path in bam_paths]
    ordered_sample_ids = [sample_id_from_bam(path) for path in bam_paths]
    sample_to_bam = dict(zip(ordered_sample_ids, ordered_bams))

    raw_count_matrix = []
    mapped_fragments = []
    for bam_key in ordered_bams:
        counts, mapped = count_results[bam_key]
        raw_count_matrix.append(counts)
        mapped_fragments.append(mapped)
    raw_count_matrix = np.vstack(raw_count_matrix)
    mapped_fragments = np.asarray(mapped_fragments, dtype=np.int64)

    for target_bin_size in target_bin_sizes:
        if target_bin_size == base_bin_size:
            bins_df = layout_to_dataframe(layout)
            usable_bins = int(layout.total_bins)
            gc_fraction = gc_fraction_base
            count_matrix = raw_count_matrix
        else:
            spec = build_target_spec(layout, gc_fraction_base, target_bin_size)
            bins_df = spec["bins_df"]
            usable_bins = int(spec["n_bins"])
            gc_fraction = spec["gc_fraction"]
            count_matrix = aggregate_count_matrix(raw_count_matrix, spec["assignments"], usable_bins)

        signal_matrix_raw = np.vstack([log2_cpm(counts) for counts in count_matrix])
        signal_matrix = []
        signal_trends = []
        signal_meta = []
        for raw_signal in signal_matrix_raw:
            corrected_signal, trend, meta = apply_gc_correction(
                signal=raw_signal,
                gc_fraction=gc_fraction,
                method=gc_config.method,
                frac=gc_config.frac,
                poly_degree=gc_config.poly_degree,
                min_valid_bins=gc_config.min_valid_bins,
                robust_iters=gc_config.robust_iters,
            )
            signal_matrix.append(corrected_signal)
            signal_trends.append(trend)
            signal_meta.append(meta)
        signal_matrix = np.vstack(signal_matrix)
        signal_trends = np.vstack(signal_trends)

        bin_outdir = outdir / binsize_label(target_bin_size) if nest_binsize_subdir else outdir
        for idx, sample_id in enumerate(ordered_sample_ids):
            target_counts = count_matrix[idx]
            target_signal_raw = signal_matrix_raw[idx]
            target_gc_trend = signal_trends[idx]
            target_signal = signal_matrix[idx]
            mask = np.arange(signal_matrix.shape[0]) != idx
            ref_signal_matrix = signal_matrix[mask]
            ref_raw_counts = count_matrix[mask]
            ref_mapped = mapped_fragments[mask]
            ref_sample_ids = [sid for j, sid in enumerate(ordered_sample_ids) if j != idx]

            ref_median, ref_sigma = compute_reference_stats(ref_signal_matrix)
            pearson_r, spearman_r = safe_corr(target_signal, ref_median)
            z_scores = (target_signal - ref_median) / ref_sigma
            abs_z = np.abs(z_scores)

            uniform_metrics = compute_uniformity_metrics(target_counts, int(mapped_fragments[idx]), target_signal)
            gc_metrics_raw = compute_gc_bias_metrics(gc_fraction, target_signal_raw)
            gc_metrics = compute_gc_bias_metrics(gc_fraction, target_signal)
            ref_metrics = {
                "pearson_r": pearson_r,
                "spearman_r": spearman_r,
                "median_abs_z": float(np.median(abs_z)),
                "outlier_frac_abs_z_gt_3": float(np.mean(abs_z > 3.0)),
                "outlier_frac_abs_z_gt_5": float(np.mean(abs_z > 5.0)),
            }

            all_metrics = {}
            all_metrics.update(uniform_metrics)
            all_metrics.update(ref_metrics)
            all_metrics.update(gc_metrics)
            decision, reason = qc_decision(all_metrics, thresholds)

            qc_row = {
                "target_bam": ordered_bams[idx],
                "bin_size": int(target_bin_size),
                "mapped_fragments": int(mapped_fragments[idx]),
                "usable_bins": int(usable_bins),
                "zero_bin_fraction": uniform_metrics["zero_bin_fraction"],
                "bin_cv": uniform_metrics["bin_cv"],
                "adjacent_diff_mad": uniform_metrics["adjacent_diff_mad"],
                "gini_coefficient": uniform_metrics["gini_coefficient"],
                "pearson_r": pearson_r,
                "spearman_r": spearman_r,
                "median_abs_z": ref_metrics["median_abs_z"],
                "outlier_frac_abs_z_gt_3": ref_metrics["outlier_frac_abs_z_gt_3"],
                "outlier_frac_abs_z_gt_5": ref_metrics["outlier_frac_abs_z_gt_5"],
                "gc_fraction_mean": gc_metrics["gc_fraction_mean"],
                "gc_signal_pearson_r": gc_metrics["gc_signal_pearson_r"],
                "gc_signal_spearman_r": gc_metrics["gc_signal_spearman_r"],
                "gc_signal_slope": gc_metrics["gc_signal_slope"],
                "raw_gc_signal_pearson_r": gc_metrics_raw["gc_signal_pearson_r"],
                "raw_gc_signal_spearman_r": gc_metrics_raw["gc_signal_spearman_r"],
                "raw_gc_signal_slope": gc_metrics_raw["gc_signal_slope"],
                "gc_correction_method": gc_config.method,
                "gc_correction_applied": int(signal_meta[idx]["gc_correction_applied"]),
                "gc_correction_valid_bins": int(signal_meta[idx]["gc_correction_valid_bins"]),
                "qc_decision": decision,
                "qc_reason": reason,
            }

            ref_rows = []
            for ref_idx, ref_sample_id in enumerate(ref_sample_ids):
                ref_counts = ref_raw_counts[ref_idx]
                ref_rows.append(
                    {
                        "ref_bam": sample_to_bam[ref_sample_id],
                        "mapped_fragments": int(ref_mapped[ref_idx]),
                        "usable_bins": int(usable_bins),
                        "zero_bin_fraction": float(np.mean(ref_counts == 0)),
                    }
                )

            write_sample_outputs(
                sample_dir=bin_outdir / sample_id,
                bins_df=bins_df,
                gc_fraction=gc_fraction,
                raw_counts=target_counts,
                raw_signal=target_signal_raw,
                gc_trend=target_gc_trend,
                signal=target_signal,
                ref_median=ref_median,
                ref_sigma=ref_sigma,
                z_scores=z_scores,
                qc_row=qc_row,
                ref_rows=ref_rows,
            )
