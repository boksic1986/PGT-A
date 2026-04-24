#!/biosoftware/miniconda/envs/snakemake_env/bin/python
# -*- coding: utf-8 -*-

import argparse
import logging
import math
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
from scipy.stats import pearsonr, spearmanr

LOGGER = logging.getLogger("bam_uniformity_qc")


@dataclass
class QcThresholds:
    min_mapped_warn: int
    min_mapped_fail: int
    max_zero_frac_warn: float
    max_zero_frac_fail: float
    max_bin_cv_warn: float
    max_bin_cv_fail: float
    max_adj_mad_warn: float
    max_adj_mad_fail: float
    max_gini_warn: float
    max_gini_fail: float
    min_pearson_warn: float
    min_pearson_fail: float
    min_spearman_warn: float
    min_spearman_fail: float
    max_median_abs_z_warn: float
    max_median_abs_z_fail: float
    max_outlier3_warn: float
    max_outlier3_fail: float


@dataclass
class BinLayout:
    chroms: List[str]
    chrom_lengths: np.ndarray
    bin_size: int
    chrom_start_bin: Dict[str, int]
    chrom_bin_counts: Dict[str, int]
    chrom_offsets: Dict[str, int]
    total_bins: int
    chrom_index_by_bin: np.ndarray
    bin_starts: np.ndarray
    bin_ends: np.ndarray


@dataclass
class BamCountTask:
    bam_path: str
    chroms: List[str]
    chrom_start_bin: Dict[str, int]
    bin_size: int
    mapq: int
    total_bins: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="BAM coverage uniformity and reference-deviation QC for low-pass WGS/CNV scenarios"
    )
    parser.add_argument("--target-bam", required=True, help="Target BAM path")
    parser.add_argument("--ref-bams", required=True, nargs="+", help="Reference BAM paths")
    parser.add_argument("--reference-fasta", required=True, help="Reference FASTA path for per-bin GC calculation")
    parser.add_argument("--bin-size", type=int, default=200000, help="Bin size in bases")
    parser.add_argument("--mapq", type=int, default=30, help="Minimum MAPQ")
    parser.add_argument("--threads", type=int, default=1, help="Parallel workers used to count target/reference BAMs")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--log", default="", help="Optional log file")
    parser.add_argument(
        "--gc-correction-method",
        choices=["none", "loess", "poly2"],
        default="loess",
        help="GC correction method applied before profile-level QC metrics",
    )
    parser.add_argument("--gc-correction-frac", type=float, default=0.2, help="LOESS neighborhood fraction")
    parser.add_argument(
        "--gc-correction-poly-degree",
        type=int,
        default=2,
        help="Polynomial degree for polynomial GC correction",
    )
    parser.add_argument(
        "--gc-correction-min-valid-bins",
        type=int,
        default=200,
        help="Minimum finite GC bins required to fit GC correction",
    )
    parser.add_argument(
        "--gc-correction-robust-iters",
        type=int,
        default=1,
        help="Additional robust reweighting iterations for LOESS GC correction",
    )

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

    return parser.parse_args()


def setup_logger(log_path: str) -> None:
    handlers: List[logging.Handler] = [logging.StreamHandler()]
    if log_path:
        log_file = Path(log_path)
        log_file.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_file, encoding="utf-8"))
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=handlers,
    )


def infer_autosome_names(bam: pysam.AlignmentFile) -> List[str]:
    ref_names = list(bam.references)
    ref_set = set(ref_names)
    if all(f"chr{i}" in ref_set for i in range(1, 23)):
        return [f"chr{i}" for i in range(1, 23)]
    if all(str(i) in ref_set for i in range(1, 23)):
        return [str(i) for i in range(1, 23)]

    autosomes = []
    for name in ref_names:
        token = name[3:] if name.startswith("chr") else name
        if token.isdigit() and 1 <= int(token) <= 22:
            autosomes.append(name)
    autosomes_sorted = sorted(autosomes, key=lambda x: int(x[3:] if x.startswith("chr") else x))
    if not autosomes_sorted:
        raise ValueError("No autosome chromosomes (1-22 / chr1-chr22) found in BAM header")
    return autosomes_sorted


def build_bin_layout(chroms: Sequence[str], chrom_lengths: Sequence[int], bin_size: int) -> BinLayout:
    chrom_list = [str(chrom) for chrom in chroms]
    lengths = np.asarray(chrom_lengths, dtype=np.int64)
    chrom_offsets: Dict[str, int] = {}
    chrom_start_bin: Dict[str, int] = {}
    chrom_bin_counts: Dict[str, int] = {}

    total_bins = 0
    total_length = int(np.sum(lengths))
    total_bins_estimate = int(sum(int(math.ceil(length / float(bin_size))) for length in lengths))
    bin_starts = np.empty(total_bins_estimate, dtype=np.int64)
    bin_ends = np.empty(total_bins_estimate, dtype=np.int64)
    chrom_index_by_bin = np.empty(total_bins_estimate, dtype=np.int32)

    offset_bp = 0
    for chrom_index, (chrom, chrom_length) in enumerate(zip(chrom_list, lengths)):
        chrom_offsets[chrom] = int(offset_bp)
        chrom_start_bin[chrom] = int(total_bins)
        n_bins = int(math.ceil(int(chrom_length) / float(bin_size)))
        chrom_bin_counts[chrom] = n_bins
        starts = np.arange(n_bins, dtype=np.int64) * int(bin_size)
        ends = np.minimum(starts + int(bin_size), int(chrom_length))
        next_total = total_bins + n_bins
        bin_starts[total_bins:next_total] = starts
        bin_ends[total_bins:next_total] = ends
        chrom_index_by_bin[total_bins:next_total] = chrom_index
        total_bins = next_total
        offset_bp += int(chrom_length)

    return BinLayout(
        chroms=chrom_list,
        chrom_lengths=lengths,
        bin_size=int(bin_size),
        chrom_start_bin=chrom_start_bin,
        chrom_bin_counts=chrom_bin_counts,
        chrom_offsets=chrom_offsets,
        total_bins=int(total_bins),
        chrom_index_by_bin=chrom_index_by_bin[:total_bins],
        bin_starts=bin_starts[:total_bins],
        bin_ends=bin_ends[:total_bins],
    )


def build_bins(bam: pysam.AlignmentFile, chroms: Sequence[str], bin_size: int) -> Tuple[pd.DataFrame, Dict[str, int]]:
    chrom_lengths = [bam.get_reference_length(chrom) for chrom in chroms]
    layout = build_bin_layout(chroms, chrom_lengths, bin_size)
    bins_df = layout_to_dataframe(layout)
    return bins_df, layout.chrom_start_bin


def layout_to_dataframe(layout: BinLayout) -> pd.DataFrame:
    chrom_array = np.asarray(layout.chroms, dtype=object)[layout.chrom_index_by_bin]
    return pd.DataFrame(
        {
            "chrom": chrom_array,
            "start": layout.bin_starts,
            "end": layout.bin_ends,
        }
    )


def compute_gc_fraction_for_layout(reference_fasta: Path, layout: BinLayout) -> np.ndarray:
    reference_fasta = Path(reference_fasta)
    if not reference_fasta.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {reference_fasta}")

    gc_values = np.full(layout.total_bins, np.nan, dtype=np.float64)
    with pysam.FastaFile(str(reference_fasta)) as fasta:
        for chrom_index, chrom in enumerate(layout.chroms):
            seq = fasta.fetch(chrom).upper().encode("ascii")
            seq_bytes = np.frombuffer(seq, dtype=np.uint8)
            valid_mask = (
                (seq_bytes == ord("A"))
                | (seq_bytes == ord("C"))
                | (seq_bytes == ord("G"))
                | (seq_bytes == ord("T"))
            )
            gc_mask = (seq_bytes == ord("G")) | (seq_bytes == ord("C"))

            valid_prefix = np.zeros(seq_bytes.size + 1, dtype=np.int64)
            gc_prefix = np.zeros(seq_bytes.size + 1, dtype=np.int64)
            np.cumsum(valid_mask.astype(np.int64), out=valid_prefix[1:])
            np.cumsum(gc_mask.astype(np.int64), out=gc_prefix[1:])

            bin_mask = layout.chrom_index_by_bin == chrom_index
            starts = layout.bin_starts[bin_mask]
            ends = layout.bin_ends[bin_mask]
            valid_counts = valid_prefix[ends] - valid_prefix[starts]
            gc_counts = gc_prefix[ends] - gc_prefix[starts]
            gc_values[bin_mask] = np.divide(
                gc_counts,
                valid_counts,
                out=np.full(valid_counts.shape, np.nan, dtype=np.float64),
                where=valid_counts > 0,
            )
    return gc_values


def compute_gc_fraction_for_bins(reference_fasta: Path, bins_df: pd.DataFrame) -> np.ndarray:
    gc_values = []
    reference_fasta = Path(reference_fasta)
    with pysam.FastaFile(str(reference_fasta)) as fasta:
        for row in bins_df.itertuples(index=False):
            seq = fasta.fetch(row.chrom, int(row.start), int(row.end)).upper()
            valid_bases = sum(base in {"A", "C", "G", "T"} for base in seq)
            if valid_bases == 0:
                gc_values.append(float("nan"))
                continue
            gc_bases = sum(base in {"G", "C"} for base in seq)
            gc_values.append(float(gc_bases) / float(valid_bases))
    return np.asarray(gc_values, dtype=np.float64)


def is_usable_read(read: pysam.AlignedSegment, mapq: int) -> bool:
    if read.is_unmapped:
        return False
    if read.is_secondary or read.is_supplementary:
        return False
    if read.is_duplicate or read.is_qcfail:
        return False
    if read.mapping_quality < mapq:
        return False
    return True


def should_count_fragment(read: pysam.AlignedSegment) -> bool:
    if read.is_paired:
        return read.is_read1
    return True


def count_bins_for_bam(
    bam_path: Path,
    chroms: Sequence[str],
    chrom_start_bin: Dict[str, int],
    bin_size: int,
    mapq: int,
    total_bins: int,
    chunk_size: int = 500000,
) -> Tuple[np.ndarray, int]:
    if not bam_path.exists():
        raise FileNotFoundError(f"BAM not found: {bam_path}")

    counts = np.zeros(int(total_bins), dtype=np.int64)
    mapped_fragments = 0

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for chrom in chroms:
            start_bin = chrom_start_bin[chrom]
            buffer = np.empty(int(chunk_size), dtype=np.int64)
            buffered = 0
            for read in bam.fetch(chrom):
                if not is_usable_read(read, mapq):
                    continue
                if not should_count_fragment(read):
                    continue
                pos = read.reference_start
                if pos < 0:
                    continue
                global_bin = start_bin + (pos // bin_size)
                if 0 <= global_bin < total_bins:
                    buffer[buffered] = global_bin
                    buffered += 1
                    mapped_fragments += 1
                    if buffered >= chunk_size:
                        counts += np.bincount(buffer[:buffered], minlength=total_bins)
                        buffered = 0
            if buffered:
                counts += np.bincount(buffer[:buffered], minlength=total_bins)

    return counts, mapped_fragments


def run_bam_count_task(task: BamCountTask) -> Tuple[str, np.ndarray, int]:
    counts, mapped = count_bins_for_bam(
        bam_path=Path(task.bam_path),
        chroms=task.chroms,
        chrom_start_bin=task.chrom_start_bin,
        bin_size=task.bin_size,
        mapq=task.mapq,
        total_bins=task.total_bins,
    )
    return task.bam_path, counts, mapped


def count_multiple_bams(
    bam_paths: Sequence[Path],
    chroms: Sequence[str],
    chrom_start_bin: Dict[str, int],
    bin_size: int,
    mapq: int,
    threads: int,
    total_bins: int,
) -> Dict[str, Tuple[np.ndarray, int]]:
    if not bam_paths:
        return {}

    workers = max(1, min(int(threads), len(bam_paths)))
    results: Dict[str, Tuple[np.ndarray, int]] = {}

    if workers == 1:
        for bam_path in bam_paths:
            key, counts, mapped = run_bam_count_task(
                BamCountTask(
                    bam_path=str(bam_path),
                    chroms=list(chroms),
                    chrom_start_bin=chrom_start_bin,
                    bin_size=int(bin_size),
                    mapq=int(mapq),
                    total_bins=int(total_bins),
                )
            )
            results[key] = (counts, mapped)
        return results

    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = [
            executor.submit(
                run_bam_count_task,
                BamCountTask(
                    bam_path=str(bam_path),
                    chroms=list(chroms),
                    chrom_start_bin=chrom_start_bin,
                    bin_size=int(bin_size),
                    mapq=int(mapq),
                    total_bins=int(total_bins),
                ),
            )
            for bam_path in bam_paths
        ]
        for future in as_completed(futures):
            key, counts, mapped = future.result()
            results[key] = (counts, mapped)
    return results


def log2_cpm(raw_counts: np.ndarray) -> np.ndarray:
    total = float(np.sum(raw_counts))
    if total <= 0:
        raise ValueError("Total fragment count is zero; cannot normalize")
    cpm = (raw_counts.astype(np.float64) / total) * 1e6
    return np.log2(cpm + 1.0)


def robust_mad(values: np.ndarray) -> float:
    if values.size == 0:
        return float("nan")
    med = float(np.median(values))
    return float(np.median(np.abs(values - med)))


def gini_coefficient(values: np.ndarray) -> float:
    x = np.asarray(values, dtype=np.float64)
    if x.size == 0:
        return float("nan")
    if np.min(x) < 0:
        x = x - np.min(x)
    total = np.sum(x)
    if total <= 0:
        return 0.0
    x = np.sort(x)
    n = x.size
    idx = np.arange(1, n + 1, dtype=np.float64)
    return float((2.0 * np.sum(idx * x) / (n * total)) - (n + 1.0) / n)


def compute_uniformity_metrics(raw_counts: np.ndarray, mapped_fragments: int, signal: np.ndarray) -> Dict[str, float]:
    usable_bins = int(raw_counts.size)
    zero_bin_fraction = float(np.mean(raw_counts == 0)) if usable_bins > 0 else float("nan")

    mean_count = float(np.mean(raw_counts)) if usable_bins > 0 else 0.0
    std_count = float(np.std(raw_counts, ddof=0)) if usable_bins > 0 else float("nan")
    bin_cv = (std_count / mean_count) if mean_count > 0 else float("inf")

    diffs = np.diff(signal)
    adj_mad = robust_mad(diffs)
    gini = gini_coefficient(raw_counts)

    return {
        "mapped_fragments": float(mapped_fragments),
        "usable_bins": float(usable_bins),
        "zero_bin_fraction": zero_bin_fraction,
        "bin_cv": float(bin_cv),
        "adjacent_diff_mad": adj_mad,
        "gini_coefficient": gini,
    }


def compute_reference_stats(ref_signals: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    ref_median = np.median(ref_signals, axis=0)
    ref_mad = np.median(np.abs(ref_signals - ref_median), axis=0)
    ref_sigma = 1.4826 * ref_mad

    positive_sigma = ref_sigma[ref_sigma > 0]
    sigma_floor = float(np.percentile(positive_sigma, 5)) if positive_sigma.size > 0 else 1e-3
    ref_sigma = np.maximum(ref_sigma, sigma_floor)
    return ref_median, ref_sigma


def safe_corr(x: np.ndarray, y: np.ndarray) -> Tuple[float, float]:
    mask = np.isfinite(x) & np.isfinite(y)
    if int(np.sum(mask)) < 3:
        return float("nan"), float("nan")
    return float(pearsonr(x[mask], y[mask])[0]), float(spearmanr(x[mask], y[mask])[0])


def compute_gc_bias_metrics(gc_fraction: np.ndarray, signal: np.ndarray) -> Dict[str, float]:
    mask = np.isfinite(gc_fraction) & np.isfinite(signal)
    if int(np.sum(mask)) < 3:
        return {
            "gc_fraction_mean": float("nan"),
            "gc_signal_pearson_r": float("nan"),
            "gc_signal_spearman_r": float("nan"),
            "gc_signal_slope": float("nan"),
        }

    gc_used = gc_fraction[mask]
    signal_used = signal[mask]
    pearson_r, spearman_r = safe_corr(gc_used, signal_used)
    slope = float(np.polyfit(gc_used, signal_used, 1)[0])
    return {
        "gc_fraction_mean": float(np.mean(gc_used)),
        "gc_signal_pearson_r": pearson_r,
        "gc_signal_spearman_r": spearman_r,
        "gc_signal_slope": slope,
    }


def tricube_kernel(distance: np.ndarray) -> np.ndarray:
    clipped = np.clip(distance, 0.0, 1.0)
    return np.power(1.0 - np.power(clipped, 3.0), 3.0)


def robust_bisquare_weights(residuals: np.ndarray, scale: float) -> np.ndarray:
    if not np.isfinite(scale) or scale <= 0:
        return np.ones_like(residuals, dtype=np.float64)
    scaled = residuals / (6.0 * scale)
    weights = np.zeros_like(scaled, dtype=np.float64)
    mask = np.abs(scaled) < 1.0
    weights[mask] = np.power(1.0 - np.power(scaled[mask], 2.0), 2.0)
    return weights


def weighted_local_linear_fit(x_window: np.ndarray, y_window: np.ndarray, weights: np.ndarray, x0: float) -> float:
    positive = weights > 0
    if int(np.sum(positive)) < 2:
        return float(np.average(y_window)) if y_window.size else float("nan")
    centered_x = x_window - x0
    design = np.column_stack([np.ones(centered_x.size, dtype=np.float64), centered_x])
    sqrt_w = np.sqrt(weights[positive])
    lhs = design[positive] * sqrt_w[:, None]
    rhs = y_window[positive] * sqrt_w
    try:
        beta, _, _, _ = np.linalg.lstsq(lhs, rhs, rcond=None)
    except np.linalg.LinAlgError:
        return float(np.average(y_window[positive], weights=weights[positive]))
    return float(beta[0])


def lowess_smooth(x: np.ndarray, y: np.ndarray, frac: float, robust_iters: int) -> np.ndarray:
    n = x.size
    if n < 3:
        return y.copy()

    order = np.argsort(x)
    x_sorted = x[order]
    y_sorted = y[order]
    fitted = np.zeros(n, dtype=np.float64)
    window = max(3, int(math.ceil(float(frac) * float(n))))
    robust = np.ones(n, dtype=np.float64)

    for iteration in range(max(0, int(robust_iters)) + 1):
        for idx in range(n):
            left = max(0, idx - window // 2)
            right = min(n, left + window)
            left = max(0, right - window)

            x_window = x_sorted[left:right]
            y_window = y_sorted[left:right]
            distances = np.abs(x_window - x_sorted[idx])
            max_distance = float(np.max(distances)) if distances.size else 0.0
            if max_distance <= 0:
                fitted[idx] = float(np.mean(y_window))
                continue

            kernel = tricube_kernel(distances / max_distance)
            weights = kernel * robust[left:right]
            fitted[idx] = weighted_local_linear_fit(x_window, y_window, weights, x_sorted[idx])

        if iteration >= int(robust_iters):
            break

        residuals = y_sorted - fitted
        scale = float(np.median(np.abs(residuals)))
        robust = robust_bisquare_weights(residuals, scale)

    restored = np.zeros(n, dtype=np.float64)
    restored[order] = fitted
    return restored


def polynomial_smooth(x: np.ndarray, y: np.ndarray, degree: int) -> np.ndarray:
    degree = max(1, int(degree))
    max_degree = min(degree, max(1, x.size - 1))
    coeff = np.polyfit(x, y, deg=max_degree)
    return np.polyval(coeff, x)


def apply_gc_correction(
    signal: np.ndarray,
    gc_fraction: np.ndarray,
    method: str,
    frac: float,
    poly_degree: int,
    min_valid_bins: int,
    robust_iters: int,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, float]]:
    corrected = signal.astype(np.float64).copy()
    trend = np.full(signal.shape, np.nan, dtype=np.float64)
    metadata = {
        "gc_correction_applied": 0.0,
        "gc_correction_valid_bins": 0.0,
    }

    method = str(method).strip().lower()
    valid = np.isfinite(signal) & np.isfinite(gc_fraction)
    valid_count = int(np.sum(valid))
    metadata["gc_correction_valid_bins"] = float(valid_count)

    if method == "none" or valid_count < max(3, int(min_valid_bins)):
        return corrected, trend, metadata

    x = gc_fraction[valid].astype(np.float64)
    y = signal[valid].astype(np.float64)

    if method == "loess":
        trend_valid = lowess_smooth(x, y, frac=max(0.05, min(float(frac), 1.0)), robust_iters=robust_iters)
    elif method == "poly2":
        trend_valid = polynomial_smooth(x, y, degree=poly_degree)
    else:
        raise ValueError(f"Unsupported GC correction method: {method}")

    center = float(np.median(y)) if y.size else 0.0
    corrected_valid = y - trend_valid + center
    corrected[valid] = corrected_valid
    trend[valid] = trend_valid
    metadata["gc_correction_applied"] = 1.0
    return corrected, trend, metadata


def qc_decision(metrics: Dict[str, float], thresholds: QcThresholds) -> Tuple[str, str]:
    # Heuristic defaults only: these thresholds are not absolute standards.
    # They should be calibrated with historical good/bad samples on this platform.
    fail_reasons: List[str] = []
    warn_reasons: List[str] = []

    checks = [
        ("mapped_fragments", metrics["mapped_fragments"], "min", thresholds.min_mapped_warn, thresholds.min_mapped_fail),
        ("zero_bin_fraction", metrics["zero_bin_fraction"], "max", thresholds.max_zero_frac_warn, thresholds.max_zero_frac_fail),
        ("bin_cv", metrics["bin_cv"], "max", thresholds.max_bin_cv_warn, thresholds.max_bin_cv_fail),
        ("adjacent_diff_mad", metrics["adjacent_diff_mad"], "max", thresholds.max_adj_mad_warn, thresholds.max_adj_mad_fail),
        ("gini_coefficient", metrics["gini_coefficient"], "max", thresholds.max_gini_warn, thresholds.max_gini_fail),
        ("pearson_r", metrics["pearson_r"], "min", thresholds.min_pearson_warn, thresholds.min_pearson_fail),
        ("spearman_r", metrics["spearman_r"], "min", thresholds.min_spearman_warn, thresholds.min_spearman_fail),
        ("median_abs_z", metrics["median_abs_z"], "max", thresholds.max_median_abs_z_warn, thresholds.max_median_abs_z_fail),
        (
            "outlier_frac_abs_z_gt_3",
            metrics["outlier_frac_abs_z_gt_3"],
            "max",
            thresholds.max_outlier3_warn,
            thresholds.max_outlier3_fail,
        ),
    ]

    for name, value, mode, warn_th, fail_th in checks:
        if not np.isfinite(value):
            fail_reasons.append(f"{name}=nan")
            continue
        if mode == "min":
            if value < fail_th:
                fail_reasons.append(f"{name}<{fail_th}")
            elif value < warn_th:
                warn_reasons.append(f"{name}<{warn_th}")
        else:
            if value > fail_th:
                fail_reasons.append(f"{name}>{fail_th}")
            elif value > warn_th:
                warn_reasons.append(f"{name}>{warn_th}")

    if fail_reasons:
        return "FAIL", ";".join(fail_reasons)
    if warn_reasons:
        return "WARN", ";".join(warn_reasons)
    return "PASS", "all_metrics_within_default_heuristic_thresholds"


def plot_profiles(
    raw_signal: np.ndarray,
    corrected_signal: np.ndarray,
    ref_median: np.ndarray,
    z_scores: np.ndarray,
    out_png: Path,
) -> None:
    x = np.arange(corrected_signal.size)
    delta = corrected_signal - ref_median

    fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)
    axes[0].plot(x, ref_median, lw=1.0, label="ref_median")
    axes[0].plot(x, raw_signal, lw=0.7, alpha=0.45, label="target_raw")
    axes[0].plot(x, corrected_signal, lw=0.9, alpha=0.95, label="target_corrected")
    axes[0].set_ylabel("log2(CPM+1)")
    axes[0].legend(loc="upper right", ncol=3)
    axes[0].set_title("Raw/corrected target vs corrected reference median profile")

    axes[1].plot(x, delta, lw=0.8, label="target-ref_median")
    axes[1].plot(x, z_scores, lw=0.6, alpha=0.7, label="z_score")
    axes[1].axhline(0, color="black", lw=0.6)
    axes[1].axhline(3, color="red", lw=0.6, ls="--")
    axes[1].axhline(-3, color="red", lw=0.6, ls="--")
    axes[1].set_xlabel("Autosome bin index")
    axes[1].set_ylabel("Delta / Z")
    axes[1].legend(loc="upper right")

    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def plot_gc_bias(gc_fraction: np.ndarray, raw_signal: np.ndarray, corrected_signal: np.ndarray, out_png: Path) -> None:
    mask_raw = np.isfinite(gc_fraction) & np.isfinite(raw_signal)
    mask_corrected = np.isfinite(gc_fraction) & np.isfinite(corrected_signal)
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 6.0), sharey=True)
    panels = [
        (axes[0], mask_raw, raw_signal, "Raw GC bias", "#2563eb"),
        (axes[1], mask_corrected, corrected_signal, "Corrected GC bias", "#059669"),
    ]
    for axis, mask, signal, title, color in panels:
        axis.scatter(gc_fraction[mask], signal[mask], s=8, alpha=0.35, color=color, edgecolors="none")
        if int(np.sum(mask)) >= 3:
            coeff = np.polyfit(gc_fraction[mask], signal[mask], 1)
            x_line = np.linspace(float(np.min(gc_fraction[mask])), float(np.max(gc_fraction[mask])), 100)
            y_line = coeff[0] * x_line + coeff[1]
            axis.plot(x_line, y_line, color="#dc2626", lw=2.0, label=f"slope={coeff[0]:.4f}")
            axis.legend(loc="best")
        axis.set_xlabel("Bin GC fraction")
        axis.set_title(title)
    axes[0].set_ylabel("log2(CPM+1)")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    setup_logger(args.log)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    target_bam = Path(args.target_bam)
    ref_bams = [Path(p) for p in args.ref_bams]
    if not ref_bams:
        raise ValueError("At least one --ref-bams path is required")

    thresholds = QcThresholds(
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

    with pysam.AlignmentFile(str(target_bam), "rb") as tb:
        chroms = infer_autosome_names(tb)
        chrom_lengths = [tb.get_reference_length(chrom) for chrom in chroms]
    layout = build_bin_layout(chroms, chrom_lengths, args.bin_size)
    bins_df = layout_to_dataframe(layout)
    gc_fraction = compute_gc_fraction_for_layout(Path(args.reference_fasta), layout)

    LOGGER.info("autosome bins: %d", len(bins_df))
    all_bams = [target_bam] + ref_bams
    count_results = count_multiple_bams(
        bam_paths=all_bams,
        chroms=layout.chroms,
        chrom_start_bin=layout.chrom_start_bin,
        bin_size=args.bin_size,
        mapq=args.mapq,
        threads=args.threads,
        total_bins=layout.total_bins,
    )

    target_counts, target_mapped = count_results[str(target_bam)]
    target_signal_raw = log2_cpm(target_counts)
    target_signal, target_gc_trend, target_gc_meta = apply_gc_correction(
        signal=target_signal_raw,
        gc_fraction=gc_fraction,
        method=args.gc_correction_method,
        frac=args.gc_correction_frac,
        poly_degree=args.gc_correction_poly_degree,
        min_valid_bins=args.gc_correction_min_valid_bins,
        robust_iters=args.gc_correction_robust_iters,
    )
    target_uniform = compute_uniformity_metrics(target_counts, target_mapped, target_signal)
    target_gc_bias_raw = compute_gc_bias_metrics(gc_fraction, target_signal_raw)
    target_gc_bias_corrected = compute_gc_bias_metrics(gc_fraction, target_signal)

    ref_rows: List[Dict[str, float]] = []
    ref_signals = []
    for ref_bam in ref_bams:
        counts, mapped = count_results[str(ref_bam)]
        signal_raw = log2_cpm(counts)
        signal, _, _ = apply_gc_correction(
            signal=signal_raw,
            gc_fraction=gc_fraction,
            method=args.gc_correction_method,
            frac=args.gc_correction_frac,
            poly_degree=args.gc_correction_poly_degree,
            min_valid_bins=args.gc_correction_min_valid_bins,
            robust_iters=args.gc_correction_robust_iters,
        )
        ref_signals.append(signal)

        ref_rows.append(
            {
                "ref_bam": str(ref_bam),
                "mapped_fragments": int(mapped),
                "usable_bins": int(len(counts)),
                "zero_bin_fraction": float(np.mean(counts == 0)),
            }
        )

    if not ref_signals:
        raise ValueError("No valid reference BAM signals generated")

    ref_matrix = np.vstack(ref_signals)
    ref_median, ref_sigma = compute_reference_stats(ref_matrix)

    pearson_r, spearman_r = safe_corr(target_signal, ref_median)
    z_scores = (target_signal - ref_median) / ref_sigma
    abs_z = np.abs(z_scores)

    ref_metrics = {
        "pearson_r": pearson_r,
        "spearman_r": spearman_r,
        "median_abs_z": float(np.median(abs_z)),
        "outlier_frac_abs_z_gt_3": float(np.mean(abs_z > 3.0)),
        "outlier_frac_abs_z_gt_5": float(np.mean(abs_z > 5.0)),
    }

    all_metrics: Dict[str, float] = {}
    all_metrics.update(target_uniform)
    all_metrics.update(ref_metrics)
    all_metrics.update(target_gc_bias_corrected)

    decision, reason = qc_decision(all_metrics, thresholds)

    profile_df = bins_df.copy()
    profile_df["gc_fraction"] = gc_fraction
    profile_df["raw_count"] = target_counts
    profile_df["raw_signal_log2cpm"] = target_signal_raw
    profile_df["gc_correction_trend"] = target_gc_trend
    profile_df["signal_log2cpm"] = target_signal
    profile_df["ref_median"] = ref_median
    profile_df["ref_sigma"] = ref_sigma
    profile_df["z_score"] = z_scores

    qc_df = pd.DataFrame(
        [
            {
                "target_bam": str(target_bam),
                "bin_size": int(args.bin_size),
                "mapped_fragments": int(target_mapped),
                "usable_bins": int(target_uniform["usable_bins"]),
                "zero_bin_fraction": target_uniform["zero_bin_fraction"],
                "bin_cv": target_uniform["bin_cv"],
                "adjacent_diff_mad": target_uniform["adjacent_diff_mad"],
                "gini_coefficient": target_uniform["gini_coefficient"],
                "pearson_r": pearson_r,
                "spearman_r": spearman_r,
                "median_abs_z": ref_metrics["median_abs_z"],
                "outlier_frac_abs_z_gt_3": ref_metrics["outlier_frac_abs_z_gt_3"],
                "outlier_frac_abs_z_gt_5": ref_metrics["outlier_frac_abs_z_gt_5"],
                "gc_fraction_mean": target_gc_bias_corrected["gc_fraction_mean"],
                "gc_signal_pearson_r": target_gc_bias_corrected["gc_signal_pearson_r"],
                "gc_signal_spearman_r": target_gc_bias_corrected["gc_signal_spearman_r"],
                "gc_signal_slope": target_gc_bias_corrected["gc_signal_slope"],
                "raw_gc_signal_pearson_r": target_gc_bias_raw["gc_signal_pearson_r"],
                "raw_gc_signal_spearman_r": target_gc_bias_raw["gc_signal_spearman_r"],
                "raw_gc_signal_slope": target_gc_bias_raw["gc_signal_slope"],
                "gc_correction_method": args.gc_correction_method,
                "gc_correction_applied": int(target_gc_meta["gc_correction_applied"]),
                "gc_correction_valid_bins": int(target_gc_meta["gc_correction_valid_bins"]),
                "qc_decision": decision,
                "qc_reason": reason,
            }
        ]
    )

    ref_df = pd.DataFrame(ref_rows)

    qc_df.to_csv(outdir / "qc_metrics.tsv", sep="\t", index=False)
    profile_df.to_csv(outdir / "target_bin_profile.tsv", sep="\t", index=False)
    ref_df.to_csv(outdir / "reference_summary.tsv", sep="\t", index=False)
    plot_profiles(target_signal_raw, target_signal, ref_median, z_scores, outdir / "target_vs_ref_profile.png")
    plot_gc_bias(gc_fraction, target_signal_raw, target_signal, outdir / "gc_bias_plot.png")

    LOGGER.info("QC done: decision=%s reason=%s", decision, reason)
    LOGGER.info("Outputs written to %s", outdir)


if __name__ == "__main__":
    main()
