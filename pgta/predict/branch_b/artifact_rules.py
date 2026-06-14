#!/biosoftware/miniconda/envs/snakemake_env/bin/python
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from pgta.predict.branch_b.common import read_bins_and_candidates, read_table, write_json, write_table
from pgta.core.logging import setup_logger


CHRX_CENTROMERE_BY_GENOME = {
    "hg19": (58527181, 61882314),
    "grch37": (58527181, 61882314),
    "hg38": (58500748, 62662843),
    "grch38": (58500748, 62662843),
}


def parse_args():
    parser = argparse.ArgumentParser(description="Apply explicit artifact rules to calibrated CNV candidate events.")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--input-bins", required=True)
    parser.add_argument("--input-candidates", required=True)
    parser.add_argument("--gender-tsv", default="")
    parser.add_argument("--output-events", required=True)
    parser.add_argument("--output-summary", required=True)
    parser.add_argument("--output-json", required=True)
    parser.add_argument("--genome-build", default="hg19")
    parser.add_argument("--par-region", action="append", default=[])
    parser.add_argument("--min-event-bins", type=int, default=3)
    parser.add_argument("--min-abs-calibrated-z", type=float, default=2.0)
    parser.add_argument("--max-chrom-fraction", type=float, default=0.35)
    parser.add_argument("--edge-bin-window", type=int, default=2)
    parser.add_argument("--max-qvalue", type=float, default=0.25)
    parser.add_argument("--keep-review", type=int, default=1)
    parser.add_argument("--high-confidence-z", type=float, default=4.0)
    parser.add_argument("--high-confidence-qvalue", type=float, default=0.05)
    parser.add_argument("--broad-support-min-abs-z", type=float, default=4.0)
    parser.add_argument("--broad-support-max-qvalue", type=float, default=0.25)
    parser.add_argument("--broad-support-min-clean-fraction", type=float, default=0.30)
    parser.add_argument("--broad-support-min-effective-bins", type=float, default=10.0)
    parser.add_argument("--broad-gain-rescue-min-abs-z", type=float, default=1.8)
    parser.add_argument("--broad-gain-rescue-min-support-fraction", type=float, default=0.35)
    parser.add_argument("--edge-review-min-priority", type=float, default=2.0)
    parser.add_argument("--ultra-pass-z", type=float, default=15.0)
    parser.add_argument("--ultra-pass-qvalue", type=float, default=0.001)
    parser.add_argument("--ultra-pass-effective-bins", type=float, default=8.0)
    parser.add_argument("--clean-review-min-support-fraction", type=float, default=0.50)
    parser.add_argument("--clean-review-max-overlap-fraction", type=float, default=0.15)
    parser.add_argument("--clean-review-max-region-risk", type=float, default=0.35)
    parser.add_argument("--focal-review-min-support-z", type=float, default=6.0)
    parser.add_argument("--focal-review-max-overlap-fraction", type=float, default=0.25)
    parser.add_argument("--focal-review-max-region-risk", type=float, default=0.20)
    parser.add_argument("--a-branch-review-min-abs-z", type=float, default=15.0)
    parser.add_argument("--a-branch-sensitive-review-min-abs-z", type=float, default=7.0)
    parser.add_argument("--a-branch-sensitive-review-min-bins", type=float, default=10.0)
    parser.add_argument("--a-branch-sensitive-review-max-high-risk-fraction", type=float, default=0.05)
    parser.add_argument("--a-branch-sensitive-review-max-region-risk", type=float, default=0.20)
    parser.add_argument("--a-branch-sensitive-review-min-same-direction-z", type=float, default=0.25)
    parser.add_argument("--a-branch-discordant-protect-min-abs-z", type=float, default=50.0)
    parser.add_argument("--branch-b-direction-min-abs-z", type=float, default=0.25)
    parser.add_argument("--recurrent-artifact-chrom", action="append", default=["chr19", "chr22"])
    parser.add_argument("--recurrent-artifact-min-same-direction-z", type=float, default=2.0)
    parser.add_argument("--recurrent-artifact-protect-min-a-abs-z", type=float, default=80.0)
    parser.add_argument("--recurrent-edge-artifact-chrom", action="append", default=["chr17"])
    parser.add_argument("--paired-event-rescue-min-a-abs-z", type=float, default=8.0)
    parser.add_argument("--paired-event-rescue-mate-min-a-abs-z", type=float, default=50.0)
    parser.add_argument("--paired-event-rescue-min-bins", type=int, default=5)
    parser.add_argument("--narrow-boundary-artifact-max-bins", type=int, default=15)
    parser.add_argument("--narrow-boundary-artifact-max-available-chrom-fraction", type=float, default=0.08)
    parser.add_argument("--narrow-boundary-artifact-protect-min-a-abs-z", type=float, default=50.0)
    parser.add_argument("--sca-xy-xgain-max-bam-x-relative", type=float, default=0.80)
    parser.add_argument("--sca-xy-xgain-focal-edge-max-bins", type=int, default=20)
    parser.add_argument("--cnvseq-reportable-min-bp", type=int, default=2_000_000)
    parser.add_argument("--cnvseq-review-min-bp", type=int, default=1_000_000)
    parser.add_argument("--cnvseq-large-event-min-bp", type=int, default=10_000_000)
    parser.add_argument("--cnvseq-boundary-max-abs-z", type=float, default=4.0)
    parser.add_argument("--cnvseq-whole-chrom-available-fraction", type=float, default=0.90)
    parser.add_argument("--log", default="")
    return parser.parse_args()


def parse_gender_tsv(path_value):
    info = parse_gender_info(path_value)
    return str(info.get("sex_call", "")).strip().upper()


def parse_gender_info(path_value):
    if not path_value:
        return {}
    path = Path(path_value)
    if not path.exists():
        return {}
    df = read_table(path, empty_ok=True)
    if df.empty or "sex_call" not in df.columns:
        return {}
    return df.iloc[0].to_dict()


def parse_par_regions(region_specs):
    parsed = {}
    for spec in region_specs:
        chrom_part, sep, range_part = str(spec).partition(":")
        start_part, sep2, end_part = range_part.partition("-")
        if not sep or not sep2:
            raise ValueError(f"Invalid --par-region specification: {spec!r}")
        chrom = chrom_part.strip()
        start = int(start_part)
        end = int(end_part)
        parsed.setdefault(chrom, []).append((start, end))
    for chrom in parsed:
        parsed[chrom] = sorted(parsed[chrom])
    return parsed


def interval_overlap(start, end, left, right):
    return max(0, min(end, right) - max(start, left))


def compute_par_overlap(chrom, start, end, par_regions):
    overlap_bp = 0
    for left, right in par_regions.get(chrom, []):
        overlap_bp += interval_overlap(start, end, left, right)
    event_length = max(int(end) - int(start), 1)
    return overlap_bp, overlap_bp / float(event_length)


def sca_region_class(chrom, start, end, par_fraction, genome_build="hg19"):
    if chrom == "chrY":
        return "sca_y"
    if chrom != "chrX":
        return ""
    if par_fraction >= 0.95:
        return "sca_x_par_only"
    if par_fraction >= 0.05:
        return "sca_x_mixed_par_nonpar"

    centromere_start, centromere_end = CHRX_CENTROMERE_BY_GENOME.get(
        str(genome_build).lower(),
        CHRX_CENTROMERE_BY_GENOME["hg19"],
    )
    if int(end) < centromere_start:
        return "sca_x_nonpar_p"
    if int(start) >= centromere_start and int(end) > centromere_end:
        return "sca_x_nonpar_q"
    return "sca_x_nonpar_centromere_crossing"


def safe_float(value, default=np.nan):
    try:
        return float(value)
    except (TypeError, ValueError):
        return float(default)


def state_direction(state):
    normalized = str(state).strip().lower()
    if normalized == "gain":
        return 1
    if normalized == "loss":
        return -1
    return 0


def has_branch_b_direction_discordance(state, evidence_values, min_abs_z=0.25):
    expected_direction = state_direction(state)
    if expected_direction == 0:
        return False
    informative_values = [
        float(value)
        for value in evidence_values
        if np.isfinite(value) and abs(float(value)) >= min_abs_z
    ]
    return len(informative_values) >= 2 and all(value * expected_direction < 0 for value in informative_values)


def append_reason(existing, reason):
    values = [item for item in str(existing or "").split(";") if item]
    if reason not in values:
        values.append(reason)
    return ";".join(values)


def remove_reason(existing, reason):
    values = [item for item in str(existing or "").split(";") if item and item != reason]
    return ";".join(values)


def append_flag(existing, flag):
    values = [item for item in str(existing or "").split(",") if item]
    if flag not in values:
        values.append(flag)
    return ",".join(values)


def recurrent_artifact_chroms(values):
    chroms = []
    for value in values or []:
        chroms.extend([item.strip() for item in str(value).split(",") if item.strip()])
    return set(chroms)


def summarize_cnvseq_event_context(row, bins_df):
    chrom = str(getattr(row, "chrom", ""))
    chrom_df = bins_df[bins_df["chrom"].astype(str).eq(chrom)].copy()
    if chrom_df.empty:
        return {
            "cnvseq_available_bin_count": 0,
            "cnvseq_chrom_available_bin_count": 0,
            "cnvseq_available_chrom_fraction": np.nan,
            "cnvseq_gap_centromere_bin_fraction": 0.0,
            "cnvseq_crosses_gap_or_centromere": 0,
            "cnvseq_region_class_transition": 0,
        }

    start_bin = int(getattr(row, "start_bin", 0))
    end_bin = int(getattr(row, "end_bin", start_bin))
    event_df = chrom_df[chrom_df["bin_index"].between(start_bin, end_bin)].copy()
    if event_df.empty:
        return {
            "cnvseq_available_bin_count": 0,
            "cnvseq_chrom_available_bin_count": int(len(chrom_df)),
            "cnvseq_available_chrom_fraction": 0.0,
            "cnvseq_gap_centromere_bin_fraction": 0.0,
            "cnvseq_crosses_gap_or_centromere": 0,
            "cnvseq_region_class_transition": 0,
        }

    chrom_mask_label = chrom_df.get("mask_label", pd.Series("pass", index=chrom_df.index)).fillna("pass").astype(str)
    event_mask_label = event_df.get("mask_label", pd.Series("pass", index=event_df.index)).fillna("pass").astype(str)
    chrom_available = int((chrom_mask_label != "hard").sum())
    event_available = int((event_mask_label != "hard").sum())
    available_fraction = event_available / float(max(chrom_available, 1))

    gap_values = pd.to_numeric(
        event_df.get("gap_centromere_telomere_overlap_fraction", pd.Series(0.0, index=event_df.index)),
        errors="coerce",
    ).fillna(0.0)
    gap_bin_fraction = float((gap_values > 0.0).mean()) if len(gap_values) else 0.0

    risk_classes = event_df.get("region_risk_class", pd.Series("", index=event_df.index)).fillna("").astype(str)
    region_class_transition = int(len([item for item in risk_classes.unique() if item]) > 1)

    return {
        "cnvseq_available_bin_count": event_available,
        "cnvseq_chrom_available_bin_count": chrom_available,
        "cnvseq_available_chrom_fraction": available_fraction,
        "cnvseq_gap_centromere_bin_fraction": gap_bin_fraction,
        "cnvseq_crosses_gap_or_centromere": int(gap_bin_fraction > 0.0),
        "cnvseq_region_class_transition": region_class_transition,
    }


def classify_cnvseq_report_tier(event_length_bp, cnvseq_available_chrom_fraction, args):
    if cnvseq_available_chrom_fraction >= float(getattr(args, "cnvseq_whole_chrom_available_fraction", 0.90)):
        return "whole_chromosome", 1
    if event_length_bp >= int(getattr(args, "cnvseq_reportable_min_bp", 2_000_000)):
        return "reportable", 1
    if event_length_bp >= int(getattr(args, "cnvseq_review_min_bp", 1_000_000)):
        return "review_1_2mb", 0
    return "subreportable_lt1mb", 0


def classify_event(row, chrom_bin_count, args, sex_call, par_regions):
    flags = []
    explanations = []
    retain_reasons = []
    downgrade_reasons = []
    filter_reasons = []
    hard_artifact = False
    review_only = False
    chrom = str(row.chrom)
    state = str(getattr(row, "state", "")).strip().lower()
    caller = str(getattr(row, "caller", ""))
    a_candidate_id = str(getattr(row, "a_candidate_id", "") or "").strip()
    a_abs_z = abs(safe_float(getattr(row, "a_abs_zscore", np.nan), default=0.0))
    a_branch_review_min_abs_z = safe_float(getattr(args, "a_branch_review_min_abs_z", 15.0), default=15.0)
    preserve_a_branch_primary_signal = (
        (caller == "wisecondorx_a_branch" or bool(a_candidate_id))
        and a_abs_z >= a_branch_review_min_abs_z
    )
    a_branch_boundary_protect_min_abs_z = safe_float(
        getattr(args, "a_branch_discordant_protect_min_abs_z", 50.0),
        default=50.0,
    )
    preserve_a_branch_boundary_signal = (
        (caller == "wisecondorx_a_branch" or bool(a_candidate_id))
        and a_abs_z >= a_branch_boundary_protect_min_abs_z
    )
    is_sex_chrom = chrom in {"chrX", "chrY"}
    chrom_fraction = float(row.n_bins) / max(chrom_bin_count, 1)
    overlap_bp, par_fraction = compute_par_overlap(chrom, int(row.start), int(row.end), par_regions)
    sca_class = sca_region_class(
        chrom,
        int(row.start),
        int(row.end),
        par_fraction,
        genome_build=getattr(args, "genome_build", "hg19"),
    )
    touches_chrom_edge = (
        int(row.start_bin) <= args.edge_bin_window
        or (chrom_bin_count - int(row.end_bin) - 1) <= args.edge_bin_window
    )
    signed_mean_z = safe_float(getattr(row, "calibrated_mean_z", np.nan))
    signed_median_z = safe_float(getattr(row, "calibrated_median_z", np.nan))
    signed_adjusted_event_z = safe_float(getattr(row, "event_corr_adjusted_z", np.nan))
    signed_segment_mean_z = safe_float(getattr(row, "segment_mean_robust_z", np.nan), default=np.nan)
    signed_segment_median_z = safe_float(getattr(row, "segment_median_robust_z", np.nan), default=np.nan)
    mean_z = abs(signed_mean_z)
    median_z = abs(signed_median_z)
    adjusted_event_z = abs(signed_adjusted_event_z)
    max_calibrated_z = max(mean_z, median_z, adjusted_event_z)
    internal_support_z = max(
        max_calibrated_z,
        abs(signed_segment_mean_z) if np.isfinite(signed_segment_mean_z) else 0.0,
        abs(signed_segment_median_z) if np.isfinite(signed_segment_median_z) else 0.0,
        abs(safe_float(getattr(row, "segment_abs_max_robust_z", np.nan), default=0.0)),
    )
    empirical_q = safe_float(getattr(row, "empirical_qvalue", np.nan))
    clean_fraction = safe_float(getattr(row, "clean_bin_fraction", np.nan), default=0.0)
    high_fraction = safe_float(getattr(row, "high_risk_bin_fraction", np.nan), default=0.0)
    moderate_fraction = safe_float(getattr(row, "moderate_risk_bin_fraction", np.nan), default=np.nan)
    if not np.isfinite(moderate_fraction):
        moderate_fraction = max(0.0, 1.0 - clean_fraction - high_fraction)
    effective_bin_count = safe_float(getattr(row, "effective_bin_count", np.nan), default=float(row.n_bins))
    region_risk_score_mean = safe_float(getattr(row, "region_risk_score_mean", np.nan), default=0.0)
    region_risk_score_max = safe_float(getattr(row, "region_risk_score_max", np.nan), default=0.0)
    event_length_bp = max(int(row.end) - int(row.start), 0)
    cnvseq_available_chrom_fraction = safe_float(
        getattr(row, "cnvseq_available_chrom_fraction", np.nan), default=chrom_fraction
    )
    cnvseq_report_tier, cnvseq_reportable = classify_cnvseq_report_tier(
        event_length_bp,
        cnvseq_available_chrom_fraction,
        args,
    )
    cnvseq_gap_centromere_bin_fraction = safe_float(
        getattr(row, "cnvseq_gap_centromere_bin_fraction", np.nan), default=0.0
    )
    cnvseq_crosses_gap_or_centromere = int(
        safe_float(getattr(row, "cnvseq_crosses_gap_or_centromere", 0), default=0.0) > 0.0
    )
    cnvseq_region_class_transition = int(
        safe_float(getattr(row, "cnvseq_region_class_transition", 0), default=0.0) > 0.0
    )
    high_risk_boundary_crossing = int(
        safe_float(getattr(row, "high_risk_boundary_crossing", 0), default=0.0) > 0.0
    )
    cnvseq_boundary_like_event = bool(cnvseq_crosses_gap_or_centromere or high_risk_boundary_crossing)
    cnvseq_segment_level_z = max(mean_z, median_z, adjusted_event_z)
    weighted_non_high_support_fraction = clean_fraction + (0.5 * max(moderate_fraction, 0.0))
    same_direction_branch_b_z = max(
        [
            value * state_direction(state)
            for value in [
                signed_mean_z,
                signed_median_z,
                signed_adjusted_event_z,
                signed_segment_mean_z,
                signed_segment_median_z,
            ]
            if np.isfinite(value)
        ]
        or [np.nan]
    )
    preserve_a_branch_sensitive_signal = (
        (caller == "wisecondorx_a_branch" or bool(a_candidate_id))
        and not is_sex_chrom
        and a_abs_z >= float(getattr(args, "a_branch_sensitive_review_min_abs_z", 7.0))
        and effective_bin_count >= float(getattr(args, "a_branch_sensitive_review_min_bins", 10.0))
        and high_fraction <= float(getattr(args, "a_branch_sensitive_review_max_high_risk_fraction", 0.05))
        and region_risk_score_max <= float(getattr(args, "a_branch_sensitive_review_max_region_risk", 0.20))
        and not cnvseq_boundary_like_event
        and (
            not np.isfinite(same_direction_branch_b_z)
            or same_direction_branch_b_z
            >= float(getattr(args, "a_branch_sensitive_review_min_same_direction_z", 0.25))
        )
    )
    preserve_broad_internal_signal = (
        chrom_fraction > args.max_chrom_fraction
        and internal_support_z >= args.broad_support_min_abs_z
        and effective_bin_count >= args.broad_support_min_effective_bins
        and clean_fraction >= args.broad_support_min_clean_fraction
        and (not np.isfinite(empirical_q) or empirical_q <= args.broad_support_max_qvalue)
    )
    preserve_broad_gain_signal = (
        chrom_fraction > args.max_chrom_fraction
        and state == "gain"
        and not is_sex_chrom
        and caller in {"chromosome_dosage_detector", "raw_chromosome_dosage_detector"}
        and internal_support_z >= args.broad_gain_rescue_min_abs_z
        and effective_bin_count >= args.broad_support_min_effective_bins
        and weighted_non_high_support_fraction >= args.broad_gain_rescue_min_support_fraction
        and high_fraction < 0.50
        and (not np.isfinite(empirical_q) or empirical_q <= args.broad_support_max_qvalue)
    )
    overlap_metrics = {
        "xtr_overlap": safe_float(getattr(row, "xtr_overlap_fraction", np.nan), default=0.0),
        "sex_homology_overlap": safe_float(getattr(row, "sex_homology_overlap_fraction", np.nan), default=0.0),
        "segmental_duplication_overlap": safe_float(
            getattr(row, "segmental_duplication_overlap_fraction", np.nan), default=0.0
        ),
        "low_mappability_overlap": safe_float(getattr(row, "low_mappability_overlap_fraction", np.nan), default=0.0),
        "gap_centromere_telomere_overlap": safe_float(
            getattr(row, "gap_centromere_telomere_overlap_fraction", np.nan), default=0.0
        ),
        "repeat_rich_overlap": safe_float(getattr(row, "repeat_rich_overlap_fraction", np.nan), default=0.0),
        "blacklist_overlap": safe_float(getattr(row, "blacklist_overlap_fraction", np.nan), default=0.0),
        "ambiguous_alignment_region": safe_float(
            getattr(row, "ambiguous_alignment_overlap_fraction", np.nan), default=0.0
        ),
    }
    max_overlap_fraction = max(overlap_metrics.values()) if overlap_metrics else 0.0
    preserve_clean_internal_signal = (
        internal_support_z >= args.high_confidence_z
        and weighted_non_high_support_fraction >= args.clean_review_min_support_fraction
        and high_fraction < 0.20
        and max_overlap_fraction <= args.clean_review_max_overlap_fraction
        and region_risk_score_max <= args.clean_review_max_region_risk
        and (not np.isfinite(empirical_q) or empirical_q <= max(args.max_qvalue, args.high_confidence_qvalue))
    )
    preserve_focal_low_risk_internal_signal = (
        chrom_fraction <= min(args.max_chrom_fraction, 0.10)
        and internal_support_z >= args.focal_review_min_support_z
        and weighted_non_high_support_fraction >= args.clean_review_min_support_fraction
        and high_fraction <= 0.05
        and max_overlap_fraction <= args.focal_review_max_overlap_fraction
        and region_risk_score_max <= args.focal_review_max_region_risk
        and overlap_metrics["xtr_overlap"] <= 0.0
        and overlap_metrics["sex_homology_overlap"] <= 0.0
        and overlap_metrics["low_mappability_overlap"] <= 0.0
        and overlap_metrics["gap_centromere_telomere_overlap"] <= 0.0
        and overlap_metrics["repeat_rich_overlap"] <= 0.0
        and overlap_metrics["blacklist_overlap"] <= 0.0
        and overlap_metrics["ambiguous_alignment_region"] <= 0.0
        and (not np.isfinite(empirical_q) or empirical_q <= max(args.max_qvalue, args.high_confidence_qvalue))
    )
    same_direction_event_z = signed_adjusted_event_z * state_direction(state)
    weak_same_direction_support = (
        not np.isfinite(same_direction_event_z)
        or same_direction_event_z < float(getattr(args, "recurrent_artifact_min_same_direction_z", 2.0))
    )
    recurrent_weak_support = (
        max_calibrated_z < args.min_abs_calibrated_z
        or (np.isfinite(empirical_q) and empirical_q > args.max_qvalue)
        or weak_same_direction_support
    )

    if int(row.n_bins) < args.min_event_bins:
        flags.append("too_few_bins")
        explanations.append("Event spans too few bins for a stable segment.")
        filter_reasons.append("bin_count_below_minimum")
        hard_artifact = True
    if max_calibrated_z < args.min_abs_calibrated_z:
        flags.append("low_calibrated_signal")
        explanations.append("Calibrated signal amplitude is below the minimum support threshold.")
        if preserve_a_branch_primary_signal:
            downgrade_reasons.append("a_branch_strong_evidence_preserved_for_review")
            review_only = True
        elif preserve_a_branch_sensitive_signal:
            downgrade_reasons.append("a_branch_sensitive_evidence_preserved_for_review")
            review_only = True
        elif preserve_broad_gain_signal:
            downgrade_reasons.append("broad_gain_preserved_by_raw_or_chromosome_support")
            review_only = True
        else:
            filter_reasons.append("signal_support_below_minimum")
            hard_artifact = True
    if chrom_fraction > args.max_chrom_fraction:
        flags.append("broad_chrom_fraction")
        if is_sex_chrom and sex_call in {"XX", "XY"}:
            explanations.append("Large sex-chromosome event requires sex-aware review.")
            downgrade_reasons.append("broad_sex_chromosome_event")
            review_only = True
        elif preserve_a_branch_primary_signal:
            explanations.append("Broad event is preserved for review because Branch A has very strong WisecondorX support.")
            downgrade_reasons.append("a_branch_strong_evidence_preserved_for_review")
            review_only = True
        elif preserve_a_branch_sensitive_signal:
            explanations.append(
                "Broad or small-chromosome event is preserved for review because Branch A has "
                "same-direction low-risk WisecondorX support."
            )
            downgrade_reasons.append("a_branch_sensitive_evidence_preserved_for_review")
            review_only = True
        elif preserve_broad_internal_signal or preserve_broad_gain_signal:
            explanations.append("Broad event is preserved for review because Branch B shows high-confidence chromosome-scale support.")
            if preserve_broad_gain_signal and not preserve_broad_internal_signal:
                downgrade_reasons.append("broad_gain_preserved_by_raw_or_chromosome_support")
            else:
                downgrade_reasons.append("broad_event_preserved_by_internal_support")
            review_only = True
        else:
            explanations.append("Event spans too much of one chromosome and is likely technical.")
            filter_reasons.append("chromosome_fraction_too_large")
            hard_artifact = True
    if touches_chrom_edge:
        flags.append("edge_event")
        explanations.append("Segment touches chromosome-edge bins and should be reviewed.")
        downgrade_reasons.append("chromosome_edge_contact")
        review_only = True
    if np.isfinite(empirical_q) and empirical_q > args.max_qvalue:
        flags.append("weak_empirical_support")
        explanations.append("Empirical null calibration indicates weak support.")
        downgrade_reasons.append("empirical_qvalue_above_threshold")
        review_only = True

    if is_sex_chrom:
        flags.append("sca_event")
        if sca_class:
            flags.append(sca_class)
            if sca_class == "sca_x_mixed_par_nonpar":
                downgrade_reasons.append("sca_mixed_par_nonpar_review")
            elif sca_class == "sca_x_par_only":
                downgrade_reasons.append("sca_par_only_review")
            elif sca_class == "sca_y":
                downgrade_reasons.append("sca_y_review")
            else:
                downgrade_reasons.append("sca_nonpar_review")
            review_only = True
        if 0.0 < par_fraction < 0.05 and sca_class not in {"sca_x_par_only", "sca_x_mixed_par_nonpar"}:
            flags.append("sca_par_boundary_overlap")
            downgrade_reasons.append("sca_par_boundary_review")
            review_only = True
        flags.append("sex_chromosome_event")
        explanations.append("Sex-chromosome event needs sex and PAR-aware interpretation.")
    if overlap_bp > 0:
        flags.append("par_overlap")
        explanations.append(f"Event overlaps a PAR interval ({par_fraction:.1%} of event span).")
        downgrade_reasons.append("par_overlap")
    if 0.0 < par_fraction < 1.0:
        flags.append("mixed_par_nonpar")
        explanations.append("Event crosses PAR and non-PAR sequence.")
        downgrade_reasons.append("mixed_par_nonpar_boundary")
        review_only = True
    if sex_call == "XX" and chrom == "chrY":
        flags.append("xx_chrY_unexpected_signal")
        explanations.append("chrY signal in an XX-routed sample is treated as artifact.")
        filter_reasons.append("xx_sample_with_chrY_signal")
        hard_artifact = True
    if sex_call == "XY" and chrom == "chrY":
        explanations.append("chrY event in XY-routed sample remains review-only by default.")
        downgrade_reasons.append("xy_chrY_default_review")
        review_only = True

    if overlap_metrics["xtr_overlap"] > 0.0:
        flags.append("xtr_overlap")
        explanations.append(f"Event overlaps XTR sequence ({overlap_metrics['xtr_overlap']:.1%}).")
        downgrade_reasons.append("xtr_overlap")
        review_only = True
    if overlap_metrics["sex_homology_overlap"] >= 0.10:
        flags.append("sex_homology_overlap")
        explanations.append("Event overlaps sex-chromosome homology sequence with elevated mapping ambiguity.")
        downgrade_reasons.append("sex_homology_overlap")
        review_only = True
    if overlap_metrics["segmental_duplication_overlap"] >= 0.25:
        flags.append("segmental_duplication_overlap")
        explanations.append("Event overlaps segmental duplication sequence.")
        if clean_fraction < args.clean_review_min_support_fraction:
            filter_reasons.append("segmental_duplication_overlap_with_limited_clean_support")
            hard_artifact = True
        else:
            downgrade_reasons.append("segmental_duplication_overlap")
            review_only = True
    if overlap_metrics["low_mappability_overlap"] >= 0.25:
        flags.append("low_mappability_overlap")
        explanations.append("Event overlaps low-mappability bins.")
        downgrade_reasons.append("low_mappability_overlap")
        review_only = True
    if overlap_metrics["repeat_rich_overlap"] >= 0.25:
        flags.append("repeat_rich_overlap")
        explanations.append("Event overlaps repeat-rich sequence.")
        downgrade_reasons.append("repeat_rich_overlap")
        review_only = True
    if overlap_metrics["gap_centromere_telomere_overlap"] > 0.0:
        flags.append("gap_centromere_telomere_overlap")
        explanations.append("Event overlaps gap / centromere / telomere-adjacent sequence.")
        downgrade_reasons.append("gap_centromere_telomere_overlap")
        review_only = True
    if overlap_metrics["blacklist_overlap"] > 0.0:
        flags.append("blacklist_overlap")
        explanations.append(
            "Event overlaps blacklisted bins; WisecondorX treats these bins as masked evidence, "
            "so overlap is review evidence rather than an event-level hard artifact."
        )
        downgrade_reasons.append("blacklist_overlap")
        review_only = True
    if overlap_metrics["ambiguous_alignment_region"] >= 0.10:
        flags.append("ambiguous_alignment_region")
        explanations.append("Event overlaps a high-risk ambiguous-alignment region.")
        downgrade_reasons.append("ambiguous_alignment_region")
        if max_calibrated_z < args.high_confidence_z:
            filter_reasons.append("ambiguous_alignment_with_non_high_confidence_signal")
            hard_artifact = True
        else:
            review_only = True
    if high_fraction >= 0.50:
        flags.append("high_risk_region_burden")
        explanations.append(f"High-risk bins cover {high_fraction:.1%} of the event.")
        downgrade_reasons.append("high_risk_bin_fraction")
        review_only = True
    if clean_fraction <= 0.20 and high_fraction >= 0.50 and max_calibrated_z < args.high_confidence_z:
        flags.append("clean_support_low")
        explanations.append("Event has limited clean-bin support relative to high-risk sequence burden.")
        filter_reasons.append("clean_bin_support_too_low")
        hard_artifact = True
    if high_risk_boundary_crossing == 1:
        flags.append("high_risk_boundary_crossing")
        explanations.append("Event crosses a clean/high-risk boundary and should be manually reviewed.")
        downgrade_reasons.append("high_risk_boundary_crossing")
        review_only = True
    if (
        high_risk_boundary_crossing == 1
        and (
            int(row.n_bins) <= int(getattr(args, "narrow_boundary_artifact_max_bins", 15))
            or cnvseq_available_chrom_fraction
            <= float(getattr(args, "narrow_boundary_artifact_max_available_chrom_fraction", 0.08))
        )
        and a_abs_z < float(getattr(args, "narrow_boundary_artifact_protect_min_a_abs_z", 50.0))
    ):
        flags.append("narrow_high_risk_boundary_artifact")
        explanations.append(
            "Small or locally sparse event crosses a high-risk sequence boundary without strong Branch A protection."
        )
        filter_reasons.append("narrow_high_risk_boundary_without_strong_a_branch_support")
        hard_artifact = True
    if (
        chrom in recurrent_artifact_chroms(getattr(args, "recurrent_artifact_chrom", []))
        and state == "gain"
        and recurrent_weak_support
        and a_abs_z < float(getattr(args, "recurrent_artifact_protect_min_a_abs_z", 80.0))
    ):
        flags.append("recurrent_artifact_region")
        explanations.append(
            "Event falls in a recurrent validation artifact region and lacks same-direction Branch B support."
        )
        filter_reasons.append("recurrent_artifact_region_weak_support")
        hard_artifact = True
    if (
        not is_sex_chrom
        and touches_chrom_edge
        and high_risk_boundary_crossing == 1
        and cnvseq_segment_level_z < args.high_confidence_z
        and not preserve_a_branch_boundary_signal
    ):
        flags.append("edge_boundary_weak_a_branch")
        explanations.append(
            "Edge-adjacent boundary event lacks strong Branch A support and is treated as a boundary artifact."
        )
        filter_reasons.append("edge_boundary_without_strong_a_branch_support")
        hard_artifact = True
    if (
        not is_sex_chrom
        and event_length_bp >= int(getattr(args, "cnvseq_large_event_min_bp", 10_000_000))
        and cnvseq_boundary_like_event
        and cnvseq_available_chrom_fraction < float(getattr(args, "cnvseq_whole_chrom_available_fraction", 0.90))
        and cnvseq_segment_level_z < float(getattr(args, "cnvseq_boundary_max_abs_z", args.high_confidence_z))
        and not preserve_a_branch_boundary_signal
    ):
        flags.append("cnvseq_subchrom_boundary_event")
        explanations.append(
            "CNVseq-style review marks this as a large subchromosomal boundary event without stable segment-level support."
        )
        filter_reasons.append("cnvseq_subchrom_boundary_weak_support")
        hard_artifact = True

    branch_b_direction_discordant = has_branch_b_direction_discordance(
        state,
        [
            signed_mean_z,
            signed_median_z,
            signed_adjusted_event_z,
            signed_segment_mean_z,
            signed_segment_median_z,
        ],
        min_abs_z=float(getattr(args, "branch_b_direction_min_abs_z", 0.25)),
    )
    a_branch_discordant_protect_min_abs_z = safe_float(
        getattr(args, "a_branch_discordant_protect_min_abs_z", 50.0),
        default=50.0,
    )
    if branch_b_direction_discordant and a_abs_z < a_branch_discordant_protect_min_abs_z:
        flags.append("branch_b_direction_discordant")
        explanations.append("Branch B signed evidence is consistently opposite to the candidate state.")
        filter_reasons.append("branch_b_direction_discordant_with_candidate_state")
        hard_artifact = True
    if (
        chrom in recurrent_artifact_chroms(getattr(args, "recurrent_edge_artifact_chrom", []))
        and touches_chrom_edge
        and max_calibrated_z < args.min_abs_calibrated_z
        and np.isfinite(empirical_q)
        and empirical_q > args.high_confidence_qvalue
        and weak_same_direction_support
        and a_abs_z < float(getattr(args, "recurrent_artifact_protect_min_a_abs_z", 80.0))
    ):
        flags.append("recurrent_edge_lowcal_artifact")
        explanations.append(
            "Edge-adjacent recurrent artifact region has weak calibrated and same-direction support."
        )
        filter_reasons.append("recurrent_edge_lowcal_weak_support")
        hard_artifact = True

    risk_penalty = (
        0.60 * high_fraction
        + 0.25 * max(overlap_metrics.values())
        + 0.20 * region_risk_score_mean
        + 0.10 * region_risk_score_max
    )
    base_priority = float(max_calibrated_z * np.log1p(max(effective_bin_count, 1.0)))
    priority_score = base_priority * max(clean_fraction, 0.10) * max(0.10, 1.0 - risk_penalty)
    if np.isfinite(empirical_q):
        priority_score *= max(0.0, 1.0 - empirical_q)
    if cnvseq_boundary_like_event:
        priority_score *= max(0.25, 1.0 - (0.50 * min(cnvseq_gap_centromere_bin_fraction, 1.0)))
    if preserve_broad_gain_signal:
        q_confidence = max(0.0, 1.0 - min(empirical_q, args.broad_support_max_qvalue)) if np.isfinite(empirical_q) else 1.0
        support_fraction = np.sqrt(min(max(weighted_non_high_support_fraction, 0.0), 1.0))
        high_risk_adjustment = max(0.50, 1.0 - high_fraction)
        broad_gain_level_z = max(
            max_calibrated_z,
            abs(safe_float(getattr(row, "segment_mean_robust_z", np.nan), default=0.0)),
            abs(safe_float(getattr(row, "segment_median_robust_z", np.nan), default=0.0)),
        )
        broad_gain_priority_floor = (
            broad_gain_level_z
            * np.log1p(max(effective_bin_count, 1.0))
            * support_fraction
            * q_confidence
            * high_risk_adjustment
            * 1.30
        )
        priority_score = max(priority_score, broad_gain_priority_floor)

    edge_only_review = review_only and set(flags) == {"edge_event"}
    if edge_only_review and priority_score < args.edge_review_min_priority:
        explanations.append("Edge-only event priority is too low to keep for manual review.")
        filter_reasons.append("edge_event_priority_below_keep_threshold")
        hard_artifact = True
        review_only = False

    ultra_pass = (
        max_calibrated_z >= args.ultra_pass_z
        and effective_bin_count >= args.ultra_pass_effective_bins
        and (not np.isfinite(empirical_q) or empirical_q <= args.ultra_pass_qvalue)
    )

    if hard_artifact:
        artifact_status = "artifact"
        keep_event = 0
        technical_confidence = "rejected"
        report_class = "technical_artifact"
    elif review_only or flags:
        artifact_status = "review"
        keep_event = int(bool(args.keep_review))
        technical_confidence = "moderate" if max_calibrated_z >= args.high_confidence_z else "low"
        report_class = "candidate_review" if keep_event else "candidate_suppressed"
    elif ultra_pass:
        artifact_status = "pass"
        keep_event = 1
        technical_confidence = "high" if (
            max_calibrated_z >= args.high_confidence_z
            and (not np.isfinite(empirical_q) or empirical_q <= args.high_confidence_qvalue)
        ) else "moderate"
        report_class = "candidate_pass"
        explanations.append("No explicit artifact rule was triggered.")
        retain_reasons.append("clean_support_and_statistical_support")
    else:
        flags.append("clean_event_below_ultra_pass")
        if preserve_clean_internal_signal or preserve_focal_low_risk_internal_signal:
            explanations.append("Low-risk event retains strong internal support but misses the ultra-pass gate, so it is kept for review.")
            if preserve_focal_low_risk_internal_signal and not preserve_clean_internal_signal:
                downgrade_reasons.append("focal_low_risk_event_preserved_below_ultra_pass")
            else:
                downgrade_reasons.append("low_risk_event_preserved_below_ultra_pass")
            artifact_status = "review"
            keep_event = int(bool(args.keep_review))
            technical_confidence = "moderate" if max_calibrated_z >= args.high_confidence_z else "low"
            report_class = "candidate_review" if keep_event else "candidate_suppressed"
        else:
            explanations.append("Clean event did not meet the ultra-pass gate and is suppressed to control false positives.")
            filter_reasons.append("clean_event_below_ultra_pass_gate")
            artifact_status = "artifact"
            keep_event = 0
            technical_confidence = "rejected"
            report_class = "technical_artifact"

    if keep_event and not cnvseq_reportable:
        flags.append(cnvseq_report_tier)
        explanations.append(
            "CNVseq-style reporting keeps this event below the routine 2 Mb reportable threshold as review-only evidence."
        )
        downgrade_reasons.append("cnvseq_subreportable_size_review")
        artifact_status = "review"
        technical_confidence = "moderate" if technical_confidence == "high" else technical_confidence
        report_class = "candidate_review"
        manual_review_recommended = 1
    else:
        manual_review_recommended = int(
            review_only or "manual" in " ".join(explanations).lower() or bool(downgrade_reasons)
        )

    biological_context = sca_class if is_sex_chrom and sca_class else ("sex_chromosome" if is_sex_chrom else "autosome")
    return {
        "artifact_status": artifact_status,
        "keep_event": keep_event,
        "artifact_flags": ",".join(flags),
        "artifact_explanations": " ".join(explanations),
        "technical_confidence": technical_confidence,
        "report_class": report_class,
        "priority_score": priority_score,
        "priority_delta": priority_score - base_priority,
        "retain_reason": ";".join(dict.fromkeys(retain_reasons)),
        "downgrade_reason": ";".join(dict.fromkeys(downgrade_reasons)),
        "filter_reason": ";".join(dict.fromkeys(filter_reasons)),
        "manual_review_recommended": manual_review_recommended,
        "biological_context": biological_context,
        "cnvseq_report_tier": cnvseq_report_tier,
        "cnvseq_reportable": int(cnvseq_reportable),
        "event_length_bp": int(event_length_bp),
    }


def apply_paired_chromosome_event_rescue(events_df, args):
    if events_df.empty:
        return events_df
    frame = events_df.copy()
    required = {
        "sample_id",
        "chrom",
        "state",
        "caller",
        "a_abs_zscore",
        "n_bins",
        "filter_reason",
        "keep_event",
    }
    if not required.issubset(frame.columns):
        return frame

    a_abs = pd.to_numeric(frame["a_abs_zscore"], errors="coerce").fillna(0.0)
    n_bins = pd.to_numeric(frame["n_bins"], errors="coerce").fillna(0.0)
    filter_text = frame["filter_reason"].fillna("").astype(str)
    paired_rescue_reason_mask = (
        filter_text.str.contains("branch_b_direction_discordant_with_candidate_state", regex=False)
        | filter_text.str.contains("signal_support_below_minimum", regex=False)
        | filter_text.str.contains("cnvseq_subchrom_boundary_weak_support", regex=False)
    )
    candidate_mask = (
        frame["keep_event"].astype(int).eq(0)
        & paired_rescue_reason_mask
        & frame["caller"].fillna("").astype(str).eq("wisecondorx_a_branch")
        & a_abs.ge(float(getattr(args, "paired_event_rescue_min_a_abs_z", 8.0)))
        & n_bins.ge(float(getattr(args, "paired_event_rescue_min_bins", 5)))
    )
    if not candidate_mask.any():
        return frame

    kept_mask = frame["keep_event"].astype(int).eq(1)
    mate_min_a = float(getattr(args, "paired_event_rescue_mate_min_a_abs_z", 50.0))
    for idx, row in frame.loc[candidate_mask].iterrows():
        mate_mask = (
            kept_mask
            & frame["sample_id"].astype(str).eq(str(row["sample_id"]))
            & frame["chrom"].astype(str).eq(str(row["chrom"]))
            & frame["state"].astype(str).str.lower().ne(str(row["state"]).lower())
            & pd.to_numeric(frame["a_abs_zscore"], errors="coerce").fillna(0.0).ge(mate_min_a)
        )
        if not mate_mask.any():
            continue
        frame.at[idx, "artifact_status"] = "review"
        frame.at[idx, "keep_event"] = int(bool(getattr(args, "keep_review", 1)))
        frame.at[idx, "technical_confidence"] = "low"
        frame.at[idx, "report_class"] = "candidate_review" if int(frame.at[idx, "keep_event"]) else "candidate_suppressed"
        frame.at[idx, "retain_reason"] = append_reason(frame.at[idx, "retain_reason"], "paired_chromosome_event_rescue")
        frame.at[idx, "downgrade_reason"] = append_reason(
            frame.at[idx, "downgrade_reason"], "paired_chromosome_event_rescue"
        )
        frame.at[idx, "filter_reason"] = remove_reason(
            frame.at[idx, "filter_reason"], "branch_b_direction_discordant_with_candidate_state"
        )
        frame.at[idx, "filter_reason"] = remove_reason(
            frame.at[idx, "filter_reason"], "signal_support_below_minimum"
        )
        frame.at[idx, "filter_reason"] = remove_reason(
            frame.at[idx, "filter_reason"], "cnvseq_subchrom_boundary_weak_support"
        )
        frame.at[idx, "manual_review_recommended"] = 1
    return frame


def apply_sca_group_context(events_df):
    if events_df.empty:
        return events_df
    frame = events_df.copy()
    required = {"sample_id", "chrom", "state", "keep_event", "artifact_flags", "biological_context"}
    if not required.issubset(frame.columns):
        return frame

    kept_chr_x = frame["keep_event"].astype(int).eq(1) & frame["chrom"].astype(str).eq("chrX")
    for (_, state), group in frame.loc[kept_chr_x].groupby(["sample_id", "state"], dropna=False):
        flags = group["artifact_flags"].fillna("").astype(str)
        has_p = flags.str.contains("sca_x_nonpar_p", regex=False)
        has_q = flags.str.contains("sca_x_nonpar_q", regex=False)
        if not has_p.any() or not has_q.any():
            continue
        idxs = group.loc[has_p | has_q].index
        for idx in idxs:
            frame.at[idx, "artifact_flags"] = append_flag(frame.at[idx, "artifact_flags"], "sca_x_nonpar_pq_pair")
            frame.at[idx, "downgrade_reason"] = append_reason(
                frame.at[idx, "downgrade_reason"], "sca_x_nonpar_pq_pair_review"
            )
            frame.at[idx, "manual_review_recommended"] = 1
            frame.at[idx, "biological_context"] = "sca_x_nonpar_pq_pair"
    return frame


def apply_sca_sex_consistency(events_df, gender_info, args):
    if events_df.empty:
        return events_df
    frame = events_df.copy()
    required = {"chrom", "state", "keep_event", "artifact_flags", "downgrade_reason", "filter_reason"}
    if not required.issubset(frame.columns):
        return frame

    sex_call = str(gender_info.get("sex_call", "") or "").strip().upper()
    if sex_call != "XY":
        return frame

    bam_x_relative = safe_float(gender_info.get("bam_x_relative_depth", np.nan))
    if not np.isfinite(bam_x_relative):
        return frame
    if bam_x_relative > float(getattr(args, "sca_xy_xgain_max_bam_x_relative", 0.80)):
        return frame

    flags = frame["artifact_flags"].fillna("").astype(str)
    chrom = frame["chrom"].fillna("").astype(str)
    state = frame["state"].fillna("").astype(str).str.lower()
    n_bins = pd.to_numeric(frame.get("n_bins", pd.Series(np.nan, index=frame.index)), errors="coerce")
    broad_or_pair = flags.str.contains("broad_chrom_fraction", regex=False) | flags.str.contains(
        "sca_x_nonpar_pq_pair", regex=False
    )
    focal_edge = (
        flags.str.contains("edge_event", regex=False)
        & (
            flags.str.contains("par_overlap", regex=False)
            | flags.str.contains("mixed_par_nonpar", regex=False)
            | flags.str.contains("sca_par_boundary_overlap", regex=False)
        )
        & n_bins.le(int(getattr(args, "sca_xy_xgain_focal_edge_max_bins", 20)))
    )
    inconsistent_mask = (
        frame["keep_event"].astype(int).eq(1)
        & chrom.eq("chrX")
        & state.eq("gain")
        & flags.str.contains("sca_event", regex=False)
        & (broad_or_pair | focal_edge)
    )
    for idx in frame.loc[inconsistent_mask].index:
        frame.at[idx, "artifact_status"] = "artifact"
        frame.at[idx, "keep_event"] = 0
        frame.at[idx, "technical_confidence"] = "rejected"
        frame.at[idx, "report_class"] = "technical_artifact"
        frame.at[idx, "artifact_flags"] = append_flag(
            frame.at[idx, "artifact_flags"], "sca_xy_xgain_without_sample_x_elevation"
        )
        frame.at[idx, "downgrade_reason"] = append_reason(
            frame.at[idx, "downgrade_reason"], "sca_xy_xgain_without_sample_x_elevation"
        )
        frame.at[idx, "filter_reason"] = append_reason(
            frame.at[idx, "filter_reason"], "sca_xy_xgain_without_sample_x_elevation"
        )
        frame.at[idx, "manual_review_recommended"] = 1
    return frame


def main():
    args = parse_args()
    logger = setup_logger("cnv_artifact_rules", args.log or None)
    bins_df, events_df = read_bins_and_candidates(
        args.input_bins,
        args.input_candidates,
        bins_required_columns=["chrom", "bin_index"],
        empty_candidates_ok=True,
    )
    gender_info = parse_gender_info(args.gender_tsv)
    sex_call = str(gender_info.get("sex_call", "") or "").strip().upper()
    par_regions = parse_par_regions(args.par_region)

    if events_df.empty:
        empty_summary = pd.DataFrame(
            [{"sample_id": args.sample_id, "sex_call": sex_call or "", "artifact_status": "none", "event_count": 0}]
        )
        write_table(args.output_events, events_df)
        write_table(args.output_summary, empty_summary)
        write_json(
            args.output_json,
            {
                "sample_id": args.sample_id,
                "sex_call": sex_call or None,
                "genome_build": args.genome_build,
                "par_regions": args.par_region,
                "events": [],
                "technical_summary": {"event_count": 0, "kept_event_count": 0},
            },
        )
        logger.info("no calibrated candidate events")
        return

    chrom_sizes = bins_df.groupby("chrom")["bin_index"].max().add(1).to_dict()
    cnvseq_context = [
        summarize_cnvseq_event_context(row, bins_df)
        for row in events_df.itertuples(index=False)
    ]
    if cnvseq_context:
        context_df = pd.DataFrame(cnvseq_context)
        for column in context_df.columns:
            events_df[column] = context_df[column]
    decisions = []
    for row in events_df.itertuples(index=False):
        decisions.append(
            classify_event(
                row=row,
                chrom_bin_count=int(chrom_sizes.get(row.chrom, 0)),
                args=args,
                sex_call=sex_call,
                par_regions=par_regions,
            )
        )

    decision_df = pd.DataFrame(decisions)
    for column in decision_df.columns:
        events_df[column] = decision_df[column]
    events_df = apply_paired_chromosome_event_rescue(events_df, args)
    events_df = apply_sca_group_context(events_df)
    events_df = apply_sca_sex_consistency(events_df, gender_info, args)
    events_df = events_df.sort_values(
        by=["keep_event", "artifact_status", "priority_score", "n_bins"],
        ascending=[False, True, False, False],
    ).reset_index(drop=True)

    summary_rows = []
    for artifact_status, frame in events_df.groupby("artifact_status", dropna=False):
        summary_rows.append(
            {
                "sample_id": args.sample_id,
                "sex_call": sex_call or "",
                "artifact_status": artifact_status,
                "event_count": int(len(frame)),
                "kept_event_count": int(frame["keep_event"].sum()),
                "top_priority_score": float(frame["priority_score"].max()),
            }
        )
    summary_df = pd.DataFrame(summary_rows)
    kept = events_df[events_df["keep_event"] == 1].copy()
    kept_preview = kept.head(5)[
        ["event_id", "chrom", "start", "end", "state", "artifact_status", "technical_confidence", "priority_score"]
    ].to_dict(orient="records")

    write_table(args.output_events, events_df)
    write_table(args.output_summary, summary_df)
    write_json(
        args.output_json,
        {
            "sample_id": args.sample_id,
            "sex_call": sex_call or None,
            "genome_build": args.genome_build,
            "par_regions": args.par_region,
            "technical_summary": {
                "event_count": int(len(events_df)),
                "kept_event_count": int(events_df["keep_event"].sum()),
                "pass_event_count": int((events_df["artifact_status"] == "pass").sum()),
                "review_event_count": int((events_df["artifact_status"] == "review").sum()),
                "artifact_event_count": int((events_df["artifact_status"] == "artifact").sum()),
            },
            "kept_event_preview": kept_preview,
            "events": events_df.to_dict(orient="records"),
        },
    )
    logger.info(
        "artifact review finished: sample=%s sex_call=%s events=%d kept=%d",
        args.sample_id,
        sex_call or "NA",
        len(events_df),
        int(events_df["keep_event"].sum()),
    )


if __name__ == "__main__":
    main()
