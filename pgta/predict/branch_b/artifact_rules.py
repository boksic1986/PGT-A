#!/biosoftware/miniconda/envs/snakemake_env/bin/python
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from pgta.predict.branch_b.common import read_bins_and_candidates, read_table, write_json, write_table
from pgta.core.logging import setup_logger


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
    parser.add_argument("--log", default="")
    return parser.parse_args()


def parse_gender_tsv(path_value):
    if not path_value:
        return ""
    path = Path(path_value)
    if not path.exists():
        return ""
    df = read_table(path, empty_ok=True)
    if df.empty or "sex_call" not in df.columns:
        return ""
    return str(df["sex_call"].iloc[0]).strip().upper()


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


def safe_float(value, default=np.nan):
    try:
        return float(value)
    except (TypeError, ValueError):
        return float(default)


def classify_event(row, chrom_bin_count, args, sex_call, par_regions):
    flags = []
    explanations = []
    retain_reasons = []
    downgrade_reasons = []
    filter_reasons = []
    hard_artifact = False
    review_only = False
    chrom = str(row.chrom)
    is_sex_chrom = chrom in {"chrX", "chrY"}
    chrom_fraction = float(row.n_bins) / max(chrom_bin_count, 1)
    overlap_bp, par_fraction = compute_par_overlap(chrom, int(row.start), int(row.end), par_regions)
    mean_z = abs(safe_float(getattr(row, "calibrated_mean_z", np.nan)))
    median_z = abs(safe_float(getattr(row, "calibrated_median_z", np.nan)))
    adjusted_event_z = abs(safe_float(getattr(row, "event_corr_adjusted_z", np.nan)))
    max_calibrated_z = max(mean_z, median_z, adjusted_event_z)
    internal_support_z = max(
        max_calibrated_z,
        abs(safe_float(getattr(row, "segment_mean_robust_z", np.nan), default=0.0)),
        abs(safe_float(getattr(row, "segment_median_robust_z", np.nan), default=0.0)),
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
    preserve_broad_internal_signal = (
        chrom_fraction > args.max_chrom_fraction
        and internal_support_z >= args.broad_support_min_abs_z
        and effective_bin_count >= args.broad_support_min_effective_bins
        and clean_fraction >= args.broad_support_min_clean_fraction
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
    weighted_non_high_support_fraction = clean_fraction + (0.5 * max(moderate_fraction, 0.0))
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

    if int(row.n_bins) < args.min_event_bins:
        flags.append("too_few_bins")
        explanations.append("Event spans too few bins for a stable segment.")
        filter_reasons.append("bin_count_below_minimum")
        hard_artifact = True
    if max_calibrated_z < args.min_abs_calibrated_z:
        flags.append("low_calibrated_signal")
        explanations.append("Calibrated signal amplitude is below the minimum support threshold.")
        filter_reasons.append("signal_support_below_minimum")
        hard_artifact = True
    if chrom_fraction > args.max_chrom_fraction:
        flags.append("broad_chrom_fraction")
        if is_sex_chrom and sex_call in {"XX", "XY"}:
            explanations.append("Large sex-chromosome event requires sex-aware review.")
            downgrade_reasons.append("broad_sex_chromosome_event")
            review_only = True
        elif preserve_broad_internal_signal:
            explanations.append("Broad event is preserved for review because Branch B shows high-confidence chromosome-scale support.")
            downgrade_reasons.append("broad_event_preserved_by_internal_support")
            review_only = True
        else:
            explanations.append("Event spans too much of one chromosome and is likely technical.")
            filter_reasons.append("chromosome_fraction_too_large")
            hard_artifact = True
    if int(row.start_bin) <= args.edge_bin_window or (chrom_bin_count - int(row.end_bin) - 1) <= args.edge_bin_window:
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
        explanations.append("Event overlaps a blacklisted genomic region.")
        filter_reasons.append("blacklist_region_overlap")
        hard_artifact = True
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
    if int(safe_float(getattr(row, "high_risk_boundary_crossing", 0), default=0.0)) == 1:
        flags.append("high_risk_boundary_crossing")
        explanations.append("Event crosses a clean/high-risk boundary and should be manually reviewed.")
        downgrade_reasons.append("high_risk_boundary_crossing")
        review_only = True

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

    biological_context = "sex_chromosome" if is_sex_chrom else "autosome"
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
        "manual_review_recommended": int(review_only or "manual" in " ".join(explanations).lower() or bool(downgrade_reasons)),
        "biological_context": biological_context,
    }


def main():
    args = parse_args()
    logger = setup_logger("cnv_artifact_rules", args.log or None)
    bins_df, events_df = read_bins_and_candidates(
        args.input_bins,
        args.input_candidates,
        bins_required_columns=["chrom", "bin_index"],
        empty_candidates_ok=True,
    )
    sex_call = parse_gender_tsv(args.gender_tsv)
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
