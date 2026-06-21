#!/biosoftware/miniconda/envs/snakemake_env/bin/python
from __future__ import annotations

import argparse
import math

import pandas as pd

from pgta.core.logging import setup_logger
from pgta.predict.branch_b.common import read_table, write_json, write_table


V2_CLASSIFIER_COLUMNS = [
    "v2_classifier_version",
    "v2_candidate_class",
    "v2_classifier_action",
    "v2_classifier_reason",
    "v2_evidence_tier",
    "v2_evidence_gate",
    "v2_review_priority",
    "v2_signal_strength_tier",
    "v2_length_tier",
    "v2_clean_support_label",
    "v2_gc_rc_context_label",
    "v2_background_context_label",
    "v2_background_context_reason",
    "v2_direction_support_label",
    "v2_direction_support_reason",
    "v2_b_signal_context_label",
    "v2_b_signal_context_reason",
    "v2_disposition",
    "v2_filter_version",
    "v2_filter_action",
    "v2_filter_reason",
    "v2_filter_scope",
    "v2_filter_hard_suppression_allowed",
    "v2_burden_reduction_version",
    "v2_burden_reduction_tier",
    "v2_burden_reduction_action",
    "v2_burden_reduction_reason",
    "v2_burden_evidence_tags",
    "v2_report_layer_version",
    "v2_report_layer_class",
    "v2_report_visibility",
    "v2_report_filter_reason",
    "v2_report_filter_rule_tags",
    "v2_report_filter_evidence_count",
    "v2_final_report_impact",
]

V2_FILTER_VERSION = "branch_b_v2_truth_safe_filter_v1"
V2_BURDEN_REDUCTION_VERSION = "branch_b_v2_burden_stratification_v1"
V2_REPORT_LAYER_VERSION = "branch_b_v2_report_layer_filter_v2"
PASS2_MULTI_REPORT_EVENT_THRESHOLD = 3
PASS2_VERY_STRONG_A_Z = 50.0


def parse_args():
    parser = argparse.ArgumentParser(description="Classify Branch A candidates with Branch B V2 shadow evidence.")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--input-ledger", required=True)
    parser.add_argument("--output-classification", required=True)
    parser.add_argument("--output-summary", required=True)
    parser.add_argument("--version", default="branch_ab_v2")
    parser.add_argument("--log", default="")
    return parser.parse_args()


def clean_text(value, default=""):
    try:
        if pd.isna(value):
            return default
    except (TypeError, ValueError):
        pass
    text = str(value if value is not None else "").strip()
    if not text or text.lower() == "nan":
        return default
    return text


def safe_float(value, default=math.nan):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return float(default)
    return number if math.isfinite(number) else float(default)


def normalized_disposition(row):
    return clean_text(row.get("final_disposition", "REVIEW_REQUIRED"), default="REVIEW_REQUIRED").upper()


def matched_background_status(row):
    return clean_text(row.get("matched_negative_background_status", ""), default="").upper()


def calibration_null_status(row):
    return clean_text(row.get("calibration_null_status", ""), default="").upper()


def has_unknown_matched_negative_background(row):
    status = matched_background_status(row)
    if status:
        return status == "UNKNOWN_BACKGROUND"
    source = clean_text(row.get("matched_negative_source", ""), default="").upper()
    return source == "UNKNOWN_BACKGROUND"


def state_direction(state):
    text = clean_text(state).lower()
    if text == "gain":
        return 1.0
    if text == "loss":
        return -1.0
    return 0.0


def is_sex_chromosome_candidate(row):
    chrom = clean_text(row.get("chrom", "")).lower().removeprefix("chr")
    return chrom in {"x", "y"}


def is_acrocentric_chromosome_candidate(row):
    chrom = clean_text(row.get("chrom", "")).lower().removeprefix("chr")
    return chrom in {"13", "14", "15", "21", "22"}


def signed_signal_support(row):
    direction = state_direction(row.get("state", ""))
    if direction == 0.0:
        return False
    for column in ["corrected_amplitude", "raw_amplitude", "a_zscore"]:
        value = safe_float(row.get(column, math.nan), default=math.nan)
        if math.isfinite(value) and abs(value) >= 2.0:
            return value * direction > 0.0
    return False


def branch_b_direction_support_label(row):
    if has_direction_conflict(row):
        return "B_SIGNAL_DISCORDANT_WITH_A_DIRECTION", "b_side_amplitude_opposite_a_direction_abs_ge_2"

    same_direction = safe_float(row.get("same_direction_fraction", math.nan), default=math.nan)
    if math.isfinite(same_direction) and same_direction >= 0.50:
        return "B_SIGNAL_SUPPORTED_A_DIRECTION", "same_direction_fraction_ge_0.50"

    direction = state_direction(row.get("state", ""))
    if direction != 0.0:
        for column in ["corrected_amplitude", "raw_amplitude"]:
            value = safe_float(row.get(column, math.nan), default=math.nan)
            if math.isfinite(value) and abs(value) >= 2.0 and value * direction > 0.0:
                return "B_SIGNAL_SUPPORTED_A_DIRECTION", f"{column}_same_direction_abs_ge_2"

    if has_positive_branch_a_support(row):
        return "A_ANCHORED_WEAK_B_SIGNAL", "positive_support_without_branch_b_signal_support"

    return "NO_POSITIVE_A_SIGNAL", "no_positive_a_signal_to_compare_b_signal"


def signal_strength_tier(row):
    a_abs_zscore = safe_float(row.get("a_abs_zscore", math.nan), default=math.nan)
    if not math.isfinite(a_abs_zscore):
        return "A_Z_MISSING"
    if a_abs_zscore >= 10.0:
        return "A_STRONG_Z_GE_10"
    if a_abs_zscore >= 5.0:
        return "A_SENSITIVE_Z_5_TO_10"
    return "A_WEAK_Z_LT_5"


def event_length_bp(row):
    start = safe_float(row.get("start", math.nan), default=math.nan)
    end = safe_float(row.get("end", math.nan), default=math.nan)
    if not math.isfinite(start) or not math.isfinite(end) or end <= start:
        return math.nan
    return float(end - start)


def length_tier(row):
    length = event_length_bp(row)
    if not math.isfinite(length):
        return "unknown_length"
    if length >= 10_000_000:
        return "large_ge10mb"
    if length >= 4_000_000:
        return "broad_review_ge4mb"
    if length >= 2_000_000:
        return "reportable_candidate_ge2mb"
    if length >= 1_000_000:
        return "review_only_ge1mb"
    return "focal_high_risk_lt1mb"


def clean_support_label(row):
    if has_ref_contract_risk(row):
        return "REF_CONTRACT_RISK"
    clean_fraction = safe_float(row.get("clean_bin_fraction", math.nan), default=math.nan)
    high_risk_fraction = safe_float(row.get("high_risk_bin_fraction", math.nan), default=math.nan)
    hard_region_fraction = safe_float(row.get("hard_region_fraction", math.nan), default=math.nan)
    if math.isfinite(clean_fraction):
        if clean_fraction < 0.20 and math.isfinite(high_risk_fraction) and high_risk_fraction >= 0.50:
            return "LOW_CLEAN_SUPPORT_HIGH_RISK"
        if clean_fraction >= 0.50 and (not math.isfinite(high_risk_fraction) or high_risk_fraction <= 0.25):
            return "CLEAN_SUPPORT_AVAILABLE"
        return "CLEAN_SUPPORT_REVIEW"
    if math.isfinite(hard_region_fraction):
        if hard_region_fraction >= 0.50:
            return "LOW_CLEAN_SUPPORT_HIGH_RISK"
        if hard_region_fraction <= 0.25:
            return "CLEAN_SUPPORT_AVAILABLE"
        return "CLEAN_SUPPORT_REVIEW"
    return "CLEAN_SUPPORT_UNKNOWN"


def gc_rc_context_label(row):
    attenuation = safe_float(row.get("attenuation_ratio", math.nan), default=math.nan)
    if not math.isfinite(attenuation):
        return "GC_RC_CONTEXT_UNKNOWN"
    if attenuation < 0.50:
        return "GC_RC_ATTENUATED_SEVERE"
    if attenuation < 0.80:
        return "GC_RC_ATTENUATED"
    if attenuation <= 1.25:
        return "GC_RC_STABLE"
    return "GC_RC_AMPLIFIED"


def background_context_label(row):
    background_status = matched_background_status(row)
    null_status = calibration_null_status(row)

    if background_status == "OK":
        percentile = matched_negative_percentile(row)
        if math.isfinite(percentile):
            return "MATCHED_NEGATIVE_BACKGROUND_INFORMATIVE", "matched_negative_percentile_available"
        return "MATCHED_NEGATIVE_BACKGROUND_WITHOUT_PERCENTILE", "matched_negative_without_percentile"

    if background_status == "SHADOW_BACKGROUND":
        if null_status == "NO_NULL_SUPPORT":
            return "SHADOW_BACKGROUND_NO_NULL_SUPPORT", "shadow_background_context_only_no_calibration_null"
        if null_status == "LIMITED_NULL_SUPPORT":
            return "SHADOW_BACKGROUND_LIMITED_NULL_SUPPORT", "shadow_background_context_only_limited_calibration_null"
        return "SHADOW_BACKGROUND_CONTEXT_ONLY", "shadow_background_context_only"

    if has_unknown_matched_negative_background(row):
        if null_status == "NO_NULL_SUPPORT":
            return "UNKNOWN_BACKGROUND_NO_NULL_SUPPORT", "no_matched_negative_and_no_calibration_null"
        if null_status == "LIMITED_NULL_SUPPORT":
            return "UNKNOWN_BACKGROUND_LIMITED_NULL_SUPPORT", "no_matched_negative_limited_calibration_null"
        if null_status == "OK":
            return "UNKNOWN_BACKGROUND_CALIBRATION_NULL_AVAILABLE", "no_matched_negative_but_calibration_null_available"
        return "UNKNOWN_BACKGROUND_CONTEXT_MISSING", "no_matched_negative_context"

    if null_status == "NO_NULL_SUPPORT":
        return "NO_MATCHED_BACKGROUND_NO_NULL_SUPPORT", "no_background_status_and_no_calibration_null"
    if null_status == "LIMITED_NULL_SUPPORT":
        return "NO_MATCHED_BACKGROUND_LIMITED_NULL_SUPPORT", "no_background_status_limited_calibration_null"
    if null_status == "OK":
        return "CALIBRATION_NULL_ONLY", "calibration_null_available_without_matched_negative"
    return "NO_BACKGROUND_CONTEXT", "no_background_or_calibration_context"


def has_direction_conflict(row):
    direction = state_direction(row.get("state", ""))
    if direction == 0.0:
        return False
    for column in ["corrected_amplitude", "raw_amplitude"]:
        value = safe_float(row.get(column, math.nan), default=math.nan)
        if math.isfinite(value) and abs(value) >= 2.0:
            return value * direction < 0.0
    return False


def has_positive_branch_a_support(row):
    a_abs_zscore = safe_float(row.get("a_abs_zscore", math.nan), default=math.nan)
    support_level = clean_text(row.get("a_support_level", "")).lower()
    same_direction = safe_float(row.get("same_direction_fraction", math.nan), default=math.nan)
    strong_a = math.isfinite(a_abs_zscore) and a_abs_zscore >= 10.0
    sensitive_a = math.isfinite(a_abs_zscore) and a_abs_zscore >= 5.0
    labeled_support = "strong" in support_level or "sensitive" in support_level
    direction_supported = (
        (math.isfinite(same_direction) and same_direction >= 0.50)
        or signed_signal_support(row)
    )
    return strong_a or labeled_support or (sensitive_a and direction_supported)


def has_ref_contract_risk(row):
    return clean_text(row.get("refmap_status", "")).upper() == "SAME_CHROM_REF_LEAKAGE"


def has_technical_review_risk(row):
    if has_ref_contract_risk(row):
        return True
    refmap_status = clean_text(row.get("refmap_status", "")).upper()
    if refmap_status in {"LOW_REF_BINS_HARD_RISK", "LOW_REF_BINS_REVIEW", "LOW_REFBIN_BURDEN"}:
        return True
    calibration_status = calibration_null_status(row)
    if calibration_status in {"NO_NULL_SUPPORT", "LIMITED_NULL_SUPPORT"}:
        return True
    hard_region_fraction = safe_float(row.get("hard_region_fraction", math.nan), default=math.nan)
    if math.isfinite(hard_region_fraction) and hard_region_fraction >= 0.50:
        return True
    return clean_text(row.get("sample_noise_status", "")).upper() == "HIGH_WAVINESS"


def matched_negative_percentile(row):
    for column in ["matched_negative_abs_percentile", "matched_negative_percentile"]:
        value = safe_float(row.get(column, math.nan), default=math.nan)
        if math.isfinite(value):
            return value
    return math.nan


def classify_evidence_tier(row):
    background_status = matched_background_status(row)
    positive_support = has_positive_branch_a_support(row)
    technical_risk = has_technical_review_risk(row)

    if has_unknown_matched_negative_background(row):
        if has_ref_contract_risk(row):
            return "UNKNOWN_BACKGROUND_REF_CONTRACT_RISK", "NO_CALL_CONTRACT_RISK", "MEDIUM"
        if positive_support:
            return "UNKNOWN_BACKGROUND_POSITIVE_SUPPORT", "NO_HARD_SUPPRESSION", "HIGH"
        if technical_risk:
            return "UNKNOWN_BACKGROUND_TECHNICAL_REVIEW", "NO_HARD_SUPPRESSION", "MEDIUM"
        return "UNKNOWN_BACKGROUND_REVIEW", "NO_HARD_SUPPRESSION", "LOW"

    if background_status == "OK":
        percentile = matched_negative_percentile(row)
        if math.isfinite(percentile) and percentile >= 0.99 and positive_support:
            return "MATCHED_NEGATIVE_OUTLIER_POSITIVE_SUPPORT", "BACKGROUND_INFORMATIVE", "HIGH"
        if math.isfinite(percentile) and percentile <= 0.95:
            return "MATCHED_NEGATIVE_BACKGROUND_COMPATIBLE", "BACKGROUND_INFORMATIVE", "LOW"
        if technical_risk:
            return "MATCHED_NEGATIVE_TECHNICAL_REVIEW", "BACKGROUND_INFORMATIVE", "MEDIUM"
        return "MATCHED_NEGATIVE_BORDERLINE_REVIEW", "BACKGROUND_INFORMATIVE", "MEDIUM"

    if background_status == "SHADOW_BACKGROUND":
        percentile = matched_negative_percentile(row)
        if math.isfinite(percentile) and percentile >= 0.99 and positive_support:
            return "SHADOW_BACKGROUND_OUTLIER_POSITIVE_SUPPORT", "SHADOW_BACKGROUND_CONTEXT", "HIGH"
        if math.isfinite(percentile) and percentile <= 0.95:
            return "SHADOW_BACKGROUND_COMPATIBLE", "SHADOW_BACKGROUND_CONTEXT", "LOW"
        if technical_risk:
            return "SHADOW_BACKGROUND_TECHNICAL_REVIEW", "SHADOW_BACKGROUND_CONTEXT", "MEDIUM"
        return "SHADOW_BACKGROUND_BORDERLINE_REVIEW", "SHADOW_BACKGROUND_CONTEXT", "MEDIUM"

    if positive_support:
        return "NO_MATCHED_BACKGROUND_POSITIVE_SUPPORT", "NO_BACKGROUND_INPUT", "HIGH"
    if technical_risk:
        return "NO_MATCHED_BACKGROUND_TECHNICAL_REVIEW", "NO_BACKGROUND_INPUT", "MEDIUM"
    return "NO_MATCHED_BACKGROUND_REVIEW", "NO_BACKGROUND_INPUT", "LOW"


def candidate_disposition(candidate_class, tier, background_label, length_label, clean_label):
    if candidate_class == "V2_SEX_CHROMOSOME_REVIEW":
        return "sca_branch_s_review"
    if candidate_class in {"V2_NO_CALL_CONTRACT_RISK", "V2_TECHNICAL_REVIEW"} or "REF_CONTRACT_RISK" in tier:
        return "technical_risk_review"
    if "POSITIVE_SUPPORT" in tier:
        if (
            background_label.startswith("UNKNOWN_BACKGROUND")
            or background_label.startswith("NO_")
            or "NO_NULL_SUPPORT" in background_label
        ):
            return "background_unknown_review"
        if length_label in {"large_ge10mb", "broad_review_ge4mb", "reportable_candidate_ge2mb"} and clean_label != "LOW_CLEAN_SUPPORT_HIGH_RISK":
            return "report_candidate"
        return "review_candidate"
    if "TECHNICAL_REVIEW" in tier:
        return "technical_risk_review"
    if (
        background_label.startswith("UNKNOWN_BACKGROUND")
        or background_label.startswith("NO_")
        or "NO_NULL_SUPPORT" in background_label
    ):
        return "background_unknown_review"
    return "review_candidate"


def candidate_filter_payload(candidate_class, tier, background_label, clean_label, disposition):
    if candidate_class == "V2_NO_CALL_CONTRACT_RISK" or "REF_CONTRACT_RISK" in tier:
        return {
            "v2_filter_version": V2_FILTER_VERSION,
            "v2_filter_action": "suppress_workflow_contract_risk",
            "v2_filter_reason": "same_chrom_ref_leakage_contract_risk",
            "v2_filter_scope": "workflow_contract_only",
            "v2_filter_hard_suppression_allowed": 1,
        }
    if candidate_class == "V2_SEX_CHROMOSOME_REVIEW":
        return {
            "v2_filter_version": V2_FILTER_VERSION,
            "v2_filter_action": "route_to_branch_s_review",
            "v2_filter_reason": "sex_chromosome_branch_s_review",
            "v2_filter_scope": "branch_s_review",
            "v2_filter_hard_suppression_allowed": 0,
        }
    if clean_label == "LOW_CLEAN_SUPPORT_HIGH_RISK" or disposition == "technical_risk_review":
        return {
            "v2_filter_version": V2_FILTER_VERSION,
            "v2_filter_action": "downgrade_to_technical_risk_review",
            "v2_filter_reason": "low_clean_support_or_technical_risk",
            "v2_filter_scope": "truth_safe_review_filter",
            "v2_filter_hard_suppression_allowed": 0,
        }
    if disposition == "report_candidate":
        return {
            "v2_filter_version": V2_FILTER_VERSION,
            "v2_filter_action": "keep_report_candidate",
            "v2_filter_reason": "positive_support_report_tier",
            "v2_filter_scope": "truth_safe_report_filter",
            "v2_filter_hard_suppression_allowed": 0,
        }
    if disposition == "background_unknown_review" or background_label.startswith(("UNKNOWN_BACKGROUND", "NO_")):
        return {
            "v2_filter_version": V2_FILTER_VERSION,
            "v2_filter_action": "keep_background_unknown_review",
            "v2_filter_reason": "background_unknown_truth_safe_review",
            "v2_filter_scope": "truth_safe_review_filter",
            "v2_filter_hard_suppression_allowed": 0,
        }
    if "BACKGROUND_COMPATIBLE" in tier:
        return {
            "v2_filter_version": V2_FILTER_VERSION,
            "v2_filter_action": "keep_background_compatible_review",
            "v2_filter_reason": "background_context_compatible_review_only",
            "v2_filter_scope": "truth_safe_review_filter",
            "v2_filter_hard_suppression_allowed": 0,
        }
    return {
        "v2_filter_version": V2_FILTER_VERSION,
        "v2_filter_action": "keep_review_candidate",
        "v2_filter_reason": "truth_safe_review_default",
        "v2_filter_scope": "truth_safe_review_filter",
        "v2_filter_hard_suppression_allowed": 0,
    }


def report_layer_rule_tags(row, length_label, clean_label, gc_label, b_signal_label):
    tags = []
    if length_label in {"focal_high_risk_lt1mb", "review_only_ge1mb"}:
        tags.append("short_or_focal")
    if clean_label == "LOW_CLEAN_SUPPORT_HIGH_RISK":
        tags.append("low_clean_high_risk")
    if b_signal_label in {
        "A_ANCHORED_WEAK_B_SIGNAL",
        "B_SIGNAL_DISCORDANT_WITH_A_DIRECTION",
        "NO_POSITIVE_A_SIGNAL",
    }:
        tags.append("b_signal_not_supportive")
    if gc_label in {"GC_RC_ATTENUATED", "GC_RC_ATTENUATED_SEVERE"}:
        tags.append("gc_rc_attenuated")
    return tags


def is_strong_a_signal(row):
    a_abs_zscore = safe_float(row.get("a_abs_zscore", math.nan), default=math.nan)
    return math.isfinite(a_abs_zscore) and a_abs_zscore >= 10.0


def is_sensitive_a_signal(row):
    a_abs_zscore = safe_float(row.get("a_abs_zscore", math.nan), default=math.nan)
    return math.isfinite(a_abs_zscore) and 5.0 <= a_abs_zscore < 10.0


def is_autosomal_chromosome(row):
    chrom = clean_text(row.get("chrom", ""))
    if chrom.startswith("chr"):
        chrom = chrom[3:]
    return chrom.isdigit() and 1 <= int(chrom) <= 22


def report_layer_payload(row, candidate_class, disposition, length_label, clean_label, gc_label, b_signal_label):
    if candidate_class == "V2_SEX_CHROMOSOME_REVIEW":
        return {
            "v2_report_layer_version": V2_REPORT_LAYER_VERSION,
            "v2_report_layer_class": "branch_s_event",
            "v2_report_visibility": "branch_s_report_section",
            "v2_report_filter_reason": "sex_chromosome_event_routed_to_branch_s",
            "v2_report_filter_rule_tags": "branch_s_route",
            "v2_report_filter_evidence_count": 1,
        }

    if candidate_class == "V2_NO_CALL_CONTRACT_RISK":
        return {
            "v2_report_layer_version": V2_REPORT_LAYER_VERSION,
            "v2_report_layer_class": "filtered_event",
            "v2_report_visibility": "audit_only",
            "v2_report_filter_reason": "workflow_reference_contract_risk",
            "v2_report_filter_rule_tags": "workflow_reference_contract_risk",
            "v2_report_filter_evidence_count": 1,
        }

    tags = report_layer_rule_tags(row, length_label, clean_label, gc_label, b_signal_label)
    strong_a_signal = is_strong_a_signal(row)
    report_layer_ready = (
        is_autosomal_chromosome(row)
        and strong_a_signal
        and b_signal_label == "B_SIGNAL_SUPPORTED_A_DIRECTION"
        and length_label in {"large_ge10mb", "broad_review_ge4mb", "reportable_candidate_ge2mb"}
        and clean_label == "CLEAN_SUPPORT_AVAILABLE"
    )
    if report_layer_ready:
        return {
            "v2_report_layer_version": V2_REPORT_LAYER_VERSION,
            "v2_report_layer_class": "report_event",
            "v2_report_visibility": "final_report",
            "v2_report_filter_reason": "strong_a_supported_report_layer_event",
            "v2_report_filter_rule_tags": "strong_a_signal;b_signal_supported;reportable_length;clean_support",
            "v2_report_filter_evidence_count": 4,
        }

    combined_filter = (
        not strong_a_signal
        and {"b_signal_not_supportive", "gc_rc_attenuated"}.issubset(tags)
        and (
            is_sensitive_a_signal(row)
            or "short_or_focal" in tags
            or "low_clean_high_risk" in tags
        )
    )
    if combined_filter:
        return {
            "v2_report_layer_version": V2_REPORT_LAYER_VERSION,
            "v2_report_layer_class": "filtered_event",
            "v2_report_visibility": "audit_only",
            "v2_report_filter_reason": "combined_sensitive_or_short_b_signal_gc_rc",
            "v2_report_filter_rule_tags": ";".join(tags),
            "v2_report_filter_evidence_count": len(tags),
        }

    if disposition == "report_candidate":
        return {
            "v2_report_layer_version": V2_REPORT_LAYER_VERSION,
            "v2_report_layer_class": "report_event",
            "v2_report_visibility": "final_report",
            "v2_report_filter_reason": "positive_support_with_report_layer_requirements",
            "v2_report_filter_rule_tags": "report_layer_pass",
            "v2_report_filter_evidence_count": 0,
        }

    return {
        "v2_report_layer_version": V2_REPORT_LAYER_VERSION,
        "v2_report_layer_class": "internal_review_event",
        "v2_report_visibility": "internal_review",
        "v2_report_filter_reason": "retained_for_internal_review",
        "v2_report_filter_rule_tags": ";".join(tags) if tags else "review_retained",
        "v2_report_filter_evidence_count": len(tags),
    }


def apply_report_layer_burden_pass2(classified_df):
    if classified_df is None or classified_df.empty or "v2_report_layer_class" not in classified_df.columns:
        return classified_df

    frame = classified_df.copy()
    report_event = frame["v2_report_layer_class"].fillna("").astype(str).eq("report_event")
    if not report_event.any():
        return frame

    report_counts = report_event.groupby(frame["sample_id"].fillna("").astype(str)).transform("sum")
    a_abs_zscore = pd.to_numeric(
        frame.get("a_abs_zscore", pd.Series(math.nan, index=frame.index)),
        errors="coerce",
    )
    background_label = frame.get("v2_background_context_label", pd.Series("", index=frame.index)).fillna("").astype(str)
    gc_label = frame.get("v2_gc_rc_context_label", pd.Series("", index=frame.index)).fillna("").astype(str)

    unknown_or_no_null_background = (
        background_label.str.startswith("UNKNOWN_BACKGROUND")
        | background_label.str.startswith("NO_")
        | background_label.str.contains("NO_NULL_SUPPORT", regex=False)
    )
    gc_rc_unstable = gc_label.isin({"GC_RC_ATTENUATED", "GC_RC_ATTENUATED_SEVERE", "GC_RC_AMPLIFIED"})
    not_very_strong_a = a_abs_zscore.lt(PASS2_VERY_STRONG_A_Z) | a_abs_zscore.isna()
    demote_mask = (
        report_event
        & report_counts.ge(PASS2_MULTI_REPORT_EVENT_THRESHOLD)
        & unknown_or_no_null_background
        & gc_rc_unstable
        & not_very_strong_a
    )

    if not demote_mask.any():
        return frame

    reason = "multi_report_unknown_background_gc_rc_unstable_internal_review"
    tags = "multi_report_sample;unknown_background_no_null;gc_rc_unstable;not_very_strong_a_signal"
    frame.loc[demote_mask, "v2_report_layer_class"] = "internal_review_event"
    frame.loc[demote_mask, "v2_report_visibility"] = "internal_review"
    frame.loc[demote_mask, "v2_report_filter_reason"] = reason
    frame.loc[demote_mask, "v2_report_filter_rule_tags"] = tags
    frame.loc[demote_mask, "v2_report_filter_evidence_count"] = 4
    frame.loc[demote_mask, "v2_filter_action"] = "downgrade_report_to_internal_review_multi_event_uncertain_context"
    frame.loc[demote_mask, "v2_filter_reason"] = reason
    frame.loc[demote_mask, "v2_filter_scope"] = "report_layer_burden_pass2"
    frame.loc[demote_mask, "v2_filter_hard_suppression_allowed"] = 0
    frame.loc[demote_mask, "v2_burden_reduction_tier"] = "background_unknown_review"
    frame.loc[demote_mask, "v2_burden_reduction_action"] = (
        "downgrade_report_to_internal_review_multi_event_uncertain_context"
    )
    frame.loc[demote_mask, "v2_burden_reduction_reason"] = reason
    return frame


def burden_evidence_tags(row, length_label, clean_label, gc_label):
    tags = [f"[CNVpro-inspired] length_tier={length_label}"]
    if is_acrocentric_chromosome_candidate(row):
        tags.append("[CNVpro-confirmed] acrocentric_qter_context_review_only")
    if is_sex_chromosome_candidate(row):
        tags.append("[CNVseq-asset] sex_homology_PAR_annotation_branch_s_context")
    if clean_label != "CLEAN_SUPPORT_UNKNOWN":
        tags.append(f"[CNVseq-asset] mask_mappability_annotation_only={clean_label}")
    if gc_label != "GC_RC_CONTEXT_UNKNOWN":
        tags.append(f"[CNVpro-like] gc_rc_context={gc_label}")
    tags.append("[Not used] CNVcalling_R_cghFLasso_not_primary_caller")
    return ";".join(tags)


def burden_reduction_payload(candidate_class, background_label, clean_label, disposition):
    if candidate_class == "V2_NO_CALL_CONTRACT_RISK":
        return {
            "v2_burden_reduction_version": V2_BURDEN_REDUCTION_VERSION,
            "v2_burden_reduction_tier": "technical_risk_review",
            "v2_burden_reduction_action": "suppress_workflow_contract_risk",
            "v2_burden_reduction_reason": "workflow_reference_contract_risk_only_hard_suppressible",
        }
    if candidate_class == "V2_SEX_CHROMOSOME_REVIEW":
        return {
            "v2_burden_reduction_version": V2_BURDEN_REDUCTION_VERSION,
            "v2_burden_reduction_tier": "branch_s_review",
            "v2_burden_reduction_action": "route_to_branch_s_review",
            "v2_burden_reduction_reason": "sex_chromosome_candidate_requires_branch_s_review",
        }
    if clean_label == "LOW_CLEAN_SUPPORT_HIGH_RISK" or disposition == "technical_risk_review":
        return {
            "v2_burden_reduction_version": V2_BURDEN_REDUCTION_VERSION,
            "v2_burden_reduction_tier": "technical_risk_review",
            "v2_burden_reduction_action": "downgrade_to_technical_risk_review",
            "v2_burden_reduction_reason": "low_clean_support_or_technical_risk_review_only",
        }
    if disposition == "report_candidate":
        return {
            "v2_burden_reduction_version": V2_BURDEN_REDUCTION_VERSION,
            "v2_burden_reduction_tier": "report_candidate",
            "v2_burden_reduction_action": "keep_report_candidate",
            "v2_burden_reduction_reason": "positive_support_with_reportable_review_tier",
        }
    if disposition == "background_unknown_review" or background_label.startswith(("UNKNOWN_BACKGROUND", "NO_")):
        return {
            "v2_burden_reduction_version": V2_BURDEN_REDUCTION_VERSION,
            "v2_burden_reduction_tier": "background_unknown_review",
            "v2_burden_reduction_action": "stratify_background_unknown_review",
            "v2_burden_reduction_reason": "background_unknown_review_preserved_for_truth_safety",
        }
    return {
        "v2_burden_reduction_version": V2_BURDEN_REDUCTION_VERSION,
        "v2_burden_reduction_tier": "review_candidate",
        "v2_burden_reduction_action": "stratify_review_candidate",
        "v2_burden_reduction_reason": "truth_safe_review_candidate",
    }


def candidate_context_payload(row, tier, gate, priority, candidate_class, action, reason):
    signal_label, signal_reason = branch_b_direction_support_label(row)
    background_label, background_reason = background_context_label(row)
    strength_label = signal_strength_tier(row)
    length_label = length_tier(row)
    clean_label = clean_support_label(row)
    gc_label = gc_rc_context_label(row)
    disposition = candidate_disposition(candidate_class, tier, background_label, length_label, clean_label)
    filter_payload = candidate_filter_payload(candidate_class, tier, background_label, clean_label, disposition)
    burden_payload = burden_reduction_payload(candidate_class, background_label, clean_label, disposition)
    evidence_tags = burden_evidence_tags(row, length_label, clean_label, gc_label)
    report_payload = report_layer_payload(
        row,
        candidate_class,
        disposition,
        length_label,
        clean_label,
        gc_label,
        signal_label,
    )
    if report_payload["v2_report_layer_class"] == "filtered_event" and not report_payload[
        "v2_report_filter_rule_tags"
    ].startswith("workflow_reference_contract_risk"):
        filter_payload = {
            "v2_filter_version": V2_FILTER_VERSION,
            "v2_filter_action": "filter_report_layer_combined_technical_risk",
            "v2_filter_reason": report_payload["v2_report_filter_reason"],
            "v2_filter_scope": "report_layer_filter",
            "v2_filter_hard_suppression_allowed": 0,
        }
    elif report_payload["v2_report_layer_class"] == "report_event":
        filter_payload = {
            "v2_filter_version": V2_FILTER_VERSION,
            "v2_filter_action": "keep_report_layer_event",
            "v2_filter_reason": report_payload["v2_report_filter_reason"],
            "v2_filter_scope": "report_layer_filter",
            "v2_filter_hard_suppression_allowed": 0,
        }
    return {
        "v2_candidate_class": candidate_class,
        "v2_classifier_action": action,
        "v2_classifier_reason": reason,
        "v2_evidence_tier": tier,
        "v2_evidence_gate": gate,
        "v2_review_priority": priority,
        "v2_signal_strength_tier": strength_label,
        "v2_length_tier": length_label,
        "v2_clean_support_label": clean_label,
        "v2_gc_rc_context_label": gc_label,
        "v2_background_context_label": background_label,
        "v2_background_context_reason": background_reason,
        "v2_direction_support_label": signal_label,
        "v2_direction_support_reason": signal_reason,
        "v2_b_signal_context_label": signal_label,
        "v2_b_signal_context_reason": signal_reason,
        "v2_disposition": disposition,
        **filter_payload,
        **burden_payload,
        "v2_burden_evidence_tags": evidence_tags,
        **report_payload,
    }


def classify_candidate_row(row):
    tier, gate, priority = classify_evidence_tier(row)

    if is_sex_chromosome_candidate(row):
        return candidate_context_payload(
            row,
            tier,
            gate,
            priority,
            "V2_SEX_CHROMOSOME_REVIEW",
            "V2_ROUTE_BRANCH_S_REVIEW",
            f"sex_chromosome_branch_s_review:{tier.lower()}",
        )
    if tier.endswith("REF_CONTRACT_RISK"):
        return candidate_context_payload(
            row,
            tier,
            gate,
            priority,
            "V2_NO_CALL_CONTRACT_RISK",
            "V2_REVIEW_NO_HARD_SUPPRESSION",
            tier.lower(),
        )
    if "POSITIVE_SUPPORT" in tier:
        return candidate_context_payload(
            row,
            tier,
            gate,
            priority,
            "V2_POSITIVE_SUPPORT_REVIEW",
            "V2_REVIEW_POSITIVE_SUPPORT",
            tier.lower(),
        )
    if "BACKGROUND_COMPATIBLE" in tier or tier == "SHADOW_BACKGROUND_COMPATIBLE":
        return candidate_context_payload(
            row,
            tier,
            gate,
            priority,
            "V2_BACKGROUND_COMPATIBLE_REVIEW",
            "V2_REVIEW_BACKGROUND_COMPATIBLE",
            tier.lower(),
        )
    if "TECHNICAL_REVIEW" in tier:
        return candidate_context_payload(
            row,
            tier,
            gate,
            priority,
            "V2_TECHNICAL_REVIEW",
            "V2_REVIEW_TECHNICAL_RISK",
            tier.lower(),
        )
    return candidate_context_payload(
        row,
        tier,
        gate,
        priority,
        "V2_REVIEW_REQUIRED",
        "V2_REVIEW_ONLY",
        tier.lower(),
    )


def classify_branch_b_v2_candidates(candidate_ledger, version="branch_ab_v2"):
    if candidate_ledger is None or candidate_ledger.empty:
        columns = list(candidate_ledger.columns if candidate_ledger is not None else []) + V2_CLASSIFIER_COLUMNS
        return pd.DataFrame(columns=columns)
    if "candidate_id" not in candidate_ledger.columns:
        raise ValueError("Branch B V2 classifier input missing candidate_id column")

    frame = candidate_ledger.copy()
    payloads = []
    for _, row in frame.iterrows():
        payload = classify_candidate_row(row)
        payload["v2_classifier_version"] = str(version)
        payload["v2_final_report_impact"] = "none_shadow_only"
        payloads.append(payload)
    payload_df = pd.DataFrame(payloads)
    classified = pd.concat([frame.reset_index(drop=True), payload_df[V2_CLASSIFIER_COLUMNS]], axis=1)
    return apply_report_layer_burden_pass2(classified)


def summarize_v2_classification(sample_id, classified_df, version="branch_ab_v2"):
    class_counts = (
        classified_df["v2_candidate_class"].fillna("UNKNOWN").astype(str).value_counts().sort_index().to_dict()
        if "v2_candidate_class" in classified_df.columns
        else {}
    )
    action_counts = (
        classified_df["v2_classifier_action"].fillna("UNKNOWN").astype(str).value_counts().sort_index().to_dict()
        if "v2_classifier_action" in classified_df.columns
        else {}
    )
    evidence_tier_counts = (
        classified_df["v2_evidence_tier"].fillna("UNKNOWN").astype(str).value_counts().sort_index().to_dict()
        if "v2_evidence_tier" in classified_df.columns
        else {}
    )
    review_priority_counts = (
        classified_df["v2_review_priority"].fillna("UNKNOWN").astype(str).value_counts().sort_index().to_dict()
        if "v2_review_priority" in classified_df.columns
        else {}
    )
    direction_support_label_counts = (
        classified_df["v2_direction_support_label"].fillna("UNKNOWN").astype(str).value_counts().sort_index().to_dict()
        if "v2_direction_support_label" in classified_df.columns
        else {}
    )
    background_context_label_counts = (
        classified_df["v2_background_context_label"].fillna("UNKNOWN").astype(str).value_counts().sort_index().to_dict()
        if "v2_background_context_label" in classified_df.columns
        else {}
    )
    b_signal_context_label_counts = (
        classified_df["v2_b_signal_context_label"].fillna("UNKNOWN").astype(str).value_counts().sort_index().to_dict()
        if "v2_b_signal_context_label" in classified_df.columns
        else {}
    )
    disposition_counts = (
        classified_df["v2_disposition"].fillna("UNKNOWN").astype(str).value_counts().sort_index().to_dict()
        if "v2_disposition" in classified_df.columns
        else {}
    )
    filter_action_counts = (
        classified_df["v2_filter_action"].fillna("UNKNOWN").astype(str).value_counts().sort_index().to_dict()
        if "v2_filter_action" in classified_df.columns
        else {}
    )
    burden_reduction_tier_counts = (
        classified_df["v2_burden_reduction_tier"].fillna("UNKNOWN").astype(str).value_counts().sort_index().to_dict()
        if "v2_burden_reduction_tier" in classified_df.columns
        else {}
    )
    burden_reduction_action_counts = (
        classified_df["v2_burden_reduction_action"].fillna("UNKNOWN").astype(str).value_counts().sort_index().to_dict()
        if "v2_burden_reduction_action" in classified_df.columns
        else {}
    )
    length_tier_counts = (
        classified_df["v2_length_tier"].fillna("UNKNOWN").astype(str).value_counts().sort_index().to_dict()
        if "v2_length_tier" in classified_df.columns
        else {}
    )
    clean_support_label_counts = (
        classified_df["v2_clean_support_label"].fillna("UNKNOWN").astype(str).value_counts().sort_index().to_dict()
        if "v2_clean_support_label" in classified_df.columns
        else {}
    )
    gc_rc_context_label_counts = (
        classified_df["v2_gc_rc_context_label"].fillna("UNKNOWN").astype(str).value_counts().sort_index().to_dict()
        if "v2_gc_rc_context_label" in classified_df.columns
        else {}
    )
    burden_evidence_tag_counts = {}
    if "v2_burden_evidence_tags" in classified_df.columns:
        for tag_text in classified_df["v2_burden_evidence_tags"].fillna("").astype(str):
            for tag in [part.strip() for part in tag_text.split(";") if part.strip()]:
                burden_evidence_tag_counts[tag] = burden_evidence_tag_counts.get(tag, 0) + 1
        burden_evidence_tag_counts = dict(sorted(burden_evidence_tag_counts.items()))
    filter_hard_suppression_allowed_count = (
        int(pd.to_numeric(classified_df["v2_filter_hard_suppression_allowed"], errors="coerce").fillna(0).sum())
        if "v2_filter_hard_suppression_allowed" in classified_df.columns
        else 0
    )
    return {
        "sample_id": str(sample_id),
        "version": str(version),
        "candidate_count": int(len(classified_df)),
        "class_counts": {str(key): int(value) for key, value in class_counts.items()},
        "action_counts": {str(key): int(value) for key, value in action_counts.items()},
        "evidence_tier_counts": {str(key): int(value) for key, value in evidence_tier_counts.items()},
        "review_priority_counts": {str(key): int(value) for key, value in review_priority_counts.items()},
        "background_context_label_counts": {
            str(key): int(value) for key, value in background_context_label_counts.items()
        },
        "b_signal_context_label_counts": {
            str(key): int(value) for key, value in b_signal_context_label_counts.items()
        },
        "disposition_counts": {
            str(key): int(value) for key, value in disposition_counts.items()
        },
        "filter_action_counts": {
            str(key): int(value) for key, value in filter_action_counts.items()
        },
        "burden_reduction_tier_counts": {
            str(key): int(value) for key, value in burden_reduction_tier_counts.items()
        },
        "burden_reduction_action_counts": {
            str(key): int(value) for key, value in burden_reduction_action_counts.items()
        },
        "length_tier_counts": {
            str(key): int(value) for key, value in length_tier_counts.items()
        },
        "clean_support_label_counts": {
            str(key): int(value) for key, value in clean_support_label_counts.items()
        },
        "gc_rc_context_label_counts": {
            str(key): int(value) for key, value in gc_rc_context_label_counts.items()
        },
        "burden_evidence_tag_counts": {
            str(key): int(value) for key, value in burden_evidence_tag_counts.items()
        },
        "filter_hard_suppression_allowed_count": filter_hard_suppression_allowed_count,
        "direction_support_label_counts": {
            str(key): int(value) for key, value in direction_support_label_counts.items()
        },
        "final_report_impact": "none_shadow_only",
        "classified_only_branch_a_candidates": True,
        "legacy_decision_fields_ignored": True,
        "ignored_legacy_decision_fields": [
            "final_disposition",
            "branch_b_keep_event",
            "branch_b_report_class",
            "branch_b_artifact_status",
        ],
    }


def main():
    args = parse_args()
    logger = setup_logger("branch_b_v2_classifier", args.log or None)
    ledger = read_table(args.input_ledger, required_columns=["candidate_id"], empty_ok=True)
    classified = classify_branch_b_v2_candidates(ledger, version=args.version)
    write_table(args.output_classification, classified)
    write_json(args.output_summary, summarize_v2_classification(args.sample_id, classified, version=args.version))
    logger.info("wrote Branch B V2 shadow classifications rows=%d to %s", len(classified), args.output_classification)


if __name__ == "__main__":
    main()
