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
    "v2_direction_support_label",
    "v2_direction_support_reason",
    "v2_final_report_impact",
]


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
        return "B_DIRECTION_CONFLICT", "amplitude_opposite_direction_abs_ge_2"

    same_direction = safe_float(row.get("same_direction_fraction", math.nan), default=math.nan)
    if math.isfinite(same_direction) and same_direction >= 0.50:
        return "B_DIRECTION_SUPPORTED", "same_direction_fraction_ge_0.50"

    direction = state_direction(row.get("state", ""))
    if direction != 0.0:
        for column in ["corrected_amplitude", "raw_amplitude"]:
            value = safe_float(row.get(column, math.nan), default=math.nan)
            if math.isfinite(value) and abs(value) >= 2.0 and value * direction > 0.0:
                return "B_DIRECTION_SUPPORTED", f"{column}_same_direction_abs_ge_2"

    if has_positive_branch_a_support(row):
        return "A_ONLY_WEAK_B_DIRECTION", "positive_support_without_branch_b_direction_support"

    return "NO_POSITIVE_SUPPORT", "no_positive_support_to_compare_direction"


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
    calibration_status = clean_text(row.get("calibration_null_status", "")).upper()
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
        if has_direction_conflict(row):
            return "UNKNOWN_BACKGROUND_DIRECTION_CONFLICT", "NO_CALL_CONTRACT_RISK", "MEDIUM"
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


def classify_candidate_row(row):
    tier, gate, priority = classify_evidence_tier(row)
    direction_label, direction_reason = branch_b_direction_support_label(row)

    if is_sex_chromosome_candidate(row):
        return {
            "v2_candidate_class": "V2_SEX_CHROMOSOME_REVIEW",
            "v2_classifier_action": "V2_ROUTE_BRANCH_S_REVIEW",
            "v2_classifier_reason": f"sex_chromosome_branch_s_review:{tier.lower()}",
            "v2_evidence_tier": tier,
            "v2_evidence_gate": gate,
            "v2_review_priority": priority,
            "v2_direction_support_label": direction_label,
            "v2_direction_support_reason": direction_reason,
        }
    if tier.endswith("REF_CONTRACT_RISK") or tier.endswith("DIRECTION_CONFLICT"):
        return {
            "v2_candidate_class": "V2_NO_CALL_CONTRACT_RISK",
            "v2_classifier_action": "V2_REVIEW_NO_HARD_SUPPRESSION",
            "v2_classifier_reason": tier.lower(),
            "v2_evidence_tier": tier,
            "v2_evidence_gate": gate,
            "v2_review_priority": priority,
            "v2_direction_support_label": direction_label,
            "v2_direction_support_reason": direction_reason,
        }
    if "POSITIVE_SUPPORT" in tier:
        return {
            "v2_candidate_class": "V2_POSITIVE_SUPPORT_REVIEW",
            "v2_classifier_action": "V2_REVIEW_POSITIVE_SUPPORT",
            "v2_classifier_reason": tier.lower(),
            "v2_evidence_tier": tier,
            "v2_evidence_gate": gate,
            "v2_review_priority": priority,
            "v2_direction_support_label": direction_label,
            "v2_direction_support_reason": direction_reason,
        }
    if "BACKGROUND_COMPATIBLE" in tier or tier == "SHADOW_BACKGROUND_COMPATIBLE":
        return {
            "v2_candidate_class": "V2_BACKGROUND_COMPATIBLE_REVIEW",
            "v2_classifier_action": "V2_REVIEW_BACKGROUND_COMPATIBLE",
            "v2_classifier_reason": tier.lower(),
            "v2_evidence_tier": tier,
            "v2_evidence_gate": gate,
            "v2_review_priority": priority,
            "v2_direction_support_label": direction_label,
            "v2_direction_support_reason": direction_reason,
        }
    if "TECHNICAL_REVIEW" in tier:
        return {
            "v2_candidate_class": "V2_TECHNICAL_REVIEW",
            "v2_classifier_action": "V2_REVIEW_TECHNICAL_RISK",
            "v2_classifier_reason": tier.lower(),
            "v2_evidence_tier": tier,
            "v2_evidence_gate": gate,
            "v2_review_priority": priority,
            "v2_direction_support_label": direction_label,
            "v2_direction_support_reason": direction_reason,
        }
    return {
        "v2_candidate_class": "V2_REVIEW_REQUIRED",
        "v2_classifier_action": "V2_REVIEW_ONLY",
        "v2_classifier_reason": tier.lower(),
        "v2_evidence_tier": tier,
        "v2_evidence_gate": gate,
        "v2_review_priority": priority,
        "v2_direction_support_label": direction_label,
        "v2_direction_support_reason": direction_reason,
    }


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
    return pd.concat([frame.reset_index(drop=True), payload_df[V2_CLASSIFIER_COLUMNS]], axis=1)


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
    return {
        "sample_id": str(sample_id),
        "version": str(version),
        "candidate_count": int(len(classified_df)),
        "class_counts": {str(key): int(value) for key, value in class_counts.items()},
        "action_counts": {str(key): int(value) for key, value in action_counts.items()},
        "evidence_tier_counts": {str(key): int(value) for key, value in evidence_tier_counts.items()},
        "review_priority_counts": {str(key): int(value) for key, value in review_priority_counts.items()},
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
