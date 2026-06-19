#!/biosoftware/miniconda/envs/snakemake_env/bin/python
from __future__ import annotations

import argparse

import pandas as pd

from pgta.core.logging import setup_logger
from pgta.predict.branch_b.common import read_table, write_json, write_table


V2_CLASSIFIER_COLUMNS = [
    "v2_classifier_version",
    "v2_candidate_class",
    "v2_classifier_action",
    "v2_classifier_reason",
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
    text = str(value if value is not None else "").strip()
    if not text or text.lower() == "nan":
        return default
    return text


def normalized_disposition(row):
    return clean_text(row.get("final_disposition", "REVIEW_REQUIRED"), default="REVIEW_REQUIRED").upper()


def matched_background_status(row):
    return clean_text(row.get("matched_negative_background_status", ""), default="").upper()


def has_unknown_matched_negative_background(row):
    if matched_background_status(row) == "UNKNOWN_BACKGROUND":
        return True
    source = clean_text(row.get("matched_negative_source", ""), default="").upper()
    return source == "UNKNOWN_BACKGROUND"


def classify_candidate_row(row):
    disposition = normalized_disposition(row)

    if has_unknown_matched_negative_background(row):
        return {
            "v2_candidate_class": "REVIEW_REQUIRED",
            "v2_classifier_action": "SHADOW_REVIEW_ONLY",
            "v2_classifier_reason": "unknown_matched_negative_background",
        }
    if disposition == "CONFIRMED":
        return {
            "v2_candidate_class": "CONFIRMED_SHADOW",
            "v2_classifier_action": "SHADOW_CONFIRM",
            "v2_classifier_reason": "legacy_confirmed_with_shadow_evidence",
        }
    if disposition == "MOSAIC_SUSPECT":
        return {
            "v2_candidate_class": "MOSAIC_SUSPECT_SHADOW",
            "v2_classifier_action": "SHADOW_MOSAIC_REVIEW",
            "v2_classifier_reason": "legacy_mosaic_suspect",
        }
    if disposition == "LIKELY_ARTIFACT":
        return {
            "v2_candidate_class": "LIKELY_ARTIFACT_SHADOW",
            "v2_classifier_action": "SHADOW_ARTIFACT_REVIEW",
            "v2_classifier_reason": "legacy_artifact_with_shadow_evidence",
        }
    if disposition == "NO_CALL":
        return {
            "v2_candidate_class": "REVIEW_REQUIRED",
            "v2_classifier_action": "SHADOW_REVIEW_ONLY",
            "v2_classifier_reason": "legacy_no_call",
        }
    return {
        "v2_candidate_class": "REVIEW_REQUIRED",
        "v2_classifier_action": "SHADOW_REVIEW_ONLY",
        "v2_classifier_reason": "legacy_review_required_or_unclassified",
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
    return {
        "sample_id": str(sample_id),
        "version": str(version),
        "candidate_count": int(len(classified_df)),
        "class_counts": {str(key): int(value) for key, value in class_counts.items()},
        "action_counts": {str(key): int(value) for key, value in action_counts.items()},
        "final_report_impact": "none_shadow_only",
        "classified_only_branch_a_candidates": True,
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
