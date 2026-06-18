#!/biosoftware/miniconda/envs/snakemake_env/bin/python
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from pgta.core.logging import setup_logger
from pgta.predict.branch_b.common import read_table, write_json, write_table


VALID_LABELS = {"N0", "N1", "N2"}


def parse_args():
    parser = argparse.ArgumentParser(description="Assign Branch A/B V2 N0/N1/N2 negative-bank labels.")
    parser.add_argument("--input-samples", required=True)
    parser.add_argument("--output-labels", required=True)
    parser.add_argument("--output-summary", required=True)
    parser.add_argument("--version", default="branch_ab_v2")
    parser.add_argument("--log", default="")
    return parser.parse_args()


def clean_text(value, default=""):
    text = str(value if value is not None else "").strip()
    if not text or text.lower() == "nan":
        return default
    return text


def normalize_label(value):
    label = clean_text(value).upper()
    return label if label in VALID_LABELS else ""


def numeric_value(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return float(default)


def classify_negative_bank_row(row):
    manual_label = normalize_label(row.get("manual_negative_bank_label", ""))
    qc_status = clean_text(row.get("qc_status", "")).upper()
    sample_role = clean_text(row.get("sample_role", "")).lower()
    review_status = clean_text(row.get("manual_review_status", "")).lower()
    kept_count = numeric_value(row.get("branch_b_kept_count", 0), default=0.0)
    reasons = []

    if manual_label:
        reasons.append(f"manual_label_{manual_label}")
        if review_status:
            reasons.append(review_status)
        return manual_label, ";".join(reasons)

    if qc_status and qc_status != "PASS":
        return "N2", "qc_not_pass"
    if sample_role in {"known_positive", "positive", "abnormal", "truth_positive"}:
        return "N2", "not_negative_sample"
    if kept_count > 0:
        return "N2", "branch_b_retained_review_event"
    if "review_event" in review_status or "hold" in review_status or "recurring" in review_status:
        return "N2", review_status or "pending_manual_review"
    if sample_role in {"locked_negative", "clean_negative", "n0"} and review_status in {"", "clean", "locked_clean"}:
        return "N0", "locked_clean_negative"
    return "N1", "presumed_negative_requires_review"


def assign_negative_bank_labels(samples_df, version="branch_ab_v2"):
    if samples_df is None or samples_df.empty:
        return pd.DataFrame(
            columns=[
                "sample_id",
                "negative_bank_version",
                "negative_bank_label",
                "matched_negative_eligible",
                "negative_bank_reason",
            ]
        )
    if "sample_id" not in samples_df.columns:
        raise ValueError("negative-bank input missing sample_id column")
    frame = samples_df.copy()
    labels = []
    reasons = []
    for _, row in frame.iterrows():
        label, reason = classify_negative_bank_row(row)
        labels.append(label)
        reasons.append(reason)
    frame["negative_bank_version"] = str(version)
    frame["negative_bank_label"] = labels
    frame["matched_negative_eligible"] = frame["negative_bank_label"].eq("N0").astype(int)
    frame["negative_bank_reason"] = reasons
    ordered = [
        "sample_id",
        "negative_bank_version",
        "negative_bank_label",
        "matched_negative_eligible",
        "negative_bank_reason",
    ]
    return frame[ordered + [column for column in frame.columns if column not in ordered]]


def summarize_negative_bank(labeled_df, version="branch_ab_v2"):
    label_counts = (
        labeled_df["negative_bank_label"].fillna("UNKNOWN").astype(str).value_counts().sort_index().to_dict()
        if "negative_bank_label" in labeled_df.columns
        else {}
    )
    matched_count = (
        int(pd.to_numeric(labeled_df["matched_negative_eligible"], errors="coerce").fillna(0).sum())
        if "matched_negative_eligible" in labeled_df.columns
        else 0
    )
    return {
        "version": str(version),
        "sample_count": int(len(labeled_df)),
        "label_counts": {str(key): int(value) for key, value in label_counts.items()},
        "matched_negative_eligible_count": matched_count,
        "n0_only_for_empirical_null": True,
    }


def main():
    args = parse_args()
    logger = setup_logger("negative_bank", args.log or None)
    samples_df = read_table(args.input_samples, required_columns=["sample_id"], empty_ok=True)
    labeled = assign_negative_bank_labels(samples_df, version=args.version)
    write_table(args.output_labels, labeled)
    write_json(args.output_summary, summarize_negative_bank(labeled, version=args.version))
    logger.info("wrote negative-bank labels rows=%d to %s", len(labeled), args.output_labels)


if __name__ == "__main__":
    main()
