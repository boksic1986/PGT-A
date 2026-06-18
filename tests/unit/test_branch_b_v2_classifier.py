from __future__ import annotations

import pandas as pd

from pgta.predict.branch_b.v2_classifier import (
    classify_branch_b_v2_candidates,
    summarize_v2_classification,
)


def _candidate_frame():
    return pd.DataFrame(
        [
            {
                "sample_id": "Y1",
                "candidate_id": "Y1_chr21_gain",
                "chrom": "chr21",
                "start": 1000000,
                "end": 6000000,
                "state": "gain",
                "a_abs_zscore": 18.0,
                "branch_b_report_class": "candidate_pass",
                "branch_b_artifact_status": "pass",
                "branch_b_keep_event": 1,
                "final_disposition": "CONFIRMED",
                "matched_negative_background_status": "OK",
                "matched_negative_action": "BACKGROUND_SUPPORTED",
            },
            {
                "sample_id": "Y1",
                "candidate_id": "Y1_chr16_loss",
                "chrom": "chr16",
                "start": 500000,
                "end": 5000000,
                "state": "loss",
                "a_abs_zscore": 29.0,
                "branch_b_report_class": "candidate_suppressed",
                "branch_b_artifact_status": "artifact",
                "branch_b_keep_event": 0,
                "final_disposition": "LIKELY_ARTIFACT",
                "matched_negative_background_status": "UNKNOWN_BACKGROUND",
                "matched_negative_action": "REVIEW_NO_CALL",
            },
        ]
    )


def test_v2_classifier_keeps_one_row_per_a_candidate_and_shadow_only():
    classified = classify_branch_b_v2_candidates(_candidate_frame())

    assert classified["candidate_id"].tolist() == ["Y1_chr21_gain", "Y1_chr16_loss"]
    assert classified["v2_final_report_impact"].tolist() == ["none_shadow_only", "none_shadow_only"]
    assert "final_disposition" in classified.columns
    assert "branch_b_report_class" in classified.columns


def test_unknown_matched_negative_background_forces_review_not_artifact():
    classified = classify_branch_b_v2_candidates(_candidate_frame())
    row = classified.loc[classified["candidate_id"].eq("Y1_chr16_loss")].iloc[0]

    assert row["final_disposition"] == "LIKELY_ARTIFACT"
    assert row["matched_negative_background_status"] == "UNKNOWN_BACKGROUND"
    assert row["v2_candidate_class"] == "REVIEW_REQUIRED"
    assert row["v2_classifier_action"] == "SHADOW_REVIEW_ONLY"
    assert "unknown_matched_negative_background" in row["v2_classifier_reason"]


def test_confirmed_legacy_candidate_is_confirmed_shadow():
    classified = classify_branch_b_v2_candidates(_candidate_frame())
    row = classified.loc[classified["candidate_id"].eq("Y1_chr21_gain")].iloc[0]

    assert row["v2_candidate_class"] == "CONFIRMED_SHADOW"
    assert row["v2_classifier_action"] == "SHADOW_CONFIRM"
    assert row["v2_classifier_reason"] == "legacy_confirmed_with_shadow_evidence"


def test_v2_summary_reports_review_burden_without_report_promotion():
    classified = classify_branch_b_v2_candidates(_candidate_frame())
    summary = summarize_v2_classification("Y1", classified)

    assert summary["sample_id"] == "Y1"
    assert summary["candidate_count"] == 2
    assert summary["class_counts"] == {"CONFIRMED_SHADOW": 1, "REVIEW_REQUIRED": 1}
    assert summary["final_report_impact"] == "none_shadow_only"
