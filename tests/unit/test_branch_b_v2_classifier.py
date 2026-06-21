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
    assert row["v2_candidate_class"] == "V2_POSITIVE_SUPPORT_REVIEW"
    assert row["v2_classifier_action"] == "V2_REVIEW_POSITIVE_SUPPORT"
    assert row["v2_evidence_tier"] == "UNKNOWN_BACKGROUND_POSITIVE_SUPPORT"


def test_phase1_unknown_matched_negative_source_forces_review_not_artifact():
    frame = pd.DataFrame(
        [
            {
                "sample_id": "H2",
                "candidate_id": "H2_chr6_tail_gain",
                "chrom": "chr6",
                "start": 168000001,
                "end": 170000000,
                "state": "gain",
                "a_abs_zscore": 17.0,
                "branch_b_report_class": "candidate_suppressed",
                "branch_b_artifact_status": "artifact",
                "branch_b_keep_event": 0,
                "final_disposition": "LIKELY_ARTIFACT",
                "matched_negative_source": "UNKNOWN_BACKGROUND",
                "matched_negative_percentile": pd.NA,
            }
        ]
    )

    row = classify_branch_b_v2_candidates(frame).iloc[0]

    assert row["final_disposition"] == "LIKELY_ARTIFACT"
    assert row["matched_negative_source"] == "UNKNOWN_BACKGROUND"
    assert row["v2_candidate_class"] == "V2_POSITIVE_SUPPORT_REVIEW"
    assert row["v2_classifier_action"] == "V2_REVIEW_POSITIVE_SUPPORT"
    assert row["v2_evidence_tier"] == "UNKNOWN_BACKGROUND_POSITIVE_SUPPORT"


def test_blank_background_status_falls_back_to_unknown_source():
    frame = pd.DataFrame(
        [
            {
                "sample_id": "H2",
                "candidate_id": "H2_chr21_gain",
                "chrom": "chr21",
                "start": 20_000_000,
                "end": 42_000_000,
                "state": "gain",
                "final_disposition": "LIKELY_ARTIFACT",
                "matched_negative_background_status": pd.NA,
                "matched_negative_source": "UNKNOWN_BACKGROUND",
            }
        ]
    )

    row = classify_branch_b_v2_candidates(frame).iloc[0]

    assert row["v2_candidate_class"] == "V2_REVIEW_REQUIRED"
    assert row["v2_classifier_action"] == "V2_REVIEW_ONLY"
    assert row["v2_evidence_tier"] == "UNKNOWN_BACKGROUND_REVIEW"


def test_legacy_final_disposition_does_not_change_v2_classification():
    base = {
        "sample_id": "Y1",
        "chrom": "chr21",
        "start": 1000000,
        "end": 6000000,
        "state": "gain",
        "a_abs_zscore": 18.0,
        "same_direction_fraction": 0.95,
        "corrected_amplitude": 5.1,
        "matched_negative_background_status": "OK",
        "matched_negative_abs_percentile": 0.995,
    }
    frame = pd.DataFrame(
        [
            {
                **base,
                "candidate_id": "legacy_confirmed",
                "final_disposition": "CONFIRMED",
                "branch_b_keep_event": 1,
            },
            {
                **base,
                "candidate_id": "legacy_artifact",
                "final_disposition": "LIKELY_ARTIFACT",
                "branch_b_keep_event": 0,
            },
        ]
    )
    classified = classify_branch_b_v2_candidates(frame).set_index("candidate_id")

    assert classified.loc["legacy_confirmed", "v2_candidate_class"] == classified.loc["legacy_artifact", "v2_candidate_class"]
    assert classified.loc["legacy_confirmed", "v2_classifier_action"] == classified.loc["legacy_artifact", "v2_classifier_action"]
    assert classified.loc["legacy_confirmed", "v2_evidence_tier"] == "MATCHED_NEGATIVE_OUTLIER_POSITIVE_SUPPORT"


def test_v2_summary_reports_review_burden_without_report_promotion():
    classified = classify_branch_b_v2_candidates(_candidate_frame())
    summary = summarize_v2_classification("Y1", classified)

    assert summary["sample_id"] == "Y1"
    assert summary["candidate_count"] == 2
    assert summary["class_counts"] == {"V2_POSITIVE_SUPPORT_REVIEW": 1, "V2_REVIEW_REQUIRED": 1}
    assert summary["final_report_impact"] == "none_shadow_only"
    assert summary["legacy_decision_fields_ignored"] is True


def test_unknown_background_strong_a_support_gets_positive_review_tier():
    frame = pd.DataFrame(
        [
            {
                "sample_id": "H6",
                "candidate_id": "H6_chr21_gain",
                "chrom": "chr21",
                "start": 15_000_001,
                "end": 42_000_000,
                "state": "gain",
                "a_abs_zscore": 7.11,
                "same_direction_fraction": 0.82,
                "corrected_amplitude": 2.4,
                "final_disposition": "LIKELY_ARTIFACT",
                "matched_negative_background_status": "UNKNOWN_BACKGROUND",
                "matched_negative_action": "REVIEW_NO_CALL",
            }
        ]
    )

    row = classify_branch_b_v2_candidates(frame).iloc[0]

    assert row["v2_candidate_class"] == "V2_POSITIVE_SUPPORT_REVIEW"
    assert row["v2_classifier_action"] == "V2_REVIEW_POSITIVE_SUPPORT"
    assert row["v2_final_report_impact"] == "none_shadow_only"
    assert row["v2_evidence_tier"] == "UNKNOWN_BACKGROUND_POSITIVE_SUPPORT"
    assert row["v2_review_priority"] == "HIGH"


def test_matched_negative_outlier_gets_positive_support_tier_shadow_only():
    frame = pd.DataFrame(
        [
            {
                "sample_id": "Y2",
                "candidate_id": "Y2_chr21_gain",
                "chrom": "chr21",
                "start": 20_000_000,
                "end": 42_000_000,
                "state": "gain",
                "a_abs_zscore": 15.0,
                "same_direction_fraction": 0.95,
                "corrected_amplitude": 5.1,
                "final_disposition": "REVIEW_REQUIRED",
                "matched_negative_background_status": "OK",
                "matched_negative_abs_percentile": 0.995,
                "matched_negative_action": "BACKGROUND_SUPPORTED",
            }
        ]
    )

    row = classify_branch_b_v2_candidates(frame).iloc[0]

    assert row["v2_candidate_class"] == "V2_POSITIVE_SUPPORT_REVIEW"
    assert row["v2_evidence_tier"] == "MATCHED_NEGATIVE_OUTLIER_POSITIVE_SUPPORT"
    assert row["v2_evidence_gate"] == "BACKGROUND_INFORMATIVE"
    assert row["v2_review_priority"] == "HIGH"


def test_v2_summary_reports_evidence_tier_counts():
    classified = classify_branch_b_v2_candidates(
        pd.DataFrame(
            [
                {
                    "sample_id": "H6",
                    "candidate_id": "H6_chr21_gain",
                    "chrom": "chr21",
                    "start": 15_000_001,
                    "end": 42_000_000,
                    "state": "gain",
                    "a_abs_zscore": 7.11,
                    "same_direction_fraction": 0.82,
                    "corrected_amplitude": 2.4,
                    "final_disposition": "LIKELY_ARTIFACT",
                    "matched_negative_background_status": "UNKNOWN_BACKGROUND",
                },
                {
                    "sample_id": "Y2",
                    "candidate_id": "Y2_chr21_gain",
                    "chrom": "chr21",
                    "start": 20_000_000,
                    "end": 42_000_000,
                    "state": "gain",
                    "a_abs_zscore": 15.0,
                    "same_direction_fraction": 0.95,
                    "corrected_amplitude": 5.1,
                    "final_disposition": "REVIEW_REQUIRED",
                    "matched_negative_background_status": "OK",
                    "matched_negative_abs_percentile": 0.995,
                },
            ]
        )
    )
    summary = summarize_v2_classification("mixed", classified)

    assert summary["evidence_tier_counts"] == {
        "MATCHED_NEGATIVE_OUTLIER_POSITIVE_SUPPORT": 1,
        "UNKNOWN_BACKGROUND_POSITIVE_SUPPORT": 1,
    }
    assert summary["review_priority_counts"] == {"HIGH": 2}


def test_shadow_background_context_does_not_become_final_filter():
    frame = pd.DataFrame(
        [
            {
                "sample_id": "JZ26125843-56-56",
                "candidate_id": "JZ26125843-56-56.A0001",
                "chrom": "chr19",
                "start": 10_000_000,
                "end": 20_000_000,
                "state": "gain",
                "a_abs_zscore": 8.0,
                "same_direction_fraction": 0.90,
                "corrected_amplitude": 2.0,
                "final_disposition": "REVIEW_REQUIRED",
                "matched_negative_background_status": "SHADOW_BACKGROUND",
                "matched_negative_source": "UNKNOWN_BACKGROUND",
                "matched_negative_abs_percentile": 0.99,
                "matched_negative_action": "SHADOW_CONTEXT_ONLY",
            }
        ]
    )

    row = classify_branch_b_v2_candidates(frame).iloc[0]

    assert row["v2_candidate_class"] == "V2_POSITIVE_SUPPORT_REVIEW"
    assert row["v2_classifier_action"] == "V2_REVIEW_POSITIVE_SUPPORT"
    assert row["v2_evidence_tier"] == "SHADOW_BACKGROUND_OUTLIER_POSITIVE_SUPPORT"
    assert row["v2_evidence_gate"] == "SHADOW_BACKGROUND_CONTEXT"
    assert row["v2_final_report_impact"] == "none_shadow_only"
