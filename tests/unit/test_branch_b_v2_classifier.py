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
    assert "v2_filter_action" in classified.columns
    assert "v2_filter_reason" in classified.columns
    assert "v2_filter_scope" in classified.columns
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


def test_sex_chromosome_candidate_routes_to_branch_s_review_not_autosomal_positive_support():
    frame = pd.DataFrame(
        [
            {
                "sample_id": "H5",
                "candidate_id": "H5_chrX_loss",
                "chrom": "chrX",
                "start": 1,
                "end": 155_000_000,
                "state": "loss",
                "a_abs_zscore": 62.0,
                "same_direction_fraction": 0.80,
                "corrected_amplitude": -2.1,
                "matched_negative_background_status": "UNKNOWN_BACKGROUND",
                "matched_negative_action": "REVIEW_NO_CALL",
            }
        ]
    )

    row = classify_branch_b_v2_candidates(frame).iloc[0]

    assert row["v2_candidate_class"] == "V2_SEX_CHROMOSOME_REVIEW"
    assert row["v2_classifier_action"] == "V2_ROUTE_BRANCH_S_REVIEW"
    assert row["v2_classifier_reason"] == "sex_chromosome_branch_s_review:unknown_background_positive_support"
    assert row["v2_evidence_tier"] == "UNKNOWN_BACKGROUND_POSITIVE_SUPPORT"
    assert row["v2_final_report_impact"] == "none_shadow_only"


def test_b_side_weak_signal_is_review_label_only_for_a_anchored_truth_like_positive():
    frame = pd.DataFrame(
        [
            {
                "sample_id": "Y2",
                "candidate_id": "Y2_chr14_gain",
                "chrom": "chr14",
                "start": 1,
                "end": 107_000_000,
                "state": "gain",
                "a_abs_zscore": 94.8,
                "same_direction_fraction": 0.15,
                "corrected_amplitude": -0.29,
                "raw_amplitude": 0.12,
                "matched_negative_background_status": "UNKNOWN_BACKGROUND",
            }
        ]
    )

    row = classify_branch_b_v2_candidates(frame).iloc[0]

    assert row["v2_candidate_class"] == "V2_POSITIVE_SUPPORT_REVIEW"
    assert row["v2_classifier_action"] == "V2_REVIEW_POSITIVE_SUPPORT"
    assert row["v2_final_report_impact"] == "none_shadow_only"
    assert row["v2_direction_support_label"] == "A_ANCHORED_WEAK_B_SIGNAL"
    assert row["v2_b_signal_context_label"] == "A_ANCHORED_WEAK_B_SIGNAL"
    assert row["v2_b_signal_context_reason"] == "positive_support_without_branch_b_signal_support"


def test_b_side_signal_support_label_marks_branch_b_supported_without_report_promotion():
    frame = pd.DataFrame(
        [
            {
                "sample_id": "H6",
                "candidate_id": "H6_chr21_gain",
                "chrom": "chr21",
                "start": 20_700_000,
                "end": 22_300_000,
                "state": "gain",
                "a_abs_zscore": 7.11,
                "same_direction_fraction": 0.92,
                "corrected_amplitude": 0.40,
                "raw_amplitude": 0.35,
                "matched_negative_background_status": "UNKNOWN_BACKGROUND",
            }
        ]
    )

    row = classify_branch_b_v2_candidates(frame).iloc[0]

    assert row["v2_candidate_class"] == "V2_POSITIVE_SUPPORT_REVIEW"
    assert row["v2_final_report_impact"] == "none_shadow_only"
    assert row["v2_direction_support_label"] == "B_SIGNAL_SUPPORTED_A_DIRECTION"
    assert row["v2_b_signal_context_label"] == "B_SIGNAL_SUPPORTED_A_DIRECTION"
    assert row["v2_b_signal_context_reason"] == "same_direction_fraction_ge_0.50"
    assert row["v2_signal_strength_tier"] == "A_SENSITIVE_Z_5_TO_10"
    assert row["v2_length_tier"] == "review_only_ge1mb"
    assert row["v2_disposition"] == "background_unknown_review"
    assert row["v2_filter_action"] == "keep_background_unknown_review"
    assert row["v2_filter_hard_suppression_allowed"] == 0


def test_v2_summary_reports_b_signal_context_label_counts():
    classified = classify_branch_b_v2_candidates(
        pd.DataFrame(
            [
                {
                    "sample_id": "Y2",
                    "candidate_id": "Y2_chr14_gain",
                    "chrom": "chr14",
                    "state": "gain",
                    "a_abs_zscore": 94.8,
                    "same_direction_fraction": 0.15,
                    "corrected_amplitude": -0.29,
                    "matched_negative_background_status": "UNKNOWN_BACKGROUND",
                },
                {
                    "sample_id": "H6",
                    "candidate_id": "H6_chr21_gain",
                    "chrom": "chr21",
                    "state": "gain",
                    "a_abs_zscore": 7.11,
                    "same_direction_fraction": 0.92,
                    "corrected_amplitude": 0.40,
                    "matched_negative_background_status": "UNKNOWN_BACKGROUND",
                },
            ]
        )
    )

    summary = summarize_v2_classification("mixed", classified)

    assert summary["b_signal_context_label_counts"] == {
        "A_ANCHORED_WEAK_B_SIGNAL": 1,
        "B_SIGNAL_SUPPORTED_A_DIRECTION": 1,
    }
    assert summary["disposition_counts"] == {
        "background_unknown_review": 2,
    }


def test_b_side_signal_discordance_is_review_context_not_no_call_contract_risk():
    frame = pd.DataFrame(
        [
            {
                "sample_id": "Y1",
                "candidate_id": "Y1_chr16_loss",
                "chrom": "chr16",
                "state": "loss",
                "a_abs_zscore": 29.0,
                "same_direction_fraction": 0.80,
                "corrected_amplitude": 2.5,
                "matched_negative_background_status": "UNKNOWN_BACKGROUND",
            }
        ]
    )

    row = classify_branch_b_v2_candidates(frame).iloc[0]

    assert row["v2_candidate_class"] == "V2_POSITIVE_SUPPORT_REVIEW"
    assert row["v2_classifier_action"] == "V2_REVIEW_POSITIVE_SUPPORT"
    assert row["v2_direction_support_label"] == "B_SIGNAL_DISCORDANT_WITH_A_DIRECTION"
    assert row["v2_b_signal_context_label"] == "B_SIGNAL_DISCORDANT_WITH_A_DIRECTION"
    assert row["v2_b_signal_context_reason"] == "b_side_amplitude_opposite_a_direction_abs_ge_2"
    assert row["v2_evidence_tier"] == "UNKNOWN_BACKGROUND_POSITIVE_SUPPORT"
    assert row["v2_disposition"] == "background_unknown_review"
    assert row["v2_filter_action"] == "keep_background_unknown_review"
    assert row["v2_filter_hard_suppression_allowed"] == 0
    assert row["v2_final_report_impact"] == "none_shadow_only"


def test_length_and_clean_support_tiers_do_not_hard_suppress_sensitive_positive():
    frame = pd.DataFrame(
        [
            {
                "sample_id": "H6",
                "candidate_id": "H6_chr21_gain",
                "chrom": "chr21",
                "start": 20_700_000,
                "end": 22_300_000,
                "state": "gain",
                "a_abs_zscore": 7.11,
                "clean_bin_fraction": 0.85,
                "high_risk_bin_fraction": 0.05,
                "attenuation_ratio": 0.42,
                "same_direction_fraction": 0.92,
                "corrected_amplitude": 0.40,
                "raw_amplitude": 0.95,
                "matched_negative_background_status": "UNKNOWN_BACKGROUND",
                "calibration_null_status": "NO_NULL_SUPPORT",
            }
        ]
    )

    row = classify_branch_b_v2_candidates(frame).iloc[0]

    assert row["v2_signal_strength_tier"] == "A_SENSITIVE_Z_5_TO_10"
    assert row["v2_length_tier"] == "review_only_ge1mb"
    assert row["v2_clean_support_label"] == "CLEAN_SUPPORT_AVAILABLE"
    assert row["v2_gc_rc_context_label"] == "GC_RC_ATTENUATED_SEVERE"
    assert row["v2_candidate_class"] == "V2_POSITIVE_SUPPORT_REVIEW"
    assert row["v2_disposition"] == "background_unknown_review"
    assert row["v2_filter_action"] == "keep_background_unknown_review"
    assert row["v2_filter_hard_suppression_allowed"] == 0


def test_unknown_background_no_null_support_is_explicit_review_context_not_filter():
    frame = pd.DataFrame(
        [
            {
                "sample_id": "H6",
                "candidate_id": "H6_chr21_gain",
                "chrom": "chr21",
                "state": "gain",
                "a_abs_zscore": 7.11,
                "same_direction_fraction": 0.92,
                "corrected_amplitude": 0.40,
                "matched_negative_background_status": "UNKNOWN_BACKGROUND",
                "calibration_null_status": "NO_NULL_SUPPORT",
            }
        ]
    )

    row = classify_branch_b_v2_candidates(frame).iloc[0]

    assert row["v2_candidate_class"] == "V2_POSITIVE_SUPPORT_REVIEW"
    assert row["v2_classifier_action"] == "V2_REVIEW_POSITIVE_SUPPORT"
    assert row["v2_final_report_impact"] == "none_shadow_only"
    assert row["v2_background_context_label"] == "UNKNOWN_BACKGROUND_NO_NULL_SUPPORT"
    assert row["v2_background_context_reason"] == "no_matched_negative_and_no_calibration_null"


def test_shadow_background_no_null_support_is_limited_context_only():
    frame = pd.DataFrame(
        [
            {
                "sample_id": "Y2",
                "candidate_id": "Y2_chr21_gain",
                "chrom": "chr21",
                "state": "gain",
                "a_abs_zscore": 15.0,
                "same_direction_fraction": 0.95,
                "corrected_amplitude": 5.1,
                "matched_negative_background_status": "SHADOW_BACKGROUND",
                "matched_negative_abs_percentile": 0.995,
                "calibration_null_status": "NO_NULL_SUPPORT",
            }
        ]
    )

    row = classify_branch_b_v2_candidates(frame).iloc[0]

    assert row["v2_candidate_class"] == "V2_POSITIVE_SUPPORT_REVIEW"
    assert row["v2_evidence_gate"] == "SHADOW_BACKGROUND_CONTEXT"
    assert row["v2_final_report_impact"] == "none_shadow_only"
    assert row["v2_background_context_label"] == "SHADOW_BACKGROUND_NO_NULL_SUPPORT"
    assert row["v2_background_context_reason"] == "shadow_background_context_only_no_calibration_null"


def test_v2_summary_reports_background_context_label_counts():
    classified = classify_branch_b_v2_candidates(
        pd.DataFrame(
            [
                {
                    "sample_id": "H6",
                    "candidate_id": "H6_chr21_gain",
                    "chrom": "chr21",
                    "state": "gain",
                    "a_abs_zscore": 7.11,
                    "same_direction_fraction": 0.92,
                    "corrected_amplitude": 0.40,
                    "matched_negative_background_status": "UNKNOWN_BACKGROUND",
                    "calibration_null_status": "NO_NULL_SUPPORT",
                },
                {
                    "sample_id": "Y2",
                    "candidate_id": "Y2_chr21_gain",
                    "chrom": "chr21",
                    "state": "gain",
                    "a_abs_zscore": 15.0,
                    "same_direction_fraction": 0.95,
                    "corrected_amplitude": 5.1,
                    "matched_negative_background_status": "SHADOW_BACKGROUND",
                    "matched_negative_abs_percentile": 0.995,
                    "calibration_null_status": "NO_NULL_SUPPORT",
                },
            ]
        )
    )

    summary = summarize_v2_classification("mixed", classified)

    assert summary["background_context_label_counts"] == {
        "SHADOW_BACKGROUND_NO_NULL_SUPPORT": 1,
        "UNKNOWN_BACKGROUND_NO_NULL_SUPPORT": 1,
    }


def test_ref_contract_risk_is_the_only_truth_safe_hard_filter_action():
    frame = pd.DataFrame(
        [
            {
                "sample_id": "Y1",
                "candidate_id": "Y1_chr1_gain_leakage",
                "chrom": "chr1",
                "start": 1_000_000,
                "end": 12_000_000,
                "state": "gain",
                "a_abs_zscore": 32.0,
                "same_direction_fraction": 0.95,
                "corrected_amplitude": 4.0,
                "matched_negative_background_status": "UNKNOWN_BACKGROUND",
                "refmap_status": "SAME_CHROM_REF_LEAKAGE",
            }
        ]
    )

    row = classify_branch_b_v2_candidates(frame).iloc[0]

    assert row["v2_candidate_class"] == "V2_NO_CALL_CONTRACT_RISK"
    assert row["v2_disposition"] == "technical_risk_review"
    assert row["v2_filter_action"] == "suppress_workflow_contract_risk"
    assert row["v2_filter_reason"] == "same_chrom_ref_leakage_contract_risk"
    assert row["v2_filter_hard_suppression_allowed"] == 1
    assert row["v2_final_report_impact"] == "none_shadow_only"


def test_low_clean_support_high_risk_downgrades_but_does_not_hard_filter():
    frame = pd.DataFrame(
        [
            {
                "sample_id": "JZ26125843-56-56",
                "candidate_id": "JZ26125843-56-56.A0001",
                "chrom": "chr16",
                "start": 500_000,
                "end": 4_500_000,
                "state": "loss",
                "a_abs_zscore": 6.0,
                "same_direction_fraction": 0.55,
                "corrected_amplitude": -0.8,
                "clean_bin_fraction": 0.10,
                "high_risk_bin_fraction": 0.80,
                "hard_region_fraction": 0.65,
                "matched_negative_background_status": "UNKNOWN_BACKGROUND",
                "calibration_null_status": "LIMITED_NULL_SUPPORT",
            }
        ]
    )

    row = classify_branch_b_v2_candidates(frame).iloc[0]

    assert row["v2_clean_support_label"] == "LOW_CLEAN_SUPPORT_HIGH_RISK"
    assert row["v2_filter_action"] == "downgrade_to_technical_risk_review"
    assert row["v2_filter_hard_suppression_allowed"] == 0
    assert row["v2_final_report_impact"] == "none_shadow_only"


def test_v2_summary_reports_filter_action_counts():
    classified = classify_branch_b_v2_candidates(
        pd.DataFrame(
            [
                {
                    "sample_id": "H6",
                    "candidate_id": "H6_chr21_gain",
                    "chrom": "chr21",
                    "state": "gain",
                    "a_abs_zscore": 7.11,
                    "same_direction_fraction": 0.92,
                    "matched_negative_background_status": "UNKNOWN_BACKGROUND",
                },
                {
                    "sample_id": "Y1",
                    "candidate_id": "Y1_chr1_gain_leakage",
                    "chrom": "chr1",
                    "state": "gain",
                    "a_abs_zscore": 32.0,
                    "refmap_status": "SAME_CHROM_REF_LEAKAGE",
                    "matched_negative_background_status": "UNKNOWN_BACKGROUND",
                },
            ]
        )
    )

    summary = summarize_v2_classification("mixed", classified)

    assert summary["filter_action_counts"] == {
        "keep_background_unknown_review": 1,
        "suppress_workflow_contract_risk": 1,
    }
    assert summary["filter_hard_suppression_allowed_count"] == 1


def test_burden_stratification_keeps_h6_chr21_review_and_marks_cnvpro_sources():
    frame = pd.DataFrame(
        [
            {
                "sample_id": "H6",
                "candidate_id": "H6_chr21_gain",
                "chrom": "chr21",
                "start": 20_700_000,
                "end": 22_300_000,
                "state": "gain",
                "a_abs_zscore": 7.11,
                "clean_bin_fraction": 0.85,
                "high_risk_bin_fraction": 0.05,
                "attenuation_ratio": 0.42,
                "same_direction_fraction": 0.92,
                "corrected_amplitude": 0.40,
                "raw_amplitude": 0.95,
                "matched_negative_background_status": "UNKNOWN_BACKGROUND",
                "calibration_null_status": "NO_NULL_SUPPORT",
            }
        ]
    )

    row = classify_branch_b_v2_candidates(frame).iloc[0]

    assert row["v2_candidate_class"] == "V2_POSITIVE_SUPPORT_REVIEW"
    assert row["v2_burden_reduction_tier"] == "background_unknown_review"
    assert row["v2_burden_reduction_action"] == "stratify_background_unknown_review"
    assert row["v2_filter_hard_suppression_allowed"] == 0
    tags = row["v2_burden_evidence_tags"]
    assert "[CNVpro-inspired] length_tier=review_only_ge1mb" in tags
    assert "[CNVpro-confirmed] acrocentric_qter_context_review_only" in tags
    assert "[CNVpro-like] gc_rc_context=GC_RC_ATTENUATED_SEVERE" in tags
    assert "[Not used] CNVcalling_R_cghFLasso_not_primary_caller" in tags


def test_burden_stratification_routes_sex_chromosome_to_branch_s_with_cnvseq_asset_tag():
    frame = pd.DataFrame(
        [
            {
                "sample_id": "H5",
                "candidate_id": "H5_chrX_loss",
                "chrom": "chrX",
                "start": 1,
                "end": 155_000_000,
                "state": "loss",
                "a_abs_zscore": 62.0,
                "same_direction_fraction": 0.80,
                "corrected_amplitude": -2.1,
                "matched_negative_background_status": "UNKNOWN_BACKGROUND",
                "calibration_null_status": "NO_NULL_SUPPORT",
            }
        ]
    )

    row = classify_branch_b_v2_candidates(frame).iloc[0]

    assert row["v2_candidate_class"] == "V2_SEX_CHROMOSOME_REVIEW"
    assert row["v2_burden_reduction_tier"] == "branch_s_review"
    assert row["v2_burden_reduction_action"] == "route_to_branch_s_review"
    assert "[CNVseq-asset] sex_homology_PAR_annotation_branch_s_context" in row["v2_burden_evidence_tags"]


def test_v2_summary_reports_burden_reduction_counts_and_tag_counts():
    classified = classify_branch_b_v2_candidates(
        pd.DataFrame(
            [
                {
                    "sample_id": "H6",
                    "candidate_id": "H6_chr21_gain",
                    "chrom": "chr21",
                    "start": 20_700_000,
                    "end": 22_300_000,
                    "state": "gain",
                    "a_abs_zscore": 7.11,
                    "same_direction_fraction": 0.92,
                    "attenuation_ratio": 0.42,
                    "matched_negative_background_status": "UNKNOWN_BACKGROUND",
                    "calibration_null_status": "NO_NULL_SUPPORT",
                },
                {
                    "sample_id": "H5",
                    "candidate_id": "H5_chrX_loss",
                    "chrom": "chrX",
                    "start": 1,
                    "end": 155_000_000,
                    "state": "loss",
                    "a_abs_zscore": 62.0,
                    "matched_negative_background_status": "UNKNOWN_BACKGROUND",
                },
            ]
        )
    )

    summary = summarize_v2_classification("mixed", classified)

    assert summary["burden_reduction_tier_counts"] == {
        "background_unknown_review": 1,
        "branch_s_review": 1,
    }
    assert summary["burden_reduction_action_counts"] == {
        "route_to_branch_s_review": 1,
        "stratify_background_unknown_review": 1,
    }
    assert summary["length_tier_counts"]["review_only_ge1mb"] == 1
    assert summary["gc_rc_context_label_counts"]["GC_RC_ATTENUATED_SEVERE"] == 1
    assert summary["burden_evidence_tag_counts"]["[CNVpro-inspired] length_tier=review_only_ge1mb"] == 1
    assert summary["burden_evidence_tag_counts"]["[CNVseq-asset] sex_homology_PAR_annotation_branch_s_context"] == 1


def test_report_layer_filters_only_combined_technical_risk_not_single_indicators():
    frame = pd.DataFrame(
        [
            {
                "sample_id": "S1",
                "candidate_id": "background_unknown_only",
                "chrom": "chr1",
                "start": 1_000_000,
                "end": 4_000_000,
                "state": "gain",
                "a_abs_zscore": 6.5,
                "same_direction_fraction": 0.90,
                "clean_bin_fraction": 0.80,
                "high_risk_bin_fraction": 0.05,
                "matched_negative_background_status": "UNKNOWN_BACKGROUND",
            },
            {
                "sample_id": "S1",
                "candidate_id": "b_discordance_only",
                "chrom": "chr2",
                "start": 1_000_000,
                "end": 4_000_000,
                "state": "gain",
                "a_abs_zscore": 6.5,
                "same_direction_fraction": 0.20,
                "corrected_amplitude": -2.4,
                "clean_bin_fraction": 0.80,
                "high_risk_bin_fraction": 0.05,
                "matched_negative_background_status": "UNKNOWN_BACKGROUND",
            },
            {
                "sample_id": "S1",
                "candidate_id": "gc_attenuation_only",
                "chrom": "chr3",
                "start": 1_000_000,
                "end": 4_000_000,
                "state": "loss",
                "a_abs_zscore": 6.5,
                "same_direction_fraction": 0.90,
                "attenuation_ratio": 0.35,
                "clean_bin_fraction": 0.80,
                "high_risk_bin_fraction": 0.05,
                "matched_negative_background_status": "UNKNOWN_BACKGROUND",
            },
            {
                "sample_id": "S1",
                "candidate_id": "short_length_only",
                "chrom": "chr4",
                "start": 1_000_000,
                "end": 1_800_000,
                "state": "gain",
                "a_abs_zscore": 6.5,
                "same_direction_fraction": 0.90,
                "clean_bin_fraction": 0.80,
                "high_risk_bin_fraction": 0.05,
                "matched_negative_background_status": "UNKNOWN_BACKGROUND",
            },
        ]
    )

    classified = classify_branch_b_v2_candidates(frame).set_index("candidate_id")

    assert set(classified["v2_report_layer_class"]) == {"internal_review_event"}
    assert set(classified["v2_report_visibility"]) == {"internal_review"}
    assert "filtered_event" not in set(classified["v2_report_layer_class"])


def test_report_layer_combined_technical_risk_filters_to_audit_only():
    frame = pd.DataFrame(
        [
            {
                "sample_id": "JZ26125843-56-56",
                "candidate_id": "combined_risk_chr16_loss",
                "chrom": "chr16",
                "start": 500_000,
                "end": 1_800_000,
                "state": "loss",
                "a_abs_zscore": 6.1,
                "same_direction_fraction": 0.10,
                "corrected_amplitude": 2.6,
                "attenuation_ratio": 0.35,
                "clean_bin_fraction": 0.10,
                "high_risk_bin_fraction": 0.85,
                "hard_region_fraction": 0.70,
                "matched_negative_background_status": "UNKNOWN_BACKGROUND",
                "calibration_null_status": "NO_NULL_SUPPORT",
            }
        ]
    )

    row = classify_branch_b_v2_candidates(frame).iloc[0]

    assert row["v2_report_layer_class"] == "filtered_event"
    assert row["v2_report_visibility"] == "audit_only"
    assert row["v2_filter_action"] == "filter_report_layer_combined_technical_risk"
    assert row["v2_filter_hard_suppression_allowed"] == 0
    assert "short_or_focal" in row["v2_report_filter_rule_tags"]
    assert "low_clean_high_risk" in row["v2_report_filter_rule_tags"]
    assert "b_signal_not_supportive" in row["v2_report_filter_rule_tags"]
    assert "gc_rc_attenuated" in row["v2_report_filter_rule_tags"]


def test_report_layer_filters_sensitive_b_unsupported_gc_attenuated_candidate():
    frame = pd.DataFrame(
        [
            {
                "sample_id": "JZ26125843-56-56",
                "candidate_id": "sensitive_b_unsupported_gc_attenuated",
                "chrom": "chr17",
                "start": 2_250_001,
                "end": 21_000_000,
                "state": "loss",
                "a_abs_zscore": 6.2,
                "same_direction_fraction": 0.20,
                "corrected_amplitude": 0.20,
                "attenuation_ratio": 0.35,
                "clean_bin_fraction": 0.85,
                "high_risk_bin_fraction": 0.05,
                "matched_negative_background_status": "UNKNOWN_BACKGROUND",
                "calibration_null_status": "NO_NULL_SUPPORT",
            }
        ]
    )

    row = classify_branch_b_v2_candidates(frame).iloc[0]

    assert row["v2_signal_strength_tier"] == "A_SENSITIVE_Z_5_TO_10"
    assert row["v2_report_layer_class"] == "filtered_event"
    assert row["v2_report_visibility"] == "audit_only"
    assert row["v2_filter_action"] == "filter_report_layer_combined_technical_risk"
    assert row["v2_filter_hard_suppression_allowed"] == 0
    assert "b_signal_not_supportive" in row["v2_report_filter_rule_tags"]
    assert "gc_rc_attenuated" in row["v2_report_filter_rule_tags"]


def test_report_layer_promotes_strong_supported_autosomal_candidate_to_report_event():
    frame = pd.DataFrame(
        [
            {
                "sample_id": "Y1",
                "candidate_id": "Y1_chr21_loss",
                "chrom": "chr21",
                "start": 1_000_000,
                "end": 42_000_000,
                "state": "loss",
                "a_abs_zscore": 112.0,
                "same_direction_fraction": 0.92,
                "attenuation_ratio": 0.72,
                "clean_bin_fraction": 0.85,
                "high_risk_bin_fraction": 0.05,
                "matched_negative_background_status": "UNKNOWN_BACKGROUND",
                "calibration_null_status": "NO_NULL_SUPPORT",
            }
        ]
    )

    row = classify_branch_b_v2_candidates(frame).iloc[0]

    assert row["v2_report_layer_class"] == "report_event"
    assert row["v2_report_visibility"] == "final_report"
    assert row["v2_report_filter_reason"] == "strong_a_supported_report_layer_event"


def test_report_layer_strong_a_b_unsupported_gc_attenuated_stays_review_not_filtered():
    frame = pd.DataFrame(
        [
            {
                "sample_id": "Y2",
                "candidate_id": "Y2_chr14_gain",
                "chrom": "chr14",
                "start": 1,
                "end": 107_000_000,
                "state": "gain",
                "a_abs_zscore": 94.8,
                "same_direction_fraction": 0.15,
                "corrected_amplitude": -0.29,
                "attenuation_ratio": 0.72,
                "clean_bin_fraction": 0.85,
                "high_risk_bin_fraction": 0.05,
                "matched_negative_background_status": "UNKNOWN_BACKGROUND",
                "calibration_null_status": "NO_NULL_SUPPORT",
            }
        ]
    )

    row = classify_branch_b_v2_candidates(frame).iloc[0]

    assert row["v2_signal_strength_tier"] == "A_STRONG_Z_GE_10"
    assert row["v2_report_layer_class"] == "internal_review_event"
    assert row["v2_report_visibility"] == "internal_review"
    assert row["v2_filter_action"] != "filter_report_layer_combined_technical_risk"


def test_report_layer_keeps_h6_chr21_truth_sensitive_candidate_visible():
    frame = pd.DataFrame(
        [
            {
                "sample_id": "H6",
                "candidate_id": "H6_chr21_gain",
                "chrom": "chr21",
                "start": 20_700_000,
                "end": 22_300_000,
                "state": "gain",
                "a_abs_zscore": 7.11,
                "same_direction_fraction": 0.92,
                "corrected_amplitude": 0.40,
                "attenuation_ratio": 0.35,
                "clean_bin_fraction": 0.85,
                "high_risk_bin_fraction": 0.05,
                "matched_negative_background_status": "UNKNOWN_BACKGROUND",
                "calibration_null_status": "NO_NULL_SUPPORT",
            }
        ]
    )

    row = classify_branch_b_v2_candidates(frame).iloc[0]

    assert row["v2_report_layer_class"] == "internal_review_event"
    assert row["v2_report_visibility"] == "internal_review"
    assert row["v2_filter_action"] != "filter_report_layer_combined_technical_risk"


def test_report_layer_matched_background_outlier_can_be_report_event():
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
                "clean_bin_fraction": 0.85,
                "high_risk_bin_fraction": 0.05,
                "matched_negative_background_status": "OK",
                "matched_negative_abs_percentile": 0.995,
            }
        ]
    )

    row = classify_branch_b_v2_candidates(frame).iloc[0]

    assert row["v2_report_layer_class"] == "report_event"
    assert row["v2_report_visibility"] == "final_report"


def test_report_layer_routes_sex_chromosome_to_branch_s_event():
    frame = pd.DataFrame(
        [
            {
                "sample_id": "H5",
                "candidate_id": "H5_chrX_loss",
                "chrom": "chrX",
                "start": 1,
                "end": 155_000_000,
                "state": "loss",
                "a_abs_zscore": 62.0,
                "same_direction_fraction": 0.80,
                "corrected_amplitude": -2.1,
                "matched_negative_background_status": "UNKNOWN_BACKGROUND",
            }
        ]
    )

    row = classify_branch_b_v2_candidates(frame).iloc[0]

    assert row["v2_report_layer_class"] == "branch_s_event"
    assert row["v2_report_visibility"] == "branch_s_report_section"
