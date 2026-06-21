from __future__ import annotations

import pandas as pd

from pgta.predict.branch_s import build_branch_s_shadow, summarize_branch_s_shadow


def test_branch_s_shadow_partitions_sex_chromosome_regions_and_never_changes_report():
    bins = pd.DataFrame(
        [
            {
                "chrom": "chrX",
                "start": 0,
                "end": 1_000_000,
                "calibrated_z": 6.0,
                "robust_z": 5.0,
                "is_PAR": False,
            },
            {
                "chrom": "chrX",
                "start": 1_000_000,
                "end": 2_000_000,
                "calibrated_z": 4.0,
                "robust_z": 3.0,
                "par_overlap_fraction": 0.75,
            },
            {
                "chrom": "chrY",
                "start": 0,
                "end": 1_000_000,
                "calibrated_z": -3.0,
                "robust_z": -2.0,
                "is_PAR": False,
            },
        ]
    )
    a_candidates = pd.DataFrame(
        [
            {"chrom": "chrX", "start": 0, "end": 1_000_000, "cnv_type": "gain"},
            {"chrom": "chr2", "start": 0, "end": 1_000_000, "cnv_type": "loss"},
        ]
    )
    gender = pd.DataFrame(
        [
            {
                "sample_id": "S1",
                "sex_call": "XY",
                "predict_gender": "M",
                "sex_call_source": "wisecondorx_bam_consensus",
            }
        ]
    )

    evidence, scores = build_branch_s_shadow(
        sample_id="S1",
        bins=bins,
        a_candidates=a_candidates,
        gender=gender,
    )
    summary = summarize_branch_s_shadow("S1", evidence, scores, gender)

    assert set(evidence["region_class"]) >= {"X_NONPAR", "X_PAR", "Y_NONPAR"}
    assert set(evidence["branch_s_action"]) == {"SHADOW_ONLY"}
    assert set(evidence["final_report_impact"]) == {"none_shadow_only"}
    assert evidence.loc[evidence["region_class"] == "X_NONPAR", "a_candidate_count"].iat[0] == 1

    assert set(scores["sca_state"]) == {"X_GAIN", "X_LOSS", "Y_GAIN", "Y_LOSS"}
    assert scores.loc[scores["sca_state"] == "X_GAIN", "state_score"].iat[0] > 0
    assert scores.loc[scores["sca_state"] == "X_LOSS", "state_score"].iat[0] < 0
    assert set(scores["final_report_impact"]) == {"none_shadow_only"}

    assert summary["sex_call"] == "XY"
    assert summary["replaces_current_sex_calling"] is False
    assert summary["replaces_final_report"] is False


def test_branch_s_shadow_outputs_empty_shadow_tables_when_no_sex_bins():
    evidence, scores = build_branch_s_shadow(
        sample_id="S2",
        bins=pd.DataFrame([{"chrom": "chr1", "start": 0, "end": 1_000_000, "calibrated_z": 2.0}]),
        a_candidates=pd.DataFrame(),
        gender=pd.DataFrame(),
    )
    summary = summarize_branch_s_shadow("S2", evidence, scores, pd.DataFrame())

    assert evidence.empty
    assert set(scores["sca_state"]) == {"X_GAIN", "X_LOSS", "Y_GAIN", "Y_LOSS"}
    assert set(scores["state_score_status"]) == {"INSUFFICIENT_EVIDENCE"}
    assert summary["region_count"] == 0
    assert summary["final_report_impact"] == "none_shadow_only"


def test_branch_s_state_scores_use_branch_a_direction_when_region_z_is_opposite():
    bins = pd.DataFrame(
        [
            {
                "chrom": "chrX",
                "start": 3_000_000,
                "end": 58_000_000,
                "calibrated_z": 2.5,
                "robust_z": 30.0,
                "is_PAR": False,
            },
            {
                "chrom": "chrX",
                "start": 61_000_000,
                "end": 155_000_000,
                "calibrated_z": 2.0,
                "robust_z": 25.0,
                "is_PAR": False,
            },
        ]
    )
    a_candidates = pd.DataFrame(
        [
            {
                "chrom": "chrX",
                "start": 3_000_000,
                "end": 58_000_000,
                "state": "loss",
                "a_zscore": -65.8,
                "a_abs_zscore": 65.8,
            },
            {
                "chrom": "chrX",
                "start": 61_000_000,
                "end": 155_000_000,
                "state": "loss",
                "a_zscore": -77.1,
                "a_abs_zscore": 77.1,
            },
        ]
    )

    _evidence, scores = build_branch_s_shadow(
        sample_id="XO",
        bins=bins,
        a_candidates=a_candidates,
        gender=pd.DataFrame(),
    )

    x_gain = scores.loc[scores["sca_state"] == "X_GAIN", "state_score"].iat[0]
    x_loss = scores.loc[scores["sca_state"] == "X_LOSS", "state_score"].iat[0]
    loss_reason = scores.loc[scores["sca_state"] == "X_LOSS", "state_score_reason"].iat[0]

    assert x_loss < 0
    assert x_gain > 0
    assert loss_reason == "branch_a_only_uncorroborated_by_nonpar_median"


def test_branch_s_does_not_call_strong_x_gain_from_uncorroborated_xy_branch_a_candidate():
    bins = pd.DataFrame(
        [
            {
                "chrom": "chrX",
                "start": 3_000_000,
                "end": 58_000_000,
                "calibrated_z": -0.2,
                "robust_z": -0.1,
                "is_PAR": False,
            },
            {
                "chrom": "chrX",
                "start": 61_000_000,
                "end": 154_500_000,
                "calibrated_z": 0.1,
                "robust_z": 0.2,
                "is_PAR": False,
            },
            {
                "chrom": "chrY",
                "start": 3_000_000,
                "end": 20_000_000,
                "calibrated_z": 0.2,
                "robust_z": 0.1,
                "is_PAR": False,
            },
        ]
    )
    a_candidates = pd.DataFrame(
        [
            {
                "chrom": "chrX",
                "start": 3_000_000,
                "end": 154_500_000,
                "state": "gain",
                "a_zscore": 131.8,
                "a_abs_zscore": 131.8,
            }
        ]
    )
    gender = pd.DataFrame(
        [
            {
                "sample_id": "XY_NORMAL",
                "sex_call": "XY",
                "predict_gender": "M",
                "sex_call_source": "wisecondorx_bam_consensus",
            }
        ]
    )

    evidence, scores = build_branch_s_shadow(
        sample_id="XY_NORMAL",
        bins=bins,
        a_candidates=a_candidates,
        gender=gender,
    )
    summary = summarize_branch_s_shadow("XY_NORMAL", evidence, scores, gender)
    x_gain = scores.loc[scores["sca_state"] == "X_GAIN"].iloc[0]

    assert x_gain["state_score"] < 5.0
    assert x_gain["state_score_reason"] == "branch_a_only_uncorroborated_by_nonpar_median"
    assert summary["branch_a_x_support"] == "present"
    assert summary["sca_candidate_state"] == "none_detected"
    assert summary["sca_confidence_tier"] == "SCA_NO_CALL"
    assert summary["sca_report_layer_class"] == "sca_filtered_or_sex_consistent_event"
    assert "branch_a_only_uncorroborated" in summary["sca_report_layer_reason"]


def test_branch_s_preserves_x_loss_when_branch_a_is_nonpar_corroborated():
    bins = pd.DataFrame(
        [
            {
                "chrom": "chrX",
                "start": 3_000_000,
                "end": 58_000_000,
                "calibrated_z": -4.5,
                "robust_z": -22.0,
                "is_PAR": False,
            },
            {
                "chrom": "chrX",
                "start": 61_000_000,
                "end": 154_500_000,
                "calibrated_z": -4.0,
                "robust_z": -18.0,
                "is_PAR": False,
            },
        ]
    )
    a_candidates = pd.DataFrame(
        [
            {
                "chrom": "chrX",
                "start": 3_000_000,
                "end": 154_500_000,
                "state": "loss",
                "a_zscore": -77.1,
                "a_abs_zscore": 77.1,
            }
        ]
    )
    gender = pd.DataFrame(
        [
            {
                "sample_id": "XO_REVIEW",
                "sex_call": "XX",
                "predict_gender": "F",
                "sex_call_source": "bam_depth_override",
            }
        ]
    )

    evidence, scores = build_branch_s_shadow(
        sample_id="XO_REVIEW",
        bins=bins,
        a_candidates=a_candidates,
        gender=gender,
    )
    summary = summarize_branch_s_shadow("XO_REVIEW", evidence, scores, gender)
    x_loss = scores.loc[scores["sca_state"] == "X_LOSS"].iloc[0]

    assert x_loss["state_score"] >= 30.0
    assert x_loss["state_score_reason"] == "branch_a_candidate_zscore_nonpar_corroborated"
    assert summary["sca_candidate_state"] == "X_LOSS"
    assert summary["sca_confidence_tier"] == "SCA_REVIEW_STRONG"
    assert summary["sca_report_layer_class"] == "sca_report_review_event"


def test_branch_s_preserves_xx_x_loss_review_when_branch_a_is_strong_but_bin_median_is_neutral():
    bins = pd.DataFrame(
        [
            {
                "chrom": "chrX",
                "start": 3_000_000,
                "end": 58_000_000,
                "calibrated_z": -0.1,
                "robust_z": -0.1,
                "is_PAR": False,
            },
            {
                "chrom": "chrX",
                "start": 61_000_000,
                "end": 154_500_000,
                "calibrated_z": 0.1,
                "robust_z": 0.1,
                "is_PAR": False,
            },
        ]
    )
    a_candidates = pd.DataFrame(
        [
            {
                "chrom": "chrX",
                "start": 3_000_000,
                "end": 58_000_000,
                "state": "loss",
                "a_zscore": -62.0,
                "a_abs_zscore": 62.0,
            },
            {
                "chrom": "chrX",
                "start": 61_000_000,
                "end": 154_500_000,
                "state": "loss",
                "a_zscore": -58.0,
                "a_abs_zscore": 58.0,
            },
        ]
    )
    gender = pd.DataFrame(
        [
            {
                "sample_id": "XX_XLOSS",
                "sex_call": "XX",
                "predict_gender": "F",
                "sex_call_source": "bam_depth_override",
            }
        ]
    )

    evidence, scores = build_branch_s_shadow(
        sample_id="XX_XLOSS",
        bins=bins,
        a_candidates=a_candidates,
        gender=gender,
    )
    summary = summarize_branch_s_shadow("XX_XLOSS", evidence, scores, gender)
    x_loss = scores.loc[scores["sca_state"] == "X_LOSS"].iloc[0]

    assert x_loss["state_score"] >= 30.0
    assert x_loss["state_score_reason"] == "branch_a_candidate_zscore_sex_call_compatible_uncorroborated_review"
    assert summary["sca_candidate_state"] == "X_LOSS"
    assert summary["sca_report_layer_class"] == "sca_report_review_event"
    assert summary["sca_report_layer_reason"] == "sca_review_strong_with_sex_call_compatible_branch_a_support"


def test_branch_s_summary_contains_p5_report_boundary_contract():
    bins = pd.DataFrame(
        [
            {
                "chrom": "chrX",
                "start": 3_000_000,
                "end": 58_000_000,
                "calibrated_z": 4.0,
                "robust_z": 22.0,
                "is_PAR": False,
            },
            {
                "chrom": "chrX",
                "start": 60_000,
                "end": 2_600_000,
                "calibrated_z": 1.0,
                "robust_z": 1.0,
                "par_overlap_fraction": 1.0,
            },
            {
                "chrom": "chrY",
                "start": 3_000_000,
                "end": 20_000_000,
                "calibrated_z": 0.5,
                "robust_z": 0.5,
                "is_PAR": False,
            },
        ]
    )
    a_candidates = pd.DataFrame(
        [
            {
                "chrom": "chrX",
                "start": 3_000_000,
                "end": 58_000_000,
                "state": "gain",
                "a_zscore": 42.0,
                "a_abs_zscore": 42.0,
            }
        ]
    )
    gender = pd.DataFrame(
        [
            {
                "sample_id": "SCA1",
                "sex_call": "XY",
                "predict_gender": "M",
                "sex_call_source": "wisecondorx_bam_consensus",
            }
        ]
    )

    evidence, scores = build_branch_s_shadow(
        sample_id="SCA1",
        bins=bins,
        a_candidates=a_candidates,
        gender=gender,
    )
    summary = summarize_branch_s_shadow("SCA1", evidence, scores, gender)

    required_keys = {
        "sample",
        "sex_call",
        "expected_x_ploidy",
        "expected_y_ploidy",
        "x_nonpar_direction",
        "x_par_context",
        "y_nonpar_direction",
        "y_par_or_homology_context",
        "branch_a_x_support",
        "branch_a_y_support",
        "sca_candidate_state",
        "sca_confidence_tier",
        "sca_output_mode",
        "sca_uncertainty_reason",
        "report_text_status",
    }
    assert required_keys.issubset(summary)
    assert summary["sample"] == "SCA1"
    assert summary["expected_x_ploidy"] == 1
    assert summary["expected_y_ploidy"] == 1
    assert summary["x_nonpar_direction"] == "gain"
    assert summary["x_par_context"] == "available"
    assert summary["branch_a_x_support"] == "present"
    assert summary["branch_a_y_support"] == "absent"
    assert summary["sca_candidate_state"] == "X_GAIN"
    assert summary["sca_confidence_tier"] == "SCA_REVIEW_STRONG"
    assert summary["sca_output_mode"] == "review_development_only"
    assert summary["report_text_status"] == "development_only_not_final_reportable"
    assert "locked_sca_truth_incomplete" in summary["sca_uncertainty_reason"]
