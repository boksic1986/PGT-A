from __future__ import annotations

import pandas as pd

from pgta.predict.branch_s import (
    build_branch_s_shadow,
    prepare_candidates,
    summarize_branch_s_shadow,
    summarize_sex_chrom_lowres_context,
)


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


def test_branch_s_uses_segment_level_nonpar_support_when_global_median_is_neutral():
    bins = pd.DataFrame(
        [
            {
                "chrom": "chrX",
                "start": idx * 1_000_000,
                "end": (idx + 1) * 1_000_000,
                "calibrated_z": 4.0 if idx in {2, 3} else 0.0,
                "robust_z": 8.0 if idx in {2, 3} else 0.0,
                "is_PAR": False,
            }
            for idx in range(10)
        ]
    )
    a_candidates = pd.DataFrame(
        [
            {
                "chrom": "chrX",
                "start": 2_000_000,
                "end": 4_000_000,
                "state": "gain",
                "a_zscore": 12.5,
                "a_abs_zscore": 12.5,
            }
        ]
    )
    gender = pd.DataFrame(
        [
            {
                "sample_id": "XY_SEGMENT",
                "sex_call": "XY",
                "predict_gender": "M",
                "sex_call_source": "wisecondorx_bam_consensus",
            }
        ]
    )

    evidence, scores = build_branch_s_shadow(
        sample_id="XY_SEGMENT",
        bins=bins,
        a_candidates=a_candidates,
        gender=gender,
    )
    summary = summarize_branch_s_shadow("XY_SEGMENT", evidence, scores, gender)
    x_gain = scores.loc[scores["sca_state"] == "X_GAIN"].iloc[0]

    assert x_gain["state_score"] == 12.5
    assert x_gain["state_score_reason"] == "branch_a_candidate_zscore_segment_nonpar_corroborated"
    assert summary["x_nonpar_direction"] == "gain"
    assert summary["sca_candidate_state"] == "X_GAIN"
    assert summary["sca_report_layer_class"] == "sca_internal_review_event"


def test_branch_s_segment_support_ignores_sparse_mean_skew_for_xy_x_gain():
    bins = pd.DataFrame(
        [
            {
                "chrom": "chrX",
                "start": idx * 1_000_000,
                "end": (idx + 1) * 1_000_000,
                "calibrated_z": 25.0 if idx in {0, 1} else 0.0,
                "robust_z": 45.0 if idx in {0, 1} else 0.0,
                "is_PAR": False,
            }
            for idx in range(100)
        ]
    )
    a_candidates = pd.DataFrame(
        [
            {
                "chrom": "chrX",
                "start": 0,
                "end": 100_000_000,
                "state": "gain",
                "a_zscore": 120.0,
                "a_abs_zscore": 120.0,
            }
        ]
    )
    gender = pd.DataFrame(
        [
            {
                "sample_id": "XY_MEAN_SKEW",
                "sex_call": "XY",
                "predict_gender": "M",
                "sex_call_source": "wisecondorx_bam_consensus",
            }
        ]
    )

    evidence, scores = build_branch_s_shadow(
        sample_id="XY_MEAN_SKEW",
        bins=bins,
        a_candidates=a_candidates,
        gender=gender,
    )
    summary = summarize_branch_s_shadow("XY_MEAN_SKEW", evidence, scores, gender)
    x_gain = scores.loc[scores["sca_state"] == "X_GAIN"].iloc[0]

    assert x_gain["state_score"] == 0.0
    assert x_gain["state_score_reason"] == "branch_a_only_uncorroborated_by_nonpar_median"
    assert summary["sca_candidate_state"] == "none_detected"
    assert summary["sca_report_layer_class"] == "sca_filtered_or_sex_consistent_event"


def test_branch_s_promotes_branch_a_supported_weak_x_loss_to_report_review():
    bins = pd.DataFrame(
        [
            {
                "chrom": "chrX",
                "start": 3_000_000,
                "end": 58_000_000,
                "calibrated_z": -0.05,
                "robust_z": -0.05,
                "is_PAR": False,
            },
            {
                "chrom": "chrX",
                "start": 61_000_000,
                "end": 154_000_000,
                "calibrated_z": 0.05,
                "robust_z": 0.05,
                "is_PAR": False,
            },
        ]
    )
    a_candidates = pd.DataFrame(
        [
            {
                "chrom": "chrX",
                "start": 3_000_000,
                "end": 154_000_000,
                "state": "loss",
                "a_zscore": -14.5,
                "a_abs_zscore": 14.5,
            }
        ]
    )
    gender = pd.DataFrame(
        [
            {
                "sample_id": "XX_MOSAIC_XLOSS",
                "sex_call": "XX",
                "predict_gender": "F",
                "sex_call_source": "wisecondorx_bam_consensus",
            }
        ]
    )

    evidence, scores = build_branch_s_shadow(
        sample_id="XX_MOSAIC_XLOSS",
        bins=bins,
        a_candidates=a_candidates,
        gender=gender,
    )
    summary = summarize_branch_s_shadow("XX_MOSAIC_XLOSS", evidence, scores, gender)

    assert summary["sca_candidate_state"] == "X_LOSS"
    assert summary["sca_confidence_tier"] == "SCA_REVIEW_WEAK"
    assert summary["sca_report_layer_class"] == "sca_report_review_event"
    assert summary["sca_report_layer_reason"] == "sca_review_weak_report_visible_with_branch_a_support"


def test_branch_s_par_only_signal_is_context_not_sca_call():
    bins = pd.DataFrame(
        [
            {
                "chrom": "chrX",
                "start": 0,
                "end": 1_000_000,
                "calibrated_z": 9.0,
                "robust_z": 12.0,
                "par_overlap_fraction": 1.0,
            },
            {
                "chrom": "chrX",
                "start": 3_000_000,
                "end": 4_000_000,
                "calibrated_z": 0.0,
                "robust_z": 0.0,
                "is_PAR": False,
            },
        ]
    )
    gender = pd.DataFrame(
        [
            {
                "sample_id": "PAR_ONLY",
                "sex_call": "XY",
                "predict_gender": "M",
                "sex_call_source": "wisecondorx_bam_consensus",
            }
        ]
    )

    evidence, scores = build_branch_s_shadow(
        sample_id="PAR_ONLY",
        bins=bins,
        a_candidates=pd.DataFrame(),
        gender=gender,
    )
    summary = summarize_branch_s_shadow("PAR_ONLY", evidence, scores, gender)

    assert summary["x_par_context"] == "available"
    assert summary["sca_candidate_state"] == "none_detected"
    assert summary["sca_report_layer_class"] == "sca_no_call"


def test_branch_s_lowres_absence_is_context_not_suppression_for_short_segment():
    bins = pd.DataFrame(
        [
            {
                "chrom": "chrX",
                "start": idx * 1_000_000,
                "end": (idx + 1) * 1_000_000,
                "calibrated_z": -3.0 if idx in {5, 6} else 0.0,
                "robust_z": -8.0 if idx in {5, 6} else 0.0,
                "is_PAR": False,
            }
            for idx in range(10)
        ]
    )
    a_candidates = pd.DataFrame(
        [
            {
                "chrom": "chrX",
                "start": 5_000_000,
                "end": 7_000_000,
                "state": "loss",
                "a_zscore": -11.0,
                "a_abs_zscore": 11.0,
            }
        ]
    )
    gender = pd.DataFrame(
        [
            {
                "sample_id": "XY_SHORT_XLOSS",
                "sex_call": "XY",
                "predict_gender": "M",
                "sex_call_source": "wisecondorx_bam_consensus",
            }
        ]
    )
    lowres_context = summarize_sex_chrom_lowres_context(
        prepare_candidates(a_candidates),
        lowres_2mb_events=pd.DataFrame(),
        lowres_3mb_events=pd.DataFrame(),
    )

    evidence, scores = build_branch_s_shadow(
        sample_id="XY_SHORT_XLOSS",
        bins=bins,
        a_candidates=a_candidates,
        gender=gender,
    )
    summary = summarize_branch_s_shadow(
        "XY_SHORT_XLOSS",
        evidence,
        scores,
        gender,
        lowres_context=lowres_context,
    )

    assert summary["sex_chrom_lowres_2mb_context"] == "not_informative_short_or_boundary_event"
    assert summary["sex_chrom_lowres_3mb_context"] == "not_informative_short_or_boundary_event"
    assert summary["sca_candidate_state"] == "X_LOSS"
    assert summary["sca_report_layer_class"] == "sca_internal_review_event"


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
    assert summary["sex_chrom_lowres_2mb_context"] == "not_configured"
    assert summary["sex_chrom_lowres_final_impact"] == "context_only_not_filter"
    assert "locked_sca_truth_incomplete" in summary["sca_uncertainty_reason"]
