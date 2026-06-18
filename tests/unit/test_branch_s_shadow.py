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
