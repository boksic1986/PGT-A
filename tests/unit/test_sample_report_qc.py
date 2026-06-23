from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from pgta.predict.sample_report_qc import build_sample_report_qc_table, classify_sample_reportability


def test_bam_qc_fail_promotes_to_no_call_recommended():
    row = classify_sample_reportability(
        {
            "sample_id": "S1",
            "bam_qc_status": "BAM_QC_FAIL",
            "cnv_qc_status": "PASS",
            "mad_log1p": 0.14,
            "cohort_size": 5,
        }
    )

    assert row["sample_report_qc_status"] == "NO_CALL_RECOMMENDED"
    assert "BAM_QC_FAIL" in row["sample_report_qc_reasons"]


def test_high_noise_with_input_pass_is_review_not_no_call_for_small_cohort():
    row = classify_sample_reportability(
        {
            "sample_id": "S60",
            "bam_qc_status": "BAM_QC_PASS",
            "cnv_qc_status": "PASS",
            "mad_log1p": 0.345,
            "cohort_mad_log1p_median": 0.14,
            "cohort_mad_log1p_robust_z": 12.0,
            "cohort_size": 5,
            "autosomal_cn_outside_1p7_2p3_fraction": 0.22,
            "report_event_count": 18,
            "multi_chromosome_report_event_count": 7,
        }
    )

    assert row["sample_report_qc_status"] == "SAMPLE_QUALITY_REVIEW"
    assert "HIGH_MAD_LOG1P_RELATIVE_OUTLIER" in row["sample_report_qc_reasons"]
    assert row["rebuild_or_resample_recommendation"] == "review_library_or_resequence_before_release"


def test_single_coherent_truth_like_cnv_does_not_trigger_review_by_itself():
    row = classify_sample_reportability(
        {
            "sample_id": "Truth1",
            "bam_qc_status": "BAM_QC_PASS",
            "cnv_qc_status": "PASS",
            "mad_log1p": 0.14,
            "cohort_size": 6,
            "report_event_count": 1,
            "multi_chromosome_report_event_count": 1,
            "autosomal_cn_outside_1p7_2p3_fraction": 0.04,
        }
    )

    assert row["sample_report_qc_status"] == "PASS_REPORTABLE"
    assert row["sample_report_qc_reasons"] == "PASS"


def test_chr_y_artifact_context_does_not_fail_sample():
    row = classify_sample_reportability(
        {
            "sample_id": "XY1",
            "bam_qc_status": "BAM_QC_PASS",
            "cnv_qc_status": "PASS",
            "mad_log1p": 0.13,
            "cohort_size": 8,
            "chrY_cn_artifact_context": "sex_chrom_ref_ratio_not_interpretable",
            "report_event_count": 0,
        }
    )

    assert row["sample_report_qc_status"] == "PASS_REPORTABLE"
    assert "chrY" in row["possible_contamination_or_mixture_context"]


def test_build_sample_report_qc_table_marks_60_like_relative_outlier_as_review():
    bam_qc = pd.DataFrame(
        [
            {"sample_id": "56", "bam_qc_status": "BAM_QC_PASS"},
            {"sample_id": "59", "bam_qc_status": "BAM_QC_PASS"},
            {"sample_id": "60", "bam_qc_status": "BAM_QC_PASS"},
            {"sample_id": "61", "bam_qc_status": "BAM_QC_PASS"},
            {"sample_id": "62", "bam_qc_status": "BAM_QC_PASS"},
        ]
    )
    cnv_qc = pd.DataFrame(
        [
            {"sample_id": "56", "status": "PASS", "mad_log1p": 0.13},
            {"sample_id": "59", "status": "PASS", "mad_log1p": 0.15},
            {"sample_id": "60", "status": "PASS", "mad_log1p": 0.345},
            {"sample_id": "61", "status": "PASS", "mad_log1p": 0.12},
            {"sample_id": "62", "status": "PASS", "mad_log1p": 0.14},
        ]
    )
    v2_sample = pd.DataFrame(
        [
            {"sample_id": "56", "v2_report_event_count": 2},
            {"sample_id": "59", "v2_report_event_count": 4},
            {"sample_id": "60", "v2_report_event_count": 18},
            {"sample_id": "61", "v2_report_event_count": 1},
            {"sample_id": "62", "v2_report_event_count": 1},
        ]
    )

    table = build_sample_report_qc_table(bam_qc=bam_qc, cnv_qc=cnv_qc, v2_sample_summary=v2_sample)
    status_by_sample = dict(zip(table["sample_id"], table["sample_report_qc_status"]))

    assert status_by_sample["60"] == "SAMPLE_QUALITY_REVIEW"
    assert status_by_sample["56"] == "PASS_REPORTABLE"
    assert status_by_sample["61"] == "PASS_REPORTABLE"
