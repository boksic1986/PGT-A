from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from pgta.predict.bam_qc import classify_predict_bam_qc, summarize_idxstats_autosomes


def test_predict_bam_qc_passes_stable_library_metrics():
    row = classify_predict_bam_qc(
        {
            "sample_id": "S1",
            "clean_reads": 12_000_000,
            "clean_gc_content": 0.42,
            "duplication_rate": 0.35,
            "insert_peak": 180,
            "mapping_rate_bam": 0.992,
            "proper_pair_rate_bam": 0.975,
            "autosome_bias_cv": 0.08,
        }
    )

    assert row["bam_qc_status"] == "BAM_QC_PASS"
    assert row["bam_qc_reasons"] == "PASS"


def test_predict_bam_qc_reviews_missing_fastp_but_valid_bam():
    row = classify_predict_bam_qc(
        {
            "sample_id": "S2",
            "mapping_rate_bam": 0.991,
            "proper_pair_rate_bam": 0.972,
            "autosome_bias_cv": 0.08,
        }
    )

    assert row["bam_qc_status"] == "BAM_QC_REVIEW"
    assert "FASTP_METRICS_MISSING" in row["bam_qc_reasons"]


def test_predict_bam_qc_fails_severe_mapping_problem():
    row = classify_predict_bam_qc(
        {
            "sample_id": "S3",
            "clean_reads": 12_000_000,
            "clean_gc_content": 0.42,
            "duplication_rate": 0.35,
            "insert_peak": 180,
            "mapping_rate_bam": 0.72,
            "proper_pair_rate_bam": 0.70,
            "autosome_bias_cv": 0.08,
        }
    )

    assert row["bam_qc_status"] == "BAM_QC_FAIL"
    assert "LOW_MAPPING_RATE_FAIL" in row["bam_qc_reasons"]


def test_idxstats_autosome_bias_summary_ignores_chr_y_artifacts():
    rows = [
        {"chrom": "chr1", "length": 100, "mapped_reads": 1000},
        {"chrom": "chr2", "length": 100, "mapped_reads": 1100},
        {"chrom": "chrX", "length": 100, "mapped_reads": 700},
        {"chrom": "chrY", "length": 100, "mapped_reads": 0},
    ]

    summary = summarize_idxstats_autosomes(rows)

    assert summary["autosome_mapped_reads"] == 2100
    assert 0.0 <= summary["autosome_bias_cv"] < 0.1
    assert summary["chrY_mapped_reads"] == 0
