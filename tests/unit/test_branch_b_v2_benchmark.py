from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from pgta.predict.branch_b.v2_benchmark import build_v2_benchmark


def test_v2_benchmark_preserves_truth_overlap_without_legacy_fields(tmp_path: Path):
    classification = tmp_path / "H6.candidate_classification.tsv"
    truth = tmp_path / "truth.tsv"
    truth_metrics = tmp_path / "truth_metrics.tsv"
    sample_summary = tmp_path / "sample_summary.tsv"
    summary = tmp_path / "summary.json"

    pd.DataFrame(
        [
            {
                "sample_id": "H6",
                "candidate_id": "H6.A0001",
                "chrom": "chr21",
                "start": 20_000_000,
                "end": 23_000_000,
                "state": "gain",
                "a_abs_zscore": 7.11,
                "final_disposition": "LIKELY_ARTIFACT",
                "branch_b_keep_event": 0,
                "v2_candidate_class": "V2_POSITIVE_SUPPORT_REVIEW",
                "v2_classifier_action": "V2_REVIEW_POSITIVE_SUPPORT",
                "v2_evidence_tier": "UNKNOWN_BACKGROUND_POSITIVE_SUPPORT",
                "v2_evidence_gate": "NO_HARD_SUPPRESSION",
                "v2_disposition": "background_unknown_review",
                "v2_length_tier": "reportable_candidate_ge2mb",
                "v2_clean_support_label": "CLEAN_SUPPORT_AVAILABLE",
                "attenuation_ratio": 0.42,
            }
        ]
    ).to_csv(classification, sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "sample_id": "H6",
                "chrom": "chr21",
                "start": 20_700_000,
                "end": 22_300_000,
                "expected_state": "gain",
            }
        ]
    ).to_csv(truth, sep="\t", index=False)

    _, _, payload = build_v2_benchmark(
        classification_paths=[str(classification)],
        truth_tsv=str(truth),
        sample_ids=["H6"],
        reference_id="h_r0_shadow_ref_20260619",
        output_truth_metrics=str(truth_metrics),
        output_sample_summary=str(sample_summary),
        output_summary=str(summary),
    )

    assert payload["truth_event_count"] == 1
    assert payload["truth_preserved_count"] == 1
    assert payload["FN_count"] == 0
    assert payload["truth_hard_suppressed_count"] == 0
    assert payload["legacy_branch_b_decision_fields_used"] is False

    truth_row = pd.read_csv(truth_metrics, sep="\t").iloc[0]
    assert truth_row["top_v2_disposition"] == "background_unknown_review"
    assert truth_row["top_v2_length_tier"] == "reportable_candidate_ge2mb"
    assert truth_row["top_v2_clean_support_label"] == "CLEAN_SUPPORT_AVAILABLE"
    assert truth_row["top_attenuation_ratio"] == 0.42

    written = json.loads(summary.read_text(encoding="utf-8"))
    assert written["ignored_legacy_decision_fields"] == [
        "final_disposition",
        "branch_b_keep_event",
        "branch_b_report_class",
        "branch_b_artifact_status",
    ]


def test_v2_benchmark_counts_hard_suppression_as_fn(tmp_path: Path):
    classification = tmp_path / "Y1.candidate_classification.tsv"
    truth = tmp_path / "truth.tsv"

    pd.DataFrame(
        [
            {
                "sample_id": "Y1",
                "candidate_id": "Y1.A0001",
                "chrom": "chr16",
                "start": 500_000,
                "end": 5_000_000,
                "state": "loss",
                "a_abs_zscore": 29.0,
                "v2_candidate_class": "V2_LEGACY_ARTIFACT_SHOULD_NOT_EXIST",
                "v2_classifier_action": "V2_SUPPRESS_ARTIFACT",
                "v2_evidence_tier": "UNKNOWN_BACKGROUND_POSITIVE_SUPPORT",
            }
        ]
    ).to_csv(classification, sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "sample_id": "Y1",
                "chrom": "chr16",
                "start": 1_000_000,
                "end": 4_000_000,
                "expected_state": "loss",
            }
        ]
    ).to_csv(truth, sep="\t", index=False)

    _, _, payload = build_v2_benchmark(
        classification_paths=[str(classification)],
        truth_tsv=str(truth),
        sample_ids=["Y1"],
        reference_id="h_r0_shadow_ref_20260619",
        output_truth_metrics=str(tmp_path / "truth_metrics.tsv"),
        output_sample_summary=str(tmp_path / "sample_summary.tsv"),
        output_summary=str(tmp_path / "summary.json"),
    )

    assert payload["truth_event_count"] == 1
    assert payload["truth_preserved_count"] == 0
    assert payload["FN_count"] == 1
    assert payload["truth_hard_suppressed_count"] == 1


def test_no_hard_suppression_review_preserves_truth_overlap(tmp_path: Path):
    classification = tmp_path / "H5.candidate_classification.tsv"
    truth = tmp_path / "truth.tsv"

    pd.DataFrame(
        [
            {
                "sample_id": "H5",
                "candidate_id": "H5.A0001",
                "chrom": "chrX",
                "start": 1,
                "end": 155_000_000,
                "state": "loss",
                "a_abs_zscore": 22.0,
                "v2_candidate_class": "V2_NO_CALL_CONTRACT_RISK",
                "v2_classifier_action": "V2_REVIEW_NO_HARD_SUPPRESSION",
                "v2_evidence_tier": "UNKNOWN_BACKGROUND_CONTRACT_RISK",
            }
        ]
    ).to_csv(classification, sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "sample_id": "H5",
                "chrom": "chrX",
                "start": 178_624,
                "end": 64_743_458,
                "expected_state": "loss",
            }
        ]
    ).to_csv(truth, sep="\t", index=False)

    _, _, payload = build_v2_benchmark(
        classification_paths=[str(classification)],
        truth_tsv=str(truth),
        sample_ids=["H5"],
        reference_id="h_r0_shadow_ref_20260619",
        output_truth_metrics=str(tmp_path / "truth_metrics.tsv"),
        output_sample_summary=str(tmp_path / "sample_summary.tsv"),
        output_summary=str(tmp_path / "summary.json"),
    )

    assert payload["truth_event_count"] == 1
    assert payload["truth_preserved_count"] == 1
    assert payload["FN_count"] == 0
    assert payload["truth_hard_suppressed_count"] == 0
