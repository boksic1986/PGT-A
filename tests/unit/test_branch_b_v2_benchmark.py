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
                "v2_filter_action": "keep_background_unknown_review",
                "v2_length_tier": "reportable_candidate_ge2mb",
                "v2_clean_support_label": "CLEAN_SUPPORT_AVAILABLE",
                "v2_burden_reduction_tier": "background_unknown_review",
                "v2_burden_reduction_action": "stratify_background_unknown_review",
                "v2_burden_reduction_reason": "background_unknown_review_preserved_for_truth_safety",
                "v2_burden_evidence_tags": "[CNVpro-inspired] length_tier=reportable_candidate_ge2mb;[CNVpro-like] gc_rc_context=GC_RC_ATTENUATED_SEVERE",
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
    assert truth_row["top_v2_filter_action"] == "keep_background_unknown_review"
    assert truth_row["top_v2_length_tier"] == "reportable_candidate_ge2mb"
    assert truth_row["top_v2_clean_support_label"] == "CLEAN_SUPPORT_AVAILABLE"
    assert truth_row["top_v2_burden_reduction_tier"] == "background_unknown_review"
    assert truth_row["top_v2_burden_reduction_action"] == "stratify_background_unknown_review"
    assert truth_row["top_attenuation_ratio"] == 0.42

    written = json.loads(summary.read_text(encoding="utf-8"))
    assert written["ignored_legacy_decision_fields"] == [
        "final_disposition",
        "branch_b_keep_event",
        "branch_b_report_class",
        "branch_b_artifact_status",
    ]
    assert written["burden_stratification_counts"]["v2_burden_reduction_tier"] == {
        "background_unknown_review": 1,
    }


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


def test_v2_filter_contract_suppression_counts_as_fn(tmp_path: Path):
    classification = tmp_path / "Y1.candidate_classification.tsv"
    truth = tmp_path / "truth.tsv"

    pd.DataFrame(
        [
            {
                "sample_id": "Y1",
                "candidate_id": "Y1.A0001",
                "chrom": "chr1",
                "start": 1_000_000,
                "end": 12_000_000,
                "state": "gain",
                "a_abs_zscore": 32.0,
                "v2_candidate_class": "V2_NO_CALL_CONTRACT_RISK",
                "v2_classifier_action": "V2_REVIEW_NO_HARD_SUPPRESSION",
                "v2_filter_action": "suppress_workflow_contract_risk",
                "v2_evidence_tier": "UNKNOWN_BACKGROUND_REF_CONTRACT_RISK",
            }
        ]
    ).to_csv(classification, sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "sample_id": "Y1",
                "chrom": "chr1",
                "start": 2_000_000,
                "end": 10_000_000,
                "expected_state": "gain",
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


def test_v2_benchmark_outputs_sample_burden_reduction_counts(tmp_path: Path):
    classification = tmp_path / "S1.candidate_classification.tsv"
    sample_summary = tmp_path / "sample_summary.tsv"
    summary = tmp_path / "summary.json"

    pd.DataFrame(
        [
            {
                "sample_id": "S1",
                "candidate_id": "S1.A0001",
                "chrom": "chr21",
                "start": 20_000_000,
                "end": 23_000_000,
                "state": "gain",
                "a_abs_zscore": 7.11,
                "v2_candidate_class": "V2_POSITIVE_SUPPORT_REVIEW",
                "v2_classifier_action": "V2_REVIEW_POSITIVE_SUPPORT",
                "v2_evidence_tier": "UNKNOWN_BACKGROUND_POSITIVE_SUPPORT",
                "v2_filter_action": "keep_background_unknown_review",
                "v2_disposition": "background_unknown_review",
                "v2_length_tier": "reportable_candidate_ge2mb",
                "v2_clean_support_label": "CLEAN_SUPPORT_AVAILABLE",
                "v2_b_signal_context_label": "B_SIGNAL_SUPPORTED_A_DIRECTION",
                "v2_gc_rc_context_label": "GC_RC_STABLE",
                "v2_burden_reduction_tier": "background_unknown_review",
                "v2_burden_reduction_action": "stratify_background_unknown_review",
            },
            {
                "sample_id": "S1",
                "candidate_id": "S1.A0002",
                "chrom": "chrX",
                "start": 1,
                "end": 155_000_000,
                "state": "loss",
                "v2_candidate_class": "V2_SEX_CHROMOSOME_REVIEW",
                "v2_filter_action": "route_to_branch_s_review",
                "v2_disposition": "sca_branch_s_review",
                "v2_burden_reduction_tier": "branch_s_review",
                "v2_burden_reduction_action": "route_to_branch_s_review",
            },
            {
                "sample_id": "S1",
                "candidate_id": "S1.A0003",
                "chrom": "chr16",
                "start": 500_000,
                "end": 4_500_000,
                "state": "loss",
                "v2_candidate_class": "V2_POSITIVE_SUPPORT_REVIEW",
                "v2_filter_action": "downgrade_to_technical_risk_review",
                "v2_disposition": "technical_risk_review",
                "v2_burden_reduction_tier": "technical_risk_review",
                "v2_burden_reduction_action": "downgrade_to_technical_risk_review",
            },
        ]
    ).to_csv(classification, sep="\t", index=False)

    _, written_sample_summary, payload = build_v2_benchmark(
        classification_paths=[str(classification)],
        truth_tsv="",
        sample_ids=["S1"],
        reference_id="h_r0_shadow_ref_20260619",
        output_truth_metrics=str(tmp_path / "truth_metrics.tsv"),
        output_sample_summary=str(sample_summary),
        output_summary=str(summary),
    )

    row = written_sample_summary.iloc[0]
    assert row["v2_report_candidate_burden_count"] == 0
    assert row["v2_review_candidate_burden_count"] == 1
    assert row["v2_background_unknown_review_burden_count"] == 1
    assert row["v2_technical_risk_burden_count"] == 1
    assert row["v2_branch_s_review_burden_count"] == 1
    assert payload["burden_stratification_counts"]["v2_burden_reduction_action"] == {
        "downgrade_to_technical_risk_review": 1,
        "route_to_branch_s_review": 1,
        "stratify_background_unknown_review": 1,
    }


def test_v2_benchmark_report_layer_counts_and_filtered_event_ledger(tmp_path: Path):
    classification = tmp_path / "S1.candidate_classification.tsv"
    sample_summary = tmp_path / "sample_summary.tsv"
    summary = tmp_path / "summary.json"
    filtered_events = tmp_path / "filtered_events.tsv"
    filtered_events_json = tmp_path / "filtered_events.json"
    report_events = tmp_path / "report_events.tsv"
    report_events_json = tmp_path / "report_events.json"

    pd.DataFrame(
        [
            {
                "sample_id": "S1",
                "candidate_id": "S1.report",
                "chrom": "chr21",
                "start": 20_000_000,
                "end": 42_000_000,
                "state": "gain",
                "v2_report_layer_class": "report_event",
                "v2_report_visibility": "final_report",
                "v2_burden_reduction_tier": "report_event",
                "v2_filter_action": "keep_report_event",
            },
            {
                "sample_id": "S1",
                "candidate_id": "S1.review",
                "chrom": "chr6",
                "start": 1_000_000,
                "end": 5_000_000,
                "state": "gain",
                "v2_report_layer_class": "internal_review_event",
                "v2_report_visibility": "internal_review",
                "v2_burden_reduction_tier": "internal_review_event",
                "v2_filter_action": "keep_internal_review_event",
            },
            {
                "sample_id": "S1",
                "candidate_id": "S1.filtered",
                "chrom": "chr16",
                "start": 500_000,
                "end": 1_800_000,
                "state": "loss",
                "v2_report_layer_class": "filtered_event",
                "v2_report_visibility": "audit_only",
                "v2_burden_reduction_tier": "filtered_event",
                "v2_filter_action": "filter_report_layer_combined_technical_risk",
                "v2_report_filter_rule_tags": "short_or_focal;low_clean_high_risk;b_signal_not_supportive;gc_rc_attenuated",
            },
            {
                "sample_id": "S1",
                "candidate_id": "S1.branch_s",
                "chrom": "chrX",
                "start": 1,
                "end": 155_000_000,
                "state": "loss",
                "v2_report_layer_class": "branch_s_event",
                "v2_report_visibility": "branch_s_report_section",
                "v2_burden_reduction_tier": "branch_s_event",
                "v2_filter_action": "route_to_branch_s_review",
            },
        ]
    ).to_csv(classification, sep="\t", index=False)

    _, written_sample_summary, payload = build_v2_benchmark(
        classification_paths=[str(classification)],
        truth_tsv="",
        sample_ids=["S1"],
        reference_id="h_r0_shadow_ref_20260619",
        output_truth_metrics=str(tmp_path / "truth_metrics.tsv"),
        output_sample_summary=str(sample_summary),
        output_summary=str(summary),
        output_filtered_events=str(filtered_events),
        output_filtered_events_json=str(filtered_events_json),
        output_report_events=str(report_events),
        output_report_events_json=str(report_events_json),
    )

    row = written_sample_summary.iloc[0]
    assert row["v2_report_event_count"] == 1
    assert row["v2_internal_review_event_count"] == 1
    assert row["v2_filtered_event_count"] == 1
    assert row["v2_branch_s_event_count"] == 1
    assert payload["v2_report_event_count"] == 1
    assert payload["v2_internal_review_event_count"] == 1
    assert payload["v2_filtered_event_count"] == 1
    assert payload["v2_branch_s_event_count"] == 1

    filtered_df = pd.read_csv(filtered_events, sep="\t")
    assert filtered_df["candidate_id"].tolist() == ["S1.filtered"]
    filtered_payload = json.loads(filtered_events_json.read_text(encoding="utf-8"))
    assert filtered_payload["filtered_event_count"] == 1
    assert filtered_payload["filtered_event_rule_counts"] == {
        "b_signal_not_supportive": 1,
        "gc_rc_attenuated": 1,
        "low_clean_high_risk": 1,
        "short_or_focal": 1,
    }

    report_df = pd.read_csv(report_events, sep="\t")
    assert report_df["candidate_id"].tolist() == ["S1.report"]
    report_payload = json.loads(report_events_json.read_text(encoding="utf-8"))
    assert report_payload["report_event_count"] == 1
    assert report_payload["report_event_ids"] == ["S1.report"]


def test_v2_benchmark_marks_sample_report_burden_without_candidate_demotion(tmp_path: Path):
    classification = tmp_path / "S1.candidate_classification.tsv"
    sample_summary = tmp_path / "sample_summary.tsv"
    summary = tmp_path / "summary.json"

    pd.DataFrame(
        [
            {
                "sample_id": "S1",
                "candidate_id": f"S1.report{i}",
                "chrom": f"chr{i}",
                "start": i * 1_000_000,
                "end": i * 1_000_000 + 20_000_000,
                "state": "gain",
                "v2_report_layer_class": "report_event",
                "v2_report_visibility": "final_report",
                "v2_filter_action": "keep_report_event",
            }
            for i in range(1, 4)
        ]
    ).to_csv(classification, sep="\t", index=False)

    _, written_sample_summary, payload = build_v2_benchmark(
        classification_paths=[str(classification)],
        truth_tsv="",
        sample_ids=["S1"],
        reference_id="h_r0_shadow_ref_20260619",
        output_truth_metrics=str(tmp_path / "truth_metrics.tsv"),
        output_sample_summary=str(sample_summary),
        output_summary=str(summary),
    )

    row = written_sample_summary.iloc[0]
    assert row["v2_report_event_count"] == 3
    assert row["sample_report_burden_flag"] == 1
    assert row["sample_report_burden_reason"] == "report_event_count_ge_3"
    assert payload["sample_report_burden_flag_count"] == 1
    assert payload["sample_report_burden_threshold"] == 3


def test_v2_report_layer_filtered_truth_overlap_counts_as_fn(tmp_path: Path):
    classification = tmp_path / "H6.candidate_classification.tsv"
    truth = tmp_path / "truth.tsv"

    pd.DataFrame(
        [
            {
                "sample_id": "H6",
                "candidate_id": "H6.filtered",
                "chrom": "chr21",
                "start": 20_000_000,
                "end": 23_000_000,
                "state": "gain",
                "a_abs_zscore": 7.11,
                "v2_candidate_class": "V2_POSITIVE_SUPPORT_REVIEW",
                "v2_report_layer_class": "filtered_event",
                "v2_report_visibility": "audit_only",
                "v2_filter_action": "filter_report_layer_combined_technical_risk",
                "v2_burden_reduction_tier": "filtered_event",
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
        output_truth_metrics=str(tmp_path / "truth_metrics.tsv"),
        output_sample_summary=str(tmp_path / "sample_summary.tsv"),
        output_summary=str(tmp_path / "summary.json"),
    )

    assert payload["truth_event_count"] == 1
    assert payload["truth_preserved_count"] == 0
    assert payload["FN_count"] == 1
    assert payload["truth_report_layer_filtered_count"] == 1
