import json

import pandas as pd

from pgta.predict.branch_b.plot_manifest_audit import _read_table, build_report_table_ablation_audit


def test_read_table_treats_whitespace_only_file_as_empty_table(tmp_path):
    support_tsv = tmp_path / "empty_support.tsv"
    support_tsv.write_text("\n\n", encoding="utf-8")

    frame = _read_table(support_tsv, columns=["candidate_id", "plot_support_class"])

    assert frame.empty
    assert frame.columns.tolist() == ["candidate_id", "plot_support_class"]


def test_report_table_ablation_downgrades_only_z_supported_cn_not_supported_non_truth(tmp_path):
    report_events = pd.DataFrame(
        [
            {
                "sample_id": "S1",
                "candidate_id": "S1.lowconf",
                "chrom": "chr4",
                "start": 10_000_000,
                "end": 20_000_000,
                "state": "gain",
                "v2_report_layer_class": "report_event",
                "v2_report_visibility": "report_strong_event",
            },
            {
                "sample_id": "S1",
                "candidate_id": "S1.truth_weak",
                "chrom": "chr21",
                "start": 20_700_000,
                "end": 22_300_000,
                "state": "gain",
                "v2_report_layer_class": "report_event",
                "v2_report_visibility": "report_weak_event",
            },
            {
                "sample_id": "S1",
                "candidate_id": "S1.cn_weak",
                "chrom": "chr7",
                "start": 50_000_000,
                "end": 55_000_000,
                "state": "loss",
                "v2_report_layer_class": "report_event",
                "v2_report_visibility": "report_weak_event",
            },
        ]
    )
    manifest = pd.DataFrame(
        [
            {
                "sample_id": "S1",
                "event_id": "evt_low",
                "candidate_id": "S1.lowconf",
                "chrom": "chr4",
                "start": 10_000_000,
                "end": 20_000_000,
                "state": "gain",
                "event_layer": "autosomal_report",
                "plot_visibility": "review_plot_only",
                "plot_layer_class": "internal_review_event_candidate",
                "plot_support_class": "Z_SUPPORTED_CN_NOT_SUPPORTED",
            },
            {
                "sample_id": "S1",
                "event_id": "evt_truth",
                "candidate_id": "S1.truth_weak",
                "chrom": "chr21",
                "start": 20_700_000,
                "end": 22_300_000,
                "state": "gain",
                "event_layer": "autosomal_report",
                "plot_visibility": "review_plot_only",
                "plot_layer_class": "internal_review_event_candidate",
                "plot_support_class": "Z_SUPPORTED_CN_NOT_SUPPORTED",
            },
            {
                "sample_id": "S1",
                "event_id": "evt_weak",
                "candidate_id": "S1.cn_weak",
                "chrom": "chr7",
                "start": 50_000_000,
                "end": 55_000_000,
                "state": "loss",
                "event_layer": "autosomal_report",
                "plot_visibility": "final_report_plot",
                "plot_layer_class": "report_event",
                "plot_support_class": "Z_SUPPORTED_CN_WEAK",
            },
        ]
    )
    support = pd.DataFrame(
        [
            {
                "sample_id": "S1",
                "event_id": "evt_low",
                "candidate_id": "S1.lowconf",
                "support_interpretation_status": "Z_DIRECTION_SUPPORTED",
                "cn_direction_consistency_status": "CN_DIRECTION_NOT_SUPPORTED",
                "same_direction_median_display_ref_z": 8.1,
                "cn_same_direction_fraction": 0.0,
            },
            {
                "sample_id": "S1",
                "event_id": "evt_truth",
                "candidate_id": "S1.truth_weak",
                "support_interpretation_status": "Z_DIRECTION_SUPPORTED",
                "cn_direction_consistency_status": "CN_DIRECTION_NOT_SUPPORTED",
                "same_direction_median_display_ref_z": 5.4,
                "cn_same_direction_fraction": 0.0,
            },
            {
                "sample_id": "S1",
                "event_id": "evt_weak",
                "candidate_id": "S1.cn_weak",
                "support_interpretation_status": "Z_DIRECTION_SUPPORTED",
                "cn_direction_consistency_status": "CN_DIRECTION_WEAK_OR_MIXED",
                "same_direction_median_display_ref_z": -6.2,
                "cn_same_direction_fraction": 0.4,
            },
        ]
    )
    truth_metrics = pd.DataFrame(
        [
            {
                "sample_id": "S1",
                "truth_index": 0,
                "top_candidate_id": "S1.truth_weak",
                "v2_preserved_count": 1,
                "v2_hard_suppressed_count": 0,
            }
        ]
    )
    sample_summary = pd.DataFrame(
        [
            {
                "sample_id": "S1",
                "truth_event_count": 1,
                "truth_preserved_count": 1,
                "FN_count": 0,
            }
        ]
    )

    audit, summary = build_report_table_ablation_audit(
        report_events=report_events,
        plot_manifests=[manifest],
        plot_supports=[support],
        truth_metrics=truth_metrics,
        sample_summary=sample_summary,
        reference_id="h_r0_shadow_ref_20260619",
    )

    actions = dict(zip(audit["candidate_id"], audit["proposed_report_table_action"]))
    assert actions["S1.lowconf"] == "downgrade_to_internal_review_candidate"
    assert actions["S1.truth_weak"] == "retain_report_event_truth_guard"
    assert actions["S1.cn_weak"] == "retain_report_event"
    assert summary["original_report_event_count"] == 3
    assert summary["proposed_internal_review_demotion_count"] == 1
    assert summary["truth_guarded_report_event_count"] == 1
    assert summary["truth_event_count"] == 1
    assert summary["FN_count"] == 0
    assert summary["truth_hard_suppressed_count"] == 0


def test_report_table_ablation_writes_tsv_json_and_markdown_without_truth_metrics(tmp_path):
    report_events = pd.DataFrame(
        [
            {
                "sample_id": "JZ26125843-56-56",
                "candidate_id": "ctx.lowconf",
                "chrom": "chr4",
                "start": 67_500_000,
                "end": 101_250_000,
                "state": "gain",
                "v2_report_layer_class": "report_event",
                "v2_report_visibility": "report_strong_event",
            }
        ]
    )
    manifest = pd.DataFrame(
        [
            {
                "sample_id": "JZ26125843-56-56",
                "event_id": "ctx_evt",
                "candidate_id": "ctx.lowconf",
                "chrom": "chr4",
                "start": 67_500_000,
                "end": 101_250_000,
                "state": "gain",
                "event_layer": "autosomal_report",
                "plot_visibility": "review_plot_only",
                "plot_support_class": "Z_SUPPORTED_CN_NOT_SUPPORTED",
            }
        ]
    )
    audit_tsv = tmp_path / "audit.tsv"
    summary_json = tmp_path / "summary.json"
    report_md = tmp_path / "audit.md"

    audit, summary = build_report_table_ablation_audit(
        report_events=report_events,
        plot_manifests=[manifest],
        plot_supports=[],
        truth_metrics=pd.DataFrame(),
        sample_summary=pd.DataFrame(),
        reference_id="h_r0_shadow_ref_20260619",
        output_audit_tsv=audit_tsv,
        output_summary_json=summary_json,
        output_report_md=report_md,
    )

    assert audit.iloc[0]["truth_guard_status"] == "context_only_no_truth"
    assert summary["status"] == "context_only_no_truth"
    assert summary["TP_FN_FP_status"] == "not_computed_no_truth"
    assert pd.read_csv(audit_tsv, sep="\t").iloc[0]["candidate_id"] == "ctx.lowconf"
    payload = json.loads(summary_json.read_text(encoding="utf-8"))
    assert payload["proposed_internal_review_demotion_count"] == 1
    md = report_md.read_text(encoding="utf-8")
    assert "0615/context cohorts remain burden-only" in md
    assert "not production filtering" in md
