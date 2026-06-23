from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
INDEX_PATH = REPO_ROOT / "docs" / "CURRENT_CONTEXT_INDEX.md"


def test_current_context_index_exists_and_pins_active_context():
    assert INDEX_PATH.exists()
    text = INDEX_PATH.read_text(encoding="utf-8")

    required_tokens = [
        "status: active_current_index",
        "active_handoff: docs/handoff/2026-06-23_2052_repo_archive_status_sync_handoff.md",
        "previous_handoff: docs/handoff/2026-06-23_1949_predict_bam_sample_reportability_qc_handoff.md",
        "active_reference_id: h_r0_shadow_ref_20260619",
        "reference_status: fixed_shadow_baseline_not_production",
        "branch_a_status: burden_phase1_gap2m_materialized_default_unchanged",
        "branch_b_status: v2_report_table_ablation_audit_development_only",
        "branch_s_status: branch_s_review_manifest_visible_not_final",
        "report_status: sample_reportability_qc_remote_validated_0615_materialized",
    ]
    for token in required_tokens:
        assert token in text


def test_current_context_index_points_to_existing_evidence_docs():
    text = INDEX_PATH.read_text(encoding="utf-8")

    expected_paths = [
        "docs/reports/repo_archive_inventory_2026-06-23.md",
        "docs/handoff/2026-06-23_2052_repo_archive_status_sync_handoff.md",
        "docs/reports/predict_bam_and_sample_reportability_qc_2026-06-23.md",
        "docs/handoff/2026-06-23_1949_predict_bam_sample_reportability_qc_handoff.md",
        "docs/reports/p1_p6_result_credibility_audit_2026-06-21.md",
        "docs/reports/branch_b_v2_report_table_ablation_audit_2026-06-22.md",
        "docs/handoff/2026-06-22_2235_report_table_ablation_audit_handoff.md",
        "docs/reports/plot_event_manifest_low_confidence_ablation_2026-06-22.md",
        "docs/handoff/2026-06-22_2145_plot_event_manifest_low_confidence_handoff.md",
        "docs/reports/sex_specific_ref_cnv_plot_2026-06-22.md",
        "docs/handoff/2026-06-22_1938_sex_specific_ref_cnv_plot_handoff.md",
        "docs/reports/h7_h16_reference_cohort_decision_2026-06-20.md",
        "docs/reports/branch_a_validation_h_r0_shadow_2026-06-20.md",
        "docs/reports/branch_a_burden_optimization_phase1_2026-06-21.md",
        "docs/reports/branch_b_v2_truth_safe_filter_2026-06-21.md",
        "docs/reports/branch_b_v2_burden_stratification_2026-06-21.md",
        "docs/reports/branch_b_v2_report_contract_2026-06-21.md",
        "docs/reports/branch_b_v2_report_layer_filter_2026-06-21.md",
        "docs/reports/branch_b_v2_g1_g8_validation_2026-06-21.md",
        "docs/reports/branch_b_v2_pass2_correction_2026-06-21.md",
        "docs/reports/branch_b_v2_lowres_ref_evidence_2026-06-21.md",
        "docs/reports/copy_number_centering_and_sca_fix_2026-06-22.md",
        "docs/reports/copy_number_ratio_cnv_plot_2026-06-22.md",
        "docs/reports/report_main_convergence_cnv_plot_2026-06-22.md",
        "docs/reports/branch_b_s_lowres_integration_2026-06-22.md",
        "docs/reports/branch_s_sca_v2_sex_aware_review_2026-06-22.md",
        "docs/reports/branch_b_v2_report_event_audit_2026-06-21.md",
        "docs/handoff/2026-06-22_1445_cn_centering_branch_s_fix_handoff.md",
        "docs/handoff/2026-06-22_1425_copy_number_ratio_cnv_plot_handoff.md",
        "docs/handoff/2026-06-22_1335_copy_number_cnv_plot_cn_threshold_scatter_handoff.md",
        "docs/handoff/2026-06-22_1305_copy_number_cnv_plot_bin_z_scatter_handoff.md",
        "docs/handoff/2026-06-22_1234_copy_number_cnv_plot_centromere_scatter_handoff.md",
        "docs/handoff/2026-06-22_1155_copy_number_cnv_plot_scatter_background_handoff.md",
        "docs/handoff/2026-06-22_1130_copy_number_cnv_plot_background_fix_handoff.md",
        "docs/handoff/2026-06-22_1100_copy_number_cnv_plot_v2_handoff.md",
        "docs/handoff/2026-06-22_1025_copy_number_cnv_plot_handoff.md",
        "docs/handoff/2026-06-22_0941_cnv_plot_wisecondor_style_handoff.md",
        "docs/handoff/2026-06-22_0437_report_main_cnv_plot_handoff.md",
        "docs/handoff/2026-06-22_0930_lowres_branch_bs_integration_handoff.md",
        "docs/handoff/2026-06-22_0101_branch_s_sca_v2_handoff.md",
        "docs/handoff/2026-06-21_2253_branch_b_v2_lowres_ref_evidence_handoff.md",
        "docs/handoff/2026-06-21_2123_branch_b_v2_pass2_correction_handoff.md",
        "docs/handoff/2026-06-21_2049_branch_b_v2_report_burden_pass2_handoff.md",
        "docs/handoff/2026-06-21_1730_branch_b_v2_report_layer_filter_handoff.md",
        "docs/reports/branch_b_v2_reference_background_and_sca_design_2026-06-20.md",
        "docs/reports/branch_s_p5_report_boundary_2026-06-20.md",
        "docs/reports/p6_report_package_contract_2026-06-20.md",
    ]

    for relative_path in expected_paths:
        assert relative_path in text
        assert (REPO_ROOT / relative_path).exists()


def test_current_context_index_blocks_known_stale_routes():
    text = INDEX_PATH.read_text(encoding="utf-8")

    stale_tokens = [
        "final_disposition",
        "branch_b_keep_event",
        "legacy/current-code Branch B",
        "old N0=0",
        "N1-only matched-negative promotion",
    ]
    for token in stale_tokens:
        assert token in text

    assert "P6 is permanently blocked" not in text
    assert "Branch S can be omitted from reports" not in text
