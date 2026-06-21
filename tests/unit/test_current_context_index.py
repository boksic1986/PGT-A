from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
INDEX_PATH = REPO_ROOT / "docs" / "CURRENT_CONTEXT_INDEX.md"


def test_current_context_index_exists_and_pins_active_context():
    assert INDEX_PATH.exists()
    text = INDEX_PATH.read_text(encoding="utf-8")

    required_tokens = [
        "status: active_current_index",
        "active_handoff: docs/handoff/2026-06-21_1810_g1_g8_current_scheme_validation_handoff.md",
        "previous_handoff: docs/handoff/2026-06-21_1730_branch_b_v2_report_layer_filter_handoff.md",
        "active_reference_id: h_r0_shadow_ref_20260619",
        "reference_status: fixed_shadow_baseline_not_production",
        "branch_a_status: burden_phase1_gap2m_materialized_default_unchanged",
        "branch_b_status: v2_report_layer_filter_materialized_development_only",
        "branch_s_status: review_reportable_with_limitations",
        "report_status: branch_b_v2_report_layer_filter_integrated_development_only",
    ]
    for token in required_tokens:
        assert token in text


def test_current_context_index_points_to_existing_evidence_docs():
    text = INDEX_PATH.read_text(encoding="utf-8")

    expected_paths = [
        "docs/reports/p1_p6_result_credibility_audit_2026-06-21.md",
        "docs/reports/h7_h16_reference_cohort_decision_2026-06-20.md",
        "docs/reports/branch_a_validation_h_r0_shadow_2026-06-20.md",
        "docs/reports/branch_a_burden_optimization_phase1_2026-06-21.md",
        "docs/reports/branch_b_v2_truth_safe_filter_2026-06-21.md",
        "docs/reports/branch_b_v2_burden_stratification_2026-06-21.md",
        "docs/reports/branch_b_v2_report_contract_2026-06-21.md",
        "docs/reports/branch_b_v2_report_layer_filter_2026-06-21.md",
        "docs/reports/branch_b_v2_g1_g8_validation_2026-06-21.md",
        "docs/handoff/2026-06-21_1810_g1_g8_current_scheme_validation_handoff.md",
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
