from __future__ import annotations

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]


def read_text(relative_path):
    return (REPO_ROOT / relative_path).read_text(encoding="utf-8")


def test_predict_dispatcher_exposes_phase1_phase2_actions():
    dispatcher = read_text("scripts/_compat_entry.py")

    assert '"cnv_evidence_ledger": "pgta.predict.branch_b.evidence_ledger"' in dispatcher
    assert '"negative_bank_labels": "pgta.predict.branch_b.negative_bank"' in dispatcher


def test_predict_workflow_has_shadow_ledger_without_report_promotion():
    layout = read_text("rules/predict_layout.smk")
    workflow = read_text("rules/predict_workflow.smk")
    snakefile = read_text("Snakefile")

    assert "CNV_B_EVIDENCE_LEDGER" in layout
    assert "CNV_B_EVIDENCE_SUMMARY" in layout
    assert "rule cnv_branch_b_evidence_ledger" in workflow
    assert "SCRIPT_CNV_EVIDENCE_LEDGER_ACTION" in workflow
    assert "CNV_B_EVIDENCE_LEDGER" in snakefile

    report_rule = workflow.split('if "cnv_report" in AVAILABLE_TARGETS:', 1)[1]
    assert "CNV_B_EVIDENCE_LEDGER" not in report_rule


def test_negative_bank_rule_is_config_gated_and_not_automatic_n0():
    layout = read_text("rules/predict_layout.smk")
    workflow = read_text("rules/predict_workflow.smk")

    assert "CNV_NEGATIVE_BANK_SAMPLES_TSV" in layout
    assert "CNV_NEGATIVE_BANK_LABELS" in layout
    assert "rule cnv_negative_bank_labels" in workflow
    assert "if CNV_NEGATIVE_BANK_SAMPLES_TSV" in workflow
    assert "SCRIPT_NEGATIVE_BANK_LABELS_ACTION" in workflow


def test_phase3_matched_negative_is_shadow_only_and_report_safe():
    dispatcher = read_text("scripts/_compat_entry.py")
    layout = read_text("rules/predict_layout.smk")
    workflow = read_text("rules/predict_workflow.smk")
    snakefile = read_text("Snakefile")

    assert '"matched_negative_percentile": "pgta.predict.branch_b.matched_negative"' in dispatcher
    assert "CNV_B_MATCHED_NEGATIVE" in layout
    assert "CNV_B_MATCHED_NEGATIVE_SUMMARY" in layout
    assert "rule cnv_branch_b_matched_negative" in workflow
    assert "SCRIPT_MATCHED_NEGATIVE_ACTION" in workflow
    assert "CNV_B_MATCHED_NEGATIVE" in snakefile

    report_rule = workflow.split('if "cnv_report" in AVAILABLE_TARGETS:', 1)[1]
    assert "CNV_B_MATCHED_NEGATIVE" not in report_rule


def test_phase4_v2_classifier_is_shadow_only_and_report_safe():
    dispatcher = read_text("scripts/_compat_entry.py")
    layout = read_text("rules/predict_layout.smk")
    workflow = read_text("rules/predict_workflow.smk")
    target_assembly = read_text("rules/target_assembly.smk")
    snakefile = read_text("Snakefile")

    assert '"branch_b_v2_classifier": "pgta.predict.branch_b.v2_classifier"' in dispatcher
    assert "CNV_B_V2_CLASSIFIER" in layout
    assert "CNV_B_V2_CLASSIFIER_SUMMARY" in layout
    assert "rule cnv_branch_b_v2_classifier" in workflow
    assert "SCRIPT_BRANCH_B_V2_CLASSIFIER_ACTION" in workflow
    assert "CNV_B_V2_CLASSIFIER" in target_assembly
    assert "CNV_B_V2_CLASSIFIER" in snakefile

    report_rule = workflow.split('if "cnv_report" in AVAILABLE_TARGETS:', 1)[1]
    assert "CNV_B_V2_CLASSIFIER" not in report_rule


def test_phase5_branch_s_is_shadow_only_and_report_safe():
    dispatcher = read_text("scripts/_compat_entry.py")
    layout = read_text("rules/predict_layout.smk")
    workflow = read_text("rules/predict_workflow.smk")
    target_assembly = read_text("rules/target_assembly.smk")
    snakefile = read_text("Snakefile")

    assert '"branch_s_shadow": "pgta.predict.branch_s"' in dispatcher
    assert "CNV_BRANCH_S_EVIDENCE" in layout
    assert "CNV_BRANCH_S_SCORES" in layout
    assert "CNV_BRANCH_S_SUMMARY" in layout
    assert "rule cnv_branch_s_shadow" in workflow
    assert "SCRIPT_BRANCH_S_SHADOW_ACTION" in workflow
    assert "CNV_BRANCH_S_EVIDENCE" in target_assembly
    assert "CNV_BRANCH_S_EVIDENCE" in snakefile

    report_rule = workflow.split('if "cnv_report" in AVAILABLE_TARGETS:', 1)[1]
    assert "CNV_BRANCH_S_EVIDENCE" not in report_rule
    assert "CNV_BRANCH_S_SCORES" not in report_rule
    assert "CNV_BRANCH_S_SUMMARY" not in report_rule
