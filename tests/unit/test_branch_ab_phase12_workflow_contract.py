from __future__ import annotations

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]


def read_text(relative_path):
    return (REPO_ROOT / relative_path).read_text(encoding="utf-8")


def test_predict_dispatcher_exposes_phase1_phase2_actions():
    dispatcher = read_text("scripts/_compat_entry.py")

    assert '"cnv_evidence_ledger": "pgta.predict.branch_b.evidence_ledger"' in dispatcher
    assert '"branch_b_v2_benchmark": "pgta.predict.branch_b.v2_benchmark"' in dispatcher
    assert '"branch_b_lowres_evidence": "pgta.predict.branch_b.lowres_evidence"' in dispatcher
    assert '"branch_b_ref_stability": "pgta.predict.branch_b.ref_stability"' in dispatcher
    assert '"negative_bank_labels": "pgta.predict.branch_b.negative_bank"' in dispatcher
    assert '"reference_candidate_audit": "pgta.reference.audit"' in dispatcher


def test_reference_audit_target_is_branch_a_level_and_report_safe():
    layout = read_text("rules/predict_layout.smk")
    workflow = read_text("rules/predict_workflow.smk")
    pipeline_modes = read_text("rules/pipeline_modes.smk")
    target_assembly = read_text("rules/target_assembly.smk")
    snakefile = read_text("Snakefile")

    assert "CNV_REFERENCE_AUDIT_TSV" in layout
    assert "CNV_REFERENCE_AUDIT_SUMMARY" in layout
    assert '"reference_audit"' in pipeline_modes
    assert "CNV_REFERENCE_AUDIT_TSV" in target_assembly
    assert "rule cnv_reference_candidate_audit" in workflow
    assert "SCRIPT_REFERENCE_CANDIDATE_AUDIT_ACTION" in workflow
    assert "rule reference_audit" in snakefile

    audit_rule = workflow.split("rule cnv_reference_candidate_audit:", 1)[1].split("if CNV_POSTPROCESS_ENABLE_BRANCH_B:", 1)[0]
    assert "CNV_A_CANDIDATES" in audit_rule
    assert "CNV_QC_TSV" in audit_rule
    assert "CNV_B_FINAL_EVENTS" not in audit_rule
    assert "CNV_B_ARTIFACT_SUMMARY" not in audit_rule
    assert "CNV_B_MATCHED_NEGATIVE" not in audit_rule

    report_rule = workflow.split('if "cnv_report" in AVAILABLE_TARGETS:', 1)[1]
    assert "CNV_REFERENCE_AUDIT_TSV" not in report_rule


def test_branch_a_validation_target_is_branch_a_only_and_report_safe():
    dispatcher = read_text("scripts/_compat_entry.py")
    layout = read_text("rules/predict_layout.smk")
    workflow = read_text("rules/predict_workflow.smk")
    pipeline_modes = read_text("rules/pipeline_modes.smk")
    target_assembly = read_text("rules/target_assembly.smk")
    snakefile = read_text("Snakefile")

    assert '"branch_a_validation": "pgta.predict.branch_a_validation"' in dispatcher
    assert "SCRIPT_BRANCH_A_VALIDATION_ACTION" in workflow
    assert "CNV_BRANCH_A_VALIDATION_SAMPLE_SUMMARY" in layout
    assert '"branch_a_validation"' in pipeline_modes
    assert "CNV_BRANCH_A_VALIDATION_SAMPLE_SUMMARY" in target_assembly
    assert "rule branch_a_validation" in snakefile

    validation_rule = workflow.split("rule cnv_branch_a_validation:", 1)[1].split("if \"reference_audit\" in AVAILABLE_TARGETS:", 1)[0]
    assert "CNV_A_CANDIDATES" in validation_rule
    assert "CNV_B_FINAL_EVENTS" not in validation_rule
    assert "CNV_B_ARTIFACT_SUMMARY" not in validation_rule
    assert "CNV_B_MATCHED_NEGATIVE" not in validation_rule

    report_rule = workflow.split('if "cnv_report" in AVAILABLE_TARGETS:', 1)[1]
    assert "CNV_BRANCH_A_VALIDATION_SAMPLE_SUMMARY" not in report_rule


def test_branch_a_assembly_exposes_configured_gap_without_branch_b_filters():
    layout = read_text("rules/predict_layout.smk")
    workflow = read_text("rules/predict_workflow.smk")

    assert "CNV_BRANCH_A_CFG" in layout
    assert "CNV_BRANCH_A_MERGE_GAP_BP" in layout
    assert "CNV_BRANCH_A_STRONG_Z" in layout
    assert "CNV_BRANCH_A_OUTPUT_DIR" in layout
    assert "CNV_BRANCH_A_VALIDATION_DIR" in layout
    assert "CNV_BRANCH_A_LOG_DIR" in layout

    assembly_rule = workflow.split("rule a_branch_candidate_assembly:", 1)[1].split("rule cnv_branch_a_validation:", 1)[0]
    assert "merge_gap_bp=CNV_BRANCH_A_MERGE_GAP_BP" in assembly_rule
    assert "strong_z=CNV_BRANCH_A_STRONG_Z" in assembly_rule
    assert '"--merge-gap-bp", str(params.merge_gap_bp)' in assembly_rule
    assert '"--strong-z", str(params.strong_z)' in assembly_rule
    assert "CNV_B_FINAL_EVENTS" not in assembly_rule
    assert "CNV_B_ARTIFACT_SUMMARY" not in assembly_rule
    assert "CNV_B_MATCHED_NEGATIVE" not in assembly_rule


def test_branch_a_output_dirs_are_configurable_without_overwriting_default_candidates():
    layout = read_text("rules/predict_layout.smk")
    workflow = read_text("rules/predict_workflow.smk")

    assert 'CNV_BRANCH_A_CFG.get("output_dir"' in layout
    assert 'CNV_BRANCH_A_CFG.get("validation_dir"' in layout
    assert 'CNV_BRANCH_A_CFG.get("log_dir"' in layout
    assert 'Path(CNV_BRANCH_A_OUTPUT_DIR) / "{sample}.candidate_events.tsv"' in layout
    assert 'Path(CNV_BRANCH_A_VALIDATION_DIR) / "summary.json"' in layout
    assert 'Path(CNV_BRANCH_A_LOG_DIR) / "{sample}.a_branch_candidates.log"' in workflow
    assert 'Path(CNV_BRANCH_A_LOG_DIR) / "branch_a_validation.log"' in workflow
    assert '"--a-candidates-dir", CNV_BRANCH_A_OUTPUT_DIR' in workflow


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
    assert "CNV_B_EVIDENCE_SUMMARY" in report_rule


def test_branch_b_evidence_target_is_p3_only_and_report_safe():
    workflow = read_text("rules/predict_workflow.smk")
    pipeline_modes = read_text("rules/pipeline_modes.smk")
    target_assembly = read_text("rules/target_assembly.smk")
    snakefile = read_text("Snakefile")

    assert '"branch_b_evidence"' in pipeline_modes
    assert "CNV_B_EVIDENCE_LEDGER" in target_assembly
    assert "rule branch_b_evidence" in snakefile

    evidence_target = target_assembly.split('if "branch_b_evidence" in REQUESTED_TARGETS', 1)[1].split('if "branch_b_v2_benchmark" in REQUESTED_TARGETS', 1)[0]
    assert "CNV_B_EVIDENCE_LEDGER" in evidence_target
    assert "CNV_B_EVIDENCE_SUMMARY" in evidence_target
    assert "CNV_B_V2_CLASSIFIER" not in evidence_target
    assert "CNV_BRANCH_S_EVIDENCE" not in evidence_target
    assert "CNV_REPORT_TSV" not in evidence_target

    report_rule = workflow.split('if "cnv_report" in AVAILABLE_TARGETS:', 1)[1]
    assert "CNV_B_EVIDENCE_LEDGER" not in report_rule
    assert "CNV_B_EVIDENCE_SUMMARY" in report_rule


def test_branch_b_v2_benchmark_is_v2_only_and_report_safe():
    dispatcher = read_text("scripts/_compat_entry.py")
    layout = read_text("rules/predict_layout.smk")
    workflow = read_text("rules/predict_workflow.smk")
    pipeline_modes = read_text("rules/pipeline_modes.smk")
    target_assembly = read_text("rules/target_assembly.smk")
    snakefile = read_text("Snakefile")

    assert '"branch_b_v2_benchmark": "pgta.predict.branch_b.v2_benchmark"' in dispatcher
    assert '"branch_b_v2_benchmark"' in pipeline_modes
    assert "CNV_POSTPROCESS_CFG.get(\"output_dir\"" in layout
    assert "CNV_B_V2_BENCHMARK_TRUTH_METRICS" in layout
    assert "CNV_B_V2_BENCHMARK_SAMPLE_SUMMARY" in layout
    assert "CNV_B_V2_BENCHMARK_FILTERED_EVENTS" in layout
    assert "CNV_B_V2_BENCHMARK_FILTERED_EVENTS_JSON" in layout
    assert "CNV_B_V2_BENCHMARK_REPORT_EVENTS" in layout
    assert "CNV_B_V2_BENCHMARK_REPORT_EVENTS_JSON" in layout
    assert "CNV_B_V2_BENCHMARK_SUMMARY" in layout
    assert "rule cnv_branch_b_v2_benchmark" in workflow
    assert "rule branch_b_v2_benchmark" in snakefile
    assert "SCRIPT_BRANCH_B_V2_BENCHMARK_ACTION" in workflow

    benchmark_rule = workflow.split("rule cnv_branch_b_v2_benchmark:", 1)[1].split("rule cnv_branch_s_shadow:", 1)[0]
    assert "CNV_B_V2_CLASSIFIER" in benchmark_rule
    assert "filtered_events=CNV_B_V2_BENCHMARK_FILTERED_EVENTS" in benchmark_rule
    assert "filtered_events_json=CNV_B_V2_BENCHMARK_FILTERED_EVENTS_JSON" in benchmark_rule
    assert "report_events=CNV_B_V2_BENCHMARK_REPORT_EVENTS" in benchmark_rule
    assert "report_events_json=CNV_B_V2_BENCHMARK_REPORT_EVENTS_JSON" in benchmark_rule
    assert "--output-filtered-events" in benchmark_rule
    assert "--output-filtered-events-json" in benchmark_rule
    assert "--output-report-events" in benchmark_rule
    assert "--output-report-events-json" in benchmark_rule
    assert "CNV_B_FINAL_EVENTS" not in benchmark_rule
    assert "CNV_B_ARTIFACT_SUMMARY" not in benchmark_rule
    assert "CNV_B_MATCHED_NEGATIVE" not in benchmark_rule

    benchmark_target = target_assembly.split('if "branch_b_v2_benchmark" in REQUESTED_TARGETS', 1)[1].split('if "reference_audit" in REQUESTED_TARGETS', 1)[0]
    assert "CNV_B_V2_CLASSIFIER" in benchmark_target
    assert "CNV_B_V2_BENCHMARK_SUMMARY" in benchmark_target
    assert "CNV_B_V2_BENCHMARK_FILTERED_EVENTS" in benchmark_target
    assert "CNV_B_V2_BENCHMARK_FILTERED_EVENTS_JSON" in benchmark_target
    assert "CNV_B_V2_BENCHMARK_REPORT_EVENTS" in benchmark_target
    assert "CNV_B_V2_BENCHMARK_REPORT_EVENTS_JSON" in benchmark_target
    assert "CNV_REPORT_TSV" not in benchmark_target

    report_rule = workflow.split('if "cnv_report" in AVAILABLE_TARGETS:', 1)[1]
    assert "CNV_B_V2_BENCHMARK_SUMMARY" in report_rule
    assert "CNV_B_V2_BENCHMARK_SAMPLE_SUMMARY" in report_rule
    assert 'if "branch_b_v2_benchmark" in AVAILABLE_TARGETS' in report_rule
    assert "--branch-b-v2-benchmark-summary" in report_rule
    assert "--branch-b-v2-sample-summary" in report_rule
    v2_report_input_block = report_rule.split("branch_b_v2_benchmark_summary", 1)[1].split("evaluation_summary", 1)[0]
    assert "REQUESTED_TARGETS" not in v2_report_input_block
    assert "CNV_B_FINAL_EVENTS" not in v2_report_input_block


def test_cnv_plot_uses_v2_report_events_and_writes_plot_bins():
    layout = read_text("rules/predict_layout.smk")
    workflow = read_text("rules/predict_workflow.smk")

    assert "CNV_B_PLOT_BINS_TSV" in layout
    assert "{sample}.plot_bins.tsv" in layout

    plot_rule = workflow.split("rule cnv_branch_ab_plot:", 1)[1].split("if CNV_NEGATIVE_BANK_SAMPLES_TSV:", 1)[0]
    assert "events=CNV_B_V2_BENCHMARK_REPORT_EVENTS" in plot_rule
    assert "CNV_B_FINAL_EVENTS" not in plot_rule
    assert "bins_tsv=CNV_B_PLOT_BINS_TSV" in plot_rule
    assert "--output-bins-tsv" in plot_rule


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
    assert "CNV_NEGATIVE_BANK_SHADOW_BACKGROUND_LABELS" in layout
    assert "--shadow-background-label" in workflow
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
    assert "CNV_BRANCH_S_SUMMARY" in report_rule


def test_branch_s_review_target_is_p5_only_and_report_safe():
    workflow = read_text("rules/predict_workflow.smk")
    pipeline_modes = read_text("rules/pipeline_modes.smk")
    target_assembly = read_text("rules/target_assembly.smk")
    snakefile = read_text("Snakefile")

    assert '"branch_s_review"' in pipeline_modes
    assert "rule branch_s_review" in snakefile

    branch_s_target = target_assembly.split('if "branch_s_review" in REQUESTED_TARGETS', 1)[1].split('if "cnv_eval" in REQUESTED_TARGETS', 1)[0]
    assert "CNV_BRANCH_S_EVIDENCE" in branch_s_target
    assert "CNV_BRANCH_S_SCORES" in branch_s_target
    assert "CNV_BRANCH_S_SUMMARY" in branch_s_target
    assert "CNV_B_V2_CLASSIFIER" not in branch_s_target
    assert "CNV_B_MATCHED_NEGATIVE" not in branch_s_target
    assert "CNV_REPORT_TSV" not in branch_s_target

    report_rule = workflow.split('if "cnv_report" in AVAILABLE_TARGETS:', 1)[1]
    assert "CNV_BRANCH_S_EVIDENCE" not in report_rule
    assert "CNV_BRANCH_S_SCORES" not in report_rule
    assert "CNV_BRANCH_S_SUMMARY" in report_rule


def test_p6_report_package_consumes_review_summaries_only():
    layout = read_text("rules/predict_layout.smk")
    workflow = read_text("rules/predict_workflow.smk")

    assert "CNV_REPORT_BRANCH_A_VALIDATION_SUMMARIES" in layout

    report_rule = workflow.split('if "cnv_report" in AVAILABLE_TARGETS:', 1)[1]
    assert "CNV_B_EVIDENCE_SUMMARY" in report_rule
    assert "CNV_BRANCH_S_SUMMARY" in report_rule
    assert "CNV_REPORT_BRANCH_A_VALIDATION_SUMMARIES" in report_rule
    assert "--reference-id" in report_rule
    assert "--wisecondorx-predict-command" in report_rule
    assert "--branch-a-validation-summary" in report_rule
    assert "--branch-b-evidence-summary" in report_rule
    assert "--branch-s-summary" in report_rule
    assert "--branch-b-v2-benchmark-summary" in report_rule
    assert "--branch-b-v2-sample-summary" in report_rule
    assert "CNV_B_EVIDENCE_LEDGER" not in report_rule
    assert "CNV_BRANCH_S_EVIDENCE" not in report_rule
    assert "CNV_BRANCH_S_SCORES" not in report_rule
    assert "CNV_B_MATCHED_NEGATIVE" not in report_rule
    assert "CNV_B_V2_CLASSIFIER" not in report_rule
    assert "final_disposition" not in report_rule
    assert "branch_b_keep_event" not in report_rule


def test_lowres_evidence_is_optional_v2_ledger_context_not_report_promotion():
    layout = read_text("rules/predict_layout.smk")
    workflow = read_text("rules/predict_workflow.smk")
    target_assembly = read_text("rules/target_assembly.smk")

    assert "CNV_LOWRES_EVIDENCE_CFG" in layout
    assert "CNV_B_REF_STABILITY_BINS" in layout
    assert "CNV_B_LOWRES_EVIDENCE" in layout
    assert "reference_npz_glob" in layout
    assert "lowres_evidence.enable=true requires lowres_evidence.reference_npz or reference_npz_glob" in layout
    assert "reference_sample_ids must match lowres_evidence.reference_npz length" in layout
    assert "SCRIPT_BRANCH_B_LOWRES_EVIDENCE_ACTION" in workflow
    assert "SCRIPT_BRANCH_B_REF_STABILITY_ACTION" in workflow
    assert "rule cnv_branch_b_ref_stability" in workflow
    assert "rule cnv_branch_b_lowres_evidence" in workflow

    lowres_rule = workflow.split("rule cnv_branch_b_lowres_evidence:", 1)[1].split("rule cnv_branch_b_v2_classifier:", 1)[0]
    assert "CNV_B_LOWRES_EVIDENCE" in lowres_rule
    assert "--lowres-2mb-events" in lowres_rule
    assert "--lowres-3mb-events" in lowres_rule
    assert "--ref-stability-events" in lowres_rule
    assert "CNV_B_FINAL_EVENTS" not in lowres_rule
    assert "final_disposition" not in lowres_rule
    assert "branch_b_keep_event" not in lowres_rule

    classifier_rule = workflow.split("rule cnv_branch_b_v2_classifier:", 1)[1].split("if \"branch_b_v2_benchmark\" in AVAILABLE_TARGETS:", 1)[0]
    assert "CNV_B_LOWRES_EVIDENCE" in classifier_rule

    branch_s_rule = workflow.split("rule cnv_branch_s_shadow:", 1)[1].split("if CNV_MOSAIC_FRACTION_TRUTH_TSV", 1)[0]
    assert "CNV_LOWRES_EVIDENCE_ENABLE" in branch_s_rule
    assert "--lowres-2mb-events" in branch_s_rule
    assert "--lowres-3mb-events" in branch_s_rule

    branch_b_target = target_assembly.split('if "branch_b_v2_benchmark" in REQUESTED_TARGETS', 1)[1].split('if "reference_audit" in REQUESTED_TARGETS', 1)[0]
    assert "CNV_B_LOWRES_EVIDENCE" in branch_b_target

    report_rule = workflow.split('if "cnv_report" in AVAILABLE_TARGETS:', 1)[1]
    assert "CNV_B_LOWRES_EVIDENCE" not in report_rule
    assert "CNV_B_REF_STABILITY_BINS" not in report_rule
