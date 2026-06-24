PIPELINE_TARGETS = set(PIPELINE_CFG.get("targets", ["mapping", "reference", "cnv"]))

if PIPELINE_MODE_RAW in {"reference", "build_reference", "build_ref"}:
    PIPELINE_MODE = "build_ref"
elif PIPELINE_MODE_RAW in {"predict", "prediction"}:
    PIPELINE_MODE = "predict"
elif PIPELINE_MODE_RAW:
    raise ValueError(
        "Unsupported pipeline.mode: "
        + PIPELINE_MODE_RAW
        + ". Use build_ref or predict."
    )
else:
    ref_targets = {"reference_qc", "reference"}
    ref_targets.add("baseline_qc")
    predict_targets = {"predict_bam_qc", "cnv_qc", "cnv", "branch_a_validation", "branch_b_evidence", "branch_b_v2_benchmark", "branch_b_v2_report_ablation", "branch_s_review", "cnv_event_annotation", "cnv_sample_report_qc", "cnv_eval", "cnv_ml", "cnv_benchmark", "cnv_report"}
    if PIPELINE_TARGETS & ref_targets and PIPELINE_TARGETS & predict_targets:
        raise ValueError(
            "pipeline.targets mixes reference and predict stages without pipeline.mode. "
            "Split configs into build_ref and predict entry YAMLs."
        )
    PIPELINE_MODE = "build_ref" if PIPELINE_TARGETS & ref_targets else "predict"

REQUESTED_TARGETS = set(PIPELINE_TARGETS)
AVAILABLE_TARGETS = {"mapping", "metadata", "baseline_qc", "reference"}
if TUNING_ENABLED:
    AVAILABLE_TARGETS.add("reference_qc")
if CNV_ENABLED:
    AVAILABLE_TARGETS.update({"predict_bam_qc", "cnv_qc", "cnv", "reference_audit", "branch_a_validation", "branch_b_evidence", "branch_b_v2_benchmark", "branch_b_v2_report_ablation", "branch_s_review", "cnv_event_annotation", "cnv_sample_report_qc", "cnv_eval", "cnv_ml", "cnv_benchmark", "cnv_report"})

UNKNOWN_TARGETS = sorted(REQUESTED_TARGETS - AVAILABLE_TARGETS)
if UNKNOWN_TARGETS:
    raise ValueError(
        "Unsupported pipeline targets: "
        + ",".join(UNKNOWN_TARGETS)
        + f". Available: {','.join(sorted(AVAILABLE_TARGETS))}"
    )

MODE_TARGETS = {
    "build_ref": {"mapping", "metadata", "baseline_qc", "reference_qc", "reference"},
    "predict": {"mapping", "metadata", "predict_bam_qc", "cnv_qc", "cnv", "reference_audit", "branch_a_validation", "branch_b_evidence", "branch_b_v2_benchmark", "branch_b_v2_report_ablation", "branch_s_review", "cnv_event_annotation", "cnv_sample_report_qc", "cnv_eval", "cnv_ml", "cnv_benchmark", "cnv_report"},
}
INVALID_MODE_TARGETS = sorted(REQUESTED_TARGETS - MODE_TARGETS[PIPELINE_MODE])
if INVALID_MODE_TARGETS:
    raise ValueError(
        "pipeline.mode="
        + PIPELINE_MODE
        + " does not allow targets: "
        + ",".join(INVALID_MODE_TARGETS)
        + ". Allowed targets: "
        + ",".join(sorted(MODE_TARGETS[PIPELINE_MODE]))
    )

if ("reference" in REQUESTED_TARGETS or "reference_qc" in REQUESTED_TARGETS) and SEX_SPECIFIC_REF_CONFIGURED and not BUILD_REF_BY_SEX_ENABLED:
    raise ValueError(
        "Sex-specific reference build requires build_reference.groups.XX and build_reference.groups.XY in the config."
    )
if "reference" in REQUESTED_TARGETS and not REF_TARGET_FILES:
    raise ValueError(
        "No reference outputs configured. Set core.wisecondorx.reference_output or configure build_reference.groups with reference_output_by_sex."
    )
if "baseline_qc" in REQUESTED_TARGETS and len(BASELINE_SAMPLE_IDS) < 2:
    raise ValueError("baseline_qc requires at least 2 baseline/reference samples.")
if ("cnv_qc" in REQUESTED_TARGETS or "cnv" in REQUESTED_TARGETS) and PREDICT_BY_SEX_ENABLED and not GENDER_REF_OUTPUT:
    raise ValueError("core.wisecondorx.gender_reference_output is required for sex-specific prediction.")
if "branch_b_evidence" in REQUESTED_TARGETS and CNV_ENABLED and not CNV_POSTPROCESS_ENABLE_BRANCH_B:
    raise ValueError("branch_b_evidence requires core.wisecondorx.cnv.postprocess.enable_branch_b=true.")
if "branch_b_v2_benchmark" in REQUESTED_TARGETS and CNV_ENABLED and not CNV_POSTPROCESS_ENABLE_BRANCH_B:
    raise ValueError("branch_b_v2_benchmark requires core.wisecondorx.cnv.postprocess.enable_branch_b=true.")
if "branch_b_v2_report_ablation" in REQUESTED_TARGETS and CNV_ENABLED and not CNV_POSTPROCESS_ENABLE_BRANCH_B:
    raise ValueError("branch_b_v2_report_ablation requires core.wisecondorx.cnv.postprocess.enable_branch_b=true.")
if "branch_s_review" in REQUESTED_TARGETS and CNV_ENABLED and not CNV_POSTPROCESS_ENABLE_BRANCH_B:
    raise ValueError("branch_s_review requires core.wisecondorx.cnv.postprocess.enable_branch_b=true.")
if "cnv_event_annotation" in REQUESTED_TARGETS and CNV_ENABLED and not CNV_POSTPROCESS_ENABLE_BRANCH_B:
    raise ValueError("cnv_event_annotation requires core.wisecondorx.cnv.postprocess.enable_branch_b=true.")
if PIPELINE_MODE == "build_ref" and "build_reference" not in config:
    raise ValueError("build_ref config must define build_reference.")
if PIPELINE_MODE == "predict" and "build_reference" in config and BUILD_REF_CFG.get("groups"):
    raise ValueError(
        "predict config should not carry build_reference.groups. Use a dedicated predict YAML generated from predict_samples.yaml."
    )
