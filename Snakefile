configfile: "build_samples.yaml"

from pathlib import Path
import sys
import yaml


def merge_config(base, override):
    merged = dict(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = merge_config(merged[key], value)
        else:
            merged[key] = value
    return merged


if "base_config" in config:
    base_cfg_path = Path(config["base_config"])
    if not base_cfg_path.is_absolute():
        base_cfg_path = (Path(workflow.basedir) / base_cfg_path).resolve()
    with open(base_cfg_path, "r", encoding="utf-8") as base_cfg_handle:
        base_cfg = yaml.safe_load(base_cfg_handle) or {}
    override_cfg = {k: v for k, v in config.items() if k != "base_config"}
    override_mode = str(override_cfg.get("pipeline", {}).get("mode", "")).strip().lower()
    if override_mode == "predict":
        # Drop default-config carryover so predict uses the dedicated base_config sample list and settings.
        override_cfg = dict(override_cfg)
        override_cfg.pop("build_reference", None)
        override_cfg.pop("samples", None)
    config = merge_config(base_cfg, override_cfg)


SAMPLES = sorted(config["samples"].keys())
assert SAMPLES, "No samples found in config['samples']."

PIPELINE_ROOT = Path(workflow.basedir).resolve()
sys.path.insert(0, str(PIPELINE_ROOT))

include: "rules/script_entrypoints.smk"

PROJECT = Path(config["core"]["project_path"])
WISE_CFG = config["core"]["wisecondorx"]
TUNE_CFG = WISE_CFG.get("tuning", {})
TUNE_PCA_CFG = TUNE_CFG.get("pca", {})
TUNE_QC_CFG = TUNE_CFG.get("qc", {})
PREFILTER_CFG = WISE_CFG.get("reference_prefilter", {})
CNV_CFG = WISE_CFG.get("cnv", {})
CNV_QC_CFG = CNV_CFG.get("qc", {})
CNV_SEX_CFG = CNV_CFG.get("sex_calling", {})
CNV_POSTPROCESS_CFG = CNV_CFG.get("postprocess", {})
CNV_POSTPROCESS_ANNOTATION_CFG = CNV_POSTPROCESS_CFG.get("annotation", {})
CNV_CORRECTION_CFG = CNV_POSTPROCESS_CFG.get("correction", {})
CNV_CALLING_CFG = CNV_POSTPROCESS_CFG.get("calling", {})
CNV_CALIBRATION_CFG = CNV_POSTPROCESS_CFG.get("calibration", {})
CNV_MOSAIC_CFG = CNV_POSTPROCESS_CFG.get("mosaic_fraction", {})
CNV_ARTIFACT_CFG = CNV_POSTPROCESS_CFG.get("artifact_rules", {})
CNV_ML_CFG = CNV_CFG.get("ml", {})
CNV_EVALUATION_CFG = CNV_CFG.get("evaluation", {})
CNV_BENCHMARK_CFG = CNV_CFG.get("benchmark", {})
CNV_REPORT_CFG = CNV_CFG.get("report", {})
PIPELINE_CFG = config.get("pipeline", {})
BUILD_REF_CFG = config.get("build_reference", {})
BASELINE_QC_CFG = config.get("baseline_qc", {})
BASELINE_QC_THRESHOLDS = BASELINE_QC_CFG.get("thresholds", {})
BASELINE_QC_GC_CORRECTION_CFG = BASELINE_QC_CFG.get("gc_correction", {})
REFERENCE_ASSETS_CFG = config.get("reference_assets", {})
REFERENCE_PACKAGE_CFG = config.get("reference_package", {})
QC_FRAMEWORK_CFG = config.get("structured_qc", {})
PIPELINE_MODE_RAW = str(PIPELINE_CFG.get("mode", "")).strip().lower()


def format_binsize_label(bin_size_bp):
    value = int(bin_size_bp)
    if value % 1000000 == 0:
        return f"{value // 1000000}mb"
    if value % 1000 == 0:
        return f"{value // 1000}kb"
    return f"{value}bp"


project_path = lambda *parts: str(PROJECT.joinpath(*parts))
resolve_path = lambda path_value: str(Path(path_value)) if Path(path_value).is_absolute() else str(PROJECT / Path(path_value))


FASTP_R1 = project_path("fastp", "{sample}.R1.clean.fastq.gz")
FASTP_R2 = project_path("fastp", "{sample}.R2.clean.fastq.gz")
FASTP_HTML = project_path("fastp", "{sample}.fastp.html")
FASTP_JSON = project_path("fastp", "{sample}.fastp.json")
SORTED_BAM = project_path("mapping", "{sample}.sorted.bam")
SORTED_BAI = project_path("mapping", "{sample}.sorted.bam.bai")

include: "rules/reference_layout.smk"
include: "rules/predict_layout.smk"
include: "rules/runtime_layout.smk"
include: "rules/pipeline_modes.smk"


def read_int_from_file(path_value):
    path = Path(path_value)
    if not path.exists():
        raise FileNotFoundError(f"Required file not found: {path}")
    return int(path.read_text(encoding="utf-8").strip())


def read_best_binsize_from_yaml(path_value, fallback_value):
    path = Path(path_value)
    if not path.exists():
        return int(fallback_value)
    payload = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    return int(payload.get("best_binsize", payload.get("binsize", fallback_value)))


def load_gender_result(gender_tsv_path):
    lines = [line.strip() for line in Path(gender_tsv_path).read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(lines) < 2:
        raise ValueError(f"Invalid gender TSV: {gender_tsv_path}")
    header = lines[0].split("\t")
    values = lines[1].split("\t")
    row = dict(zip(header, values))
    sex_call = row.get("sex_call", "").strip().upper()
    wise_gender = row.get("wise_gender", "").strip().upper()
    if sex_call not in {"XX", "XY"}:
        raise ValueError(f"Unsupported sex_call in {gender_tsv_path}: {sex_call}")
    if wise_gender not in {"F", "M"}:
        raise ValueError(f"Unsupported wise_gender in {gender_tsv_path}: {wise_gender}")
    return row


def select_predict_reference(gender_tsv_path):
    row = load_gender_result(gender_tsv_path)
    return REF_OUTPUTS_BY_SEX[row["sex_call"]]


def select_predict_gender(gender_tsv_path):
    row = load_gender_result(gender_tsv_path)
    predict_gender = row.get("predict_gender", "").strip().upper()
    if predict_gender in {"F", "M"}:
        return predict_gender
    return "F" if row["sex_call"] == "XX" else "M"

include: "rules/target_assembly.smk"


rule all:
    input:
        ALL_TARGET_FILES


rule mapping:
    input:
        expand(SORTED_BAM, sample=SAMPLES),
        expand(SORTED_BAI, sample=SAMPLES)


rule reference:
    input:
        REF_TARGET_FILES


rule baseline_qc:
    input:
        BASELINE_QC_TARGET_FILES


rule reference_qc:
    input:
        REFERENCE_QC_TARGET_FILES


rule cnv_qc:
    input:
        (
            expand(CNV_QC_TSV, sample=SAMPLES)
            + expand(CNV_QC_PLOT, sample=SAMPLES)
        ) if CNV_ENABLED else []


rule cnv:
    input:
        (
            (expand(CNV_DONE, sample=SAMPLES) if CNV_POSTPROCESS_PRESERVE_BRANCH_A else [])
            + (expand(CNV_B_FINAL_EVENTS, sample=SAMPLES) if CNV_POSTPROCESS_ENABLE_BRANCH_B else [])
            + (expand(CNV_B_ARTIFACT_SUMMARY, sample=SAMPLES) if CNV_POSTPROCESS_ENABLE_BRANCH_B else [])
            + (expand(CNV_B_FINAL_JSON, sample=SAMPLES) if CNV_POSTPROCESS_ENABLE_BRANCH_B else [])
        ) if CNV_ENABLED else []


rule cnv_eval:
    input:
        [CNV_EVAL_SAMPLE_METRICS, CNV_EVAL_EVENT_METRICS, CNV_EVAL_CALIBRATION, CNV_EVAL_SUMMARY] if CNV_ENABLED else []


rule cnv_ml:
    input:
        [CNV_ML_FEATURES, CNV_ML_CV_METRICS, CNV_ML_CALIBRATION, CNV_ML_IMPORTANCE, CNV_ML_PREDICTIONS, CNV_ML_SUMMARY] if CNV_ENABLED else []


rule cnv_benchmark:
    input:
        [CNV_BENCHMARK_SIMULATION, CNV_BENCHMARK_ADMIXTURE, CNV_BENCHMARK_SUMMARY] if CNV_ENABLED else []


rule cnv_report:
    input:
        [CNV_REPORT_TSV, CNV_REPORT_JSON, CNV_REPORT_MD, CNV_REPORT_HTML] if CNV_ENABLED else []


include: "rules/common_preprocess.smk"
include: "rules/qc_workflow.smk"
include: "rules/reference_workflow.smk"
include: "rules/predict_workflow.smk"
include: "rules/runtime_tracking.smk"
include: "rules/reference_assets.smk"
include: "rules/qc_framework.smk"
