configfile: "config.yaml"

from pathlib import Path


SAMPLES = sorted(config["samples"].keys())
assert SAMPLES, "No samples found in config['samples']."
PIPELINE_ROOT = Path(workflow.basedir).resolve()
SCRIPT_COLLECT_RUN_METADATA = str(PIPELINE_ROOT / "scripts" / "collect_run_metadata.py")
SCRIPT_REFERENCE_PREFILTER_QC = str(PIPELINE_ROOT / "scripts" / "reference_prefilter_qc.py")
SCRIPT_TUNE_WISECONDORX = str(PIPELINE_ROOT / "scripts" / "tune_wisecondorx_bin_pca.py")
SCRIPT_BUILD_REF_FROM_TUNING = str(PIPELINE_ROOT / "scripts" / "build_reference_from_tuning.py")
SCRIPT_CNV_QC = str(PIPELINE_ROOT / "scripts" / "cnv_qc.py")

PROJECT = Path(config["core"]["project_path"])
WISE_CFG = config["core"]["wisecondorx"]
TUNE_CFG = WISE_CFG.get("tuning", {})
TUNE_PCA_CFG = TUNE_CFG.get("pca", {})
TUNE_QC_CFG = TUNE_CFG.get("qc", {})
PREFILTER_CFG = WISE_CFG.get("reference_prefilter", {})
CNV_CFG = WISE_CFG.get("cnv", {})
CNV_QC_CFG = CNV_CFG.get("qc", {})
PIPELINE_CFG = config.get("pipeline", {})
BUILD_REF_CFG = config.get("build_reference", {})


project_path = lambda *parts: str(PROJECT.joinpath(*parts))
resolve_path = lambda path_value: str(Path(path_value)) if Path(path_value).is_absolute() else str(PROJECT / Path(path_value))


FASTP_R1 = project_path("fastp", "{sample}.R1.clean.fastq.gz")
FASTP_R2 = project_path("fastp", "{sample}.R2.clean.fastq.gz")
FASTP_HTML = project_path("fastp", "{sample}.fastp.html")
FASTP_JSON = project_path("fastp", "{sample}.fastp.json")
SORTED_BAM = project_path("mapping", "{sample}.sorted.bam")
SORTED_BAI = project_path("mapping", "{sample}.sorted.bam.bai")
REF_OUTPUT = resolve_path(WISE_CFG["reference_output"])
REF_OUTPUT_BY_SEX_CFG = WISE_CFG.get("reference_output_by_sex", {})
REF_MODEL_ROOT = resolve_path(WISE_CFG.get("reference_model_root", "reference"))
REF_GROUPS = BUILD_REF_CFG.get("groups", {})
REF_SEXES = [
    sex
    for sex in ("XX", "XY")
    if REF_OUTPUT_BY_SEX_CFG.get(sex) and REF_GROUPS.get(sex)
]
REF_OUTPUTS_BY_SEX = {
    sex: resolve_path(REF_OUTPUT_BY_SEX_CFG[sex])
    for sex in REF_SEXES
}
REF_PREFILTER_DIR_BY_SEX = {sex: str(Path(REF_MODEL_ROOT) / sex / "prefilter") for sex in REF_SEXES}
REF_TUNING_DIR_BY_SEX = {sex: str(Path(REF_MODEL_ROOT) / sex / "tuning") for sex in REF_SEXES}
REF_TARGET_FILES = [REF_OUTPUTS_BY_SEX[sex] for sex in REF_SEXES] if REF_SEXES else [REF_OUTPUT]
REF_SAMPLE_IDS_BY_SEX = {}
for sex in REF_SEXES:
    ids = [sample_id for sample_id in REF_GROUPS.get(sex, []) if sample_id in SAMPLES]
    if not ids:
        raise ValueError(f"build_reference.groups[{sex}] has no valid samples in config['samples']")
    REF_SAMPLE_IDS_BY_SEX[sex] = sorted(ids)
TUNING_ENABLED = bool(TUNE_CFG.get("enable", False))
TUNING_BIN_SIZES = [int(item) for item in TUNE_CFG.get("bin_sizes", [int(WISE_CFG["binsize"])])]
TUNING_WORKDIR = project_path("wisecondorx", "tuning")
TUNING_SUMMARY = project_path("wisecondorx", "tuning", "bin_pca_grid.tsv")
TUNING_BEST = project_path("wisecondorx", "tuning", "best_params.yaml")
TUNING_QC = project_path("wisecondorx", "tuning", "reference_sample_qc.tsv")
TUNING_PLOT = project_path("wisecondorx", "tuning", "best_bin_pca_elbow.svg")
TUNING_QC_STATS_PLOT = project_path("wisecondorx", "tuning", "reference_qc_metrics.svg")
TUNING_INLIERS = project_path("wisecondorx", "tuning", "reference_inlier_samples.txt")
TUNING_PCA_MIN = int(TUNE_PCA_CFG.get("min_components", 2))
TUNING_PCA_MAX = int(TUNE_PCA_CFG.get("max_components", 20))
TUNING_MIN_VAR = float(TUNE_PCA_CFG.get("min_explained_variance", 0.0))
TUNING_MIN_REF_SAMPLES = int(TUNE_QC_CFG.get("min_reference_samples", 8))
TUNING_MAX_OUTLIER_FRAC = float(TUNE_QC_CFG.get("max_outlier_fraction", 0.25))
TUNING_MIN_READS = float(TUNE_QC_CFG.get("min_reads_per_sample", 3000000))
TUNING_MIN_CORR = float(TUNE_QC_CFG.get("min_corr_to_median", 0.9))
TUNING_MAX_RECON_Z = float(TUNE_QC_CFG.get("max_reconstruction_error_z", 3.5))
TUNING_MAX_NOISE_Z = float(TUNE_QC_CFG.get("max_noise_mad_z", 3.5))
PREFILTER_BINSIZE = int(PREFILTER_CFG.get("binsize", 100000))
PREFILTER_MAX_ITER = int(PREFILTER_CFG.get("max_iterations", 3))
CNV_ENABLED = bool(CNV_CFG.get("enable", False))
CNV_DIR = resolve_path(CNV_CFG.get("output_dir", "wisecondorx/cnv"))
CNV_QC_DIR = resolve_path(CNV_CFG.get("qc_dir", "wisecondorx/cnv/qc"))
CNV_PREDICT_DIR = resolve_path(CNV_CFG.get("predict_dir", "wisecondorx/cnv/predict"))
CNV_CONVERT_BINSIZE = int(CNV_CFG.get("convert_binsize", 100000))
CNV_ZSCORE = float(CNV_CFG.get("zscore", 8))
CNV_ALPHA = float(CNV_CFG.get("alpha", 0.001))
CNV_MASKREPEATS = int(CNV_CFG.get("maskrepeats", 5))
CNV_MINREFBINS = int(CNV_CFG.get("minrefbins", 150))
CNV_QC_MIN_TOTAL = float(CNV_QC_CFG.get("min_total_counts", 1000000))
CNV_QC_MIN_NONZERO = float(CNV_QC_CFG.get("min_nonzero_fraction", 0.4))
CNV_QC_MAX_MAD = float(CNV_QC_CFG.get("max_mad_log1p", 1.2))
CNV_NPZ = str(Path(CNV_DIR) / "npz" / "{sample}.npz")
CNV_QC_TSV = str(Path(CNV_QC_DIR) / "{sample}.qc.tsv")
CNV_QC_PLOT = str(Path(CNV_QC_DIR) / "{sample}.qc.svg")
CNV_QC_PASS = str(Path(CNV_QC_DIR) / "{sample}.pass")
CNV_DONE = str(Path(CNV_PREDICT_DIR) / "{sample}.done")
RUN_METADATA = project_path("logs", "run_metadata.tsv")
PIPELINE_TARGETS = set(PIPELINE_CFG.get("targets", ["mapping", "reference", "cnv"]))

include: "rules/pipeline.smk"
