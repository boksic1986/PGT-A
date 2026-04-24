REF_OUTPUT_CFG = WISE_CFG.get("reference_output", "")
REF_OUTPUT = resolve_path(REF_OUTPUT_CFG) if REF_OUTPUT_CFG else ""
REF_OUTPUT_BY_SEX_CFG = WISE_CFG.get("reference_output_by_sex", {})
SEX_SPECIFIC_REF_CONFIGURED = any(REF_OUTPUT_BY_SEX_CFG.get(sex) for sex in ("XX", "XY"))
if SEX_SPECIFIC_REF_CONFIGURED:
    missing_ref_outputs = [sex for sex in ("XX", "XY") if not REF_OUTPUT_BY_SEX_CFG.get(sex)]
    if missing_ref_outputs:
        raise ValueError(
            "core.wisecondorx.reference_output_by_sex must define both XX and XY when sex-specific references are enabled: "
            + ",".join(missing_ref_outputs)
        )
REF_MODEL_ROOT = resolve_path(WISE_CFG.get("reference_model_root", "reference"))
GENDER_REF_OUTPUT = resolve_path(WISE_CFG.get("gender_reference_output", str(Path(REF_MODEL_ROOT) / "gender" / "result" / "ref_gender_best.npz")))
COMMON_REF_BINSIZE = resolve_path(WISE_CFG.get("common_reference_binsize_output", str(Path(REF_MODEL_ROOT) / "gender" / "common_best_binsize.txt")))
REF_GROUPS = BUILD_REF_CFG.get("groups", {})
AVAILABLE_REF_SEXES = [sex for sex in ("XX", "XY") if REF_OUTPUT_BY_SEX_CFG.get(sex)]
REF_OUTPUTS_BY_SEX = {sex: resolve_path(REF_OUTPUT_BY_SEX_CFG[sex]) for sex in AVAILABLE_REF_SEXES}
REF_GROUPS_BY_SEX_AVAILABLE = all(REF_GROUPS.get(sex) for sex in ("XX", "XY"))
BUILD_REF_BY_SEX_ENABLED = SEX_SPECIFIC_REF_CONFIGURED and REF_GROUPS_BY_SEX_AVAILABLE
PREDICT_BY_SEX_ENABLED = SEX_SPECIFIC_REF_CONFIGURED
REF_SEXES = [sex for sex in AVAILABLE_REF_SEXES if REF_GROUPS.get(sex)] if BUILD_REF_BY_SEX_ENABLED else []
REF_SAMPLE_IDS_BY_SEX = {}
for sex in REF_SEXES:
    ids = [sample_id for sample_id in REF_GROUPS.get(sex, []) if sample_id in SAMPLES]
    if not ids:
        raise ValueError(f"build_reference.groups[{sex}] has no valid samples in config['samples']")
    REF_SAMPLE_IDS_BY_SEX[sex] = sorted(ids)

REF_SAMPLE_IDS = []
if REF_SEXES:
    seen_samples = set()
    for sex in REF_SEXES:
        for sample_id in REF_SAMPLE_IDS_BY_SEX[sex]:
            if sample_id not in seen_samples:
                REF_SAMPLE_IDS.append(sample_id)
                seen_samples.add(sample_id)
elif REF_GROUPS:
    seen_samples = set()
    for group_name in sorted(REF_GROUPS):
        for sample_id in REF_GROUPS.get(group_name, []):
            if sample_id in SAMPLES and sample_id not in seen_samples:
                REF_SAMPLE_IDS.append(sample_id)
                seen_samples.add(sample_id)
    if not REF_SAMPLE_IDS:
        raise ValueError("build_reference.groups does not contain any valid samples present in config['samples'].")
else:
    REF_SAMPLE_IDS = list(SAMPLES)
REF_SAMPLE_IDS = sorted(REF_SAMPLE_IDS)
BASELINE_SAMPLE_IDS = list(REF_SAMPLE_IDS)

REF_PREFILTER_DIR = str(Path(REF_MODEL_ROOT) / "prefilter")
REF_TUNING_DIR = str(Path(REF_MODEL_ROOT) / "tuning")
REF_PREFILTER_DIR_BY_SEX = {sex: str(Path(REF_MODEL_ROOT) / sex / "prefilter") for sex in REF_SEXES}
REF_TUNING_DIR_BY_SEX = {sex: str(Path(REF_MODEL_ROOT) / sex / "tuning") for sex in REF_SEXES}
REF_PREFILTER_MERGED_INLIERS = str(Path(REF_PREFILTER_DIR) / "reference_inlier_samples.txt")
BASELINE_QC_DIR = project_path("qc", "baseline")
BASELINE_QC_TSV = project_path("qc", "baseline", "{sample}", "qc_metrics.tsv")
BASELINE_QC_PROFILE_TSV = project_path("qc", "baseline", "{sample}", "target_bin_profile.tsv")
BASELINE_QC_REF_SUMMARY_TSV = project_path("qc", "baseline", "{sample}", "reference_summary.tsv")
BASELINE_QC_PLOT = project_path("qc", "baseline", "{sample}", "target_vs_ref_profile.png")
BASELINE_QC_GC_PLOT = project_path("qc", "baseline", "{sample}", "gc_bias_plot.png")
BASELINE_QC_BIN_SIZE = int(BASELINE_QC_CFG.get("bin_size", WISE_CFG.get("binsize", 200000)))
BASELINE_QC_THREADS = int(BASELINE_QC_CFG.get("threads", 96))
BASELINE_QC_SAMTOOLS_THREADS = int(BASELINE_QC_CFG.get("samtools_threads", 4))
BASELINE_QC_BIN_SIZES = sorted({int(item) for item in BASELINE_QC_CFG.get("bin_sizes", [BASELINE_QC_BIN_SIZE])})
BASELINE_QC_BIN_LABELS = [format_binsize_label(item) for item in BASELINE_QC_BIN_SIZES]
BASELINE_QC_GC_CORRECTION_METHOD = str(BASELINE_QC_GC_CORRECTION_CFG.get("method", "loess")).strip().lower()
BASELINE_QC_GC_CORRECTION_FRAC = float(BASELINE_QC_GC_CORRECTION_CFG.get("loess_frac", 0.2))
BASELINE_QC_GC_CORRECTION_POLY_DEGREE = int(BASELINE_QC_GC_CORRECTION_CFG.get("poly_degree", 2))
BASELINE_QC_GC_CORRECTION_MIN_VALID_BINS = int(BASELINE_QC_GC_CORRECTION_CFG.get("min_valid_bins", 200))
BASELINE_QC_GC_CORRECTION_ROBUST_ITERS = int(BASELINE_QC_GC_CORRECTION_CFG.get("robust_iters", 1))
BASELINE_QC_SUMMARY = project_path("qc", "baseline", "baseline_qc_summary.tsv")
BASELINE_QC_PASS_SAMPLES = project_path("qc", "baseline", "baseline_qc_pass_samples.txt")
BASELINE_QC_RETAINED_SAMPLES = project_path("qc", "baseline", "retained_baseline_samples.tsv")
BASELINE_QC_OUTLIER_SAMPLES = project_path("qc", "baseline", "outlier_samples.tsv")
BASELINE_QC_REPORT_MD = project_path("qc", "baseline", "baseline_qc_report.md")
BASELINE_QC_FIGURES_DIR = project_path("qc", "baseline", "figures")
BASELINE_QC_FIG_DECISION = project_path("qc", "baseline", "figures", "qc_decision_counts.png")
BASELINE_QC_FIG_MAPPED = project_path("qc", "baseline", "figures", "mapped_fragments_by_sample.png")
BASELINE_QC_FIG_METRICS = project_path("qc", "baseline", "figures", "metric_scatter.png")
BASELINE_QC_FIG_BATCH = project_path("qc", "baseline", "figures", "batch_metric_boxplots.png")
BASELINE_QC_FIG_CORR = project_path("qc", "baseline", "figures", "sample_correlation_heatmap.png")
BASELINE_MULTI_DIR = project_path("qc", "baseline_multiscale")
BASELINE_MULTI_TSV = project_path("qc", "baseline_multiscale", "{bin_label}", "{sample}", "qc_metrics.tsv")
BASELINE_MULTI_PROFILE_TSV = project_path("qc", "baseline_multiscale", "{bin_label}", "{sample}", "target_bin_profile.tsv")
BASELINE_MULTI_REF_SUMMARY_TSV = project_path("qc", "baseline_multiscale", "{bin_label}", "{sample}", "reference_summary.tsv")
BASELINE_MULTI_PLOT = project_path("qc", "baseline_multiscale", "{bin_label}", "{sample}", "target_vs_ref_profile.png")
BASELINE_MULTI_GC_PLOT = project_path("qc", "baseline_multiscale", "{bin_label}", "{sample}", "gc_bias_plot.png")
BASELINE_MULTI_SUMMARY = project_path("qc", "baseline_multiscale", "{bin_label}", "baseline_qc_summary.tsv")
BASELINE_MULTI_PASS_SAMPLES = project_path("qc", "baseline_multiscale", "{bin_label}", "baseline_qc_pass_samples.txt")
BASELINE_MULTI_RETAINED_SAMPLES = project_path("qc", "baseline_multiscale", "{bin_label}", "retained_baseline_samples.tsv")
BASELINE_MULTI_OUTLIER_SAMPLES = project_path("qc", "baseline_multiscale", "{bin_label}", "outlier_samples.tsv")
BASELINE_MULTI_REPORT_MD = project_path("qc", "baseline_multiscale", "{bin_label}", "baseline_qc_report.md")
BASELINE_MULTI_FIGURES_DIR = project_path("qc", "baseline_multiscale", "{bin_label}", "figures")
BASELINE_MULTI_FIG_DECISION = project_path("qc", "baseline_multiscale", "{bin_label}", "figures", "qc_decision_counts.png")
BASELINE_MULTI_FIG_MAPPED = project_path("qc", "baseline_multiscale", "{bin_label}", "figures", "mapped_fragments_by_sample.png")
BASELINE_MULTI_FIG_METRICS = project_path("qc", "baseline_multiscale", "{bin_label}", "figures", "metric_scatter.png")
BASELINE_MULTI_FIG_BATCH = project_path("qc", "baseline_multiscale", "{bin_label}", "figures", "batch_metric_boxplots.png")
BASELINE_MULTI_FIG_CORR = project_path("qc", "baseline_multiscale", "{bin_label}", "figures", "sample_correlation_heatmap.png")
BASELINE_DIAG_DIR = project_path("qc", "baseline_diagnostics")
BASELINE_FLAGSTAT = project_path("qc", "baseline_diagnostics", "{sample}", "samtools.flagstat.txt")
BASELINE_IDXSTATS = project_path("qc", "baseline_diagnostics", "{sample}", "samtools.idxstats.tsv")
BASELINE_REPORT_DIR = project_path("report")
BASELINE_REPORT_MD = project_path("report", "qc_report.md")
BASELINE_MULTI_BINSIZE_SUMMARY = project_path("qc", "baseline_multiscale", "multiscale_binsize_summary.tsv")
BASELINE_MULTI_SAMPLE_MATRIX = project_path("qc", "baseline_multiscale", "multiscale_sample_decisions.tsv")
BASELINE_MULTI_STABLE_RETAINED = project_path("qc", "baseline_multiscale", "stable_retained_baseline_samples.tsv")
BASELINE_MULTI_FAILED_DIAG = project_path("qc", "baseline_multiscale", "failed_sample_diagnostics.tsv")
BASELINE_MULTI_FIG_DECISION_TREND = project_path("qc", "baseline_multiscale", "figures", "binsize_decision_trend.png")
BASELINE_MULTI_FIG_SAMPLE_HEATMAP = project_path("qc", "baseline_multiscale", "figures", "sample_decision_heatmap.png")
BASELINE_MULTI_FIG_LIBRARY = project_path("qc", "baseline_multiscale", "figures", "library_metric_overview.png")
BASELINE_MULTI_FIG_FAIL_CHROM = project_path("qc", "baseline_multiscale", "figures", "failed_sample_chromosome_bias.png")
REFERENCE_ASSETS_DIR = str(Path(REF_MODEL_ROOT) / "assets")
REFERENCE_BINS_DIR = str(Path(REFERENCE_ASSETS_DIR) / "bins")
REFERENCE_MASK_DIR = str(Path(REFERENCE_ASSETS_DIR) / "mask")
REFERENCE_ATOMIC_BINS = str(Path(REFERENCE_BINS_DIR) / "atomic_bins.tsv")
REFERENCE_ANALYSIS_BINS = str(Path(REFERENCE_BINS_DIR) / "analysis_bins.tsv")
REFERENCE_QC_BINS = str(Path(REFERENCE_BINS_DIR) / "qc_bins.tsv")
REFERENCE_ATOMIC_BIN_ANNOTATIONS = str(Path(REFERENCE_BINS_DIR) / "atomic_bin_annotations.tsv")
REFERENCE_ANALYSIS_BIN_ANNOTATIONS = str(Path(REFERENCE_BINS_DIR) / "analysis_bin_annotations.tsv")
REFERENCE_QC_BIN_ANNOTATIONS = str(Path(REFERENCE_BINS_DIR) / "qc_bin_annotations.tsv")
REFERENCE_BIN_SUMMARY_JSON = str(Path(REFERENCE_BINS_DIR) / "bin_summary.json")
REFERENCE_HARD_MASK_TSV = str(Path(REFERENCE_MASK_DIR) / "hard_mask.tsv")
REFERENCE_SOFT_MASK_TSV = str(Path(REFERENCE_MASK_DIR) / "soft_mask.tsv")
REFERENCE_DYNAMIC_MASK_TSV = str(Path(REFERENCE_MASK_DIR) / "dynamic_mask.tsv")
REFERENCE_COMBINED_MASK_TSV = str(Path(REFERENCE_MASK_DIR) / "combined_mask.tsv")
REFERENCE_HARD_MASK_JSON = str(Path(REFERENCE_MASK_DIR) / "hard_mask.json")
REFERENCE_SOFT_MASK_JSON = str(Path(REFERENCE_MASK_DIR) / "soft_mask.json")
REFERENCE_DYNAMIC_MASK_JSON = str(Path(REFERENCE_MASK_DIR) / "dynamic_mask.json")
REFERENCE_COMBINED_MASK_JSON = str(Path(REFERENCE_MASK_DIR) / "combined_mask.json")
REFERENCE_MASK_SUMMARY_JSON = str(Path(REFERENCE_MASK_DIR) / "mask_summary.json")
REFERENCE_PACKAGE_NAME = str(REFERENCE_PACKAGE_CFG.get("name", "frozen_reference_v1"))
REFERENCE_PACKAGE_DIR = str(Path(REF_MODEL_ROOT) / "package" / REFERENCE_PACKAGE_NAME)
REFERENCE_PACKAGE_MANIFEST = str(Path(REFERENCE_PACKAGE_DIR) / "manifest.json")
REFERENCE_PACKAGE_INVENTORY = str(Path(REFERENCE_PACKAGE_DIR) / "inventory.tsv")
REFERENCE_PACKAGE_README = str(Path(REFERENCE_PACKAGE_DIR) / "README.md")
REFERENCE_PACKAGE_DONE = str(Path(REFERENCE_PACKAGE_DIR) / "package.done")
QC_FRAMEWORK_DIR = project_path("qc", "framework")
RUN_QC_TSV = project_path("qc", "framework", "run_qc.tsv")
RUN_QC_JSON = project_path("qc", "framework", "run_qc.json")
SAMPLE_QC_TSV = project_path("qc", "framework", "sample_qc.tsv")
SAMPLE_QC_JSON = project_path("qc", "framework", "sample_qc.json")
BIN_QC_TSV = project_path("qc", "framework", "bin_qc.tsv")
BIN_QC_JSON = project_path("qc", "framework", "bin_qc.json")
EVENT_QC_TSV = project_path("qc", "framework", "event_qc.tsv")
EVENT_QC_JSON = project_path("qc", "framework", "event_qc.json")
QC_FRAMEWORK_FIG_SAMPLE = project_path("qc", "framework", "figures", "sample_status_counts.png")
QC_FRAMEWORK_FIG_EVENT = project_path("qc", "framework", "figures", "event_status_counts.png")
QC_FRAMEWORK_FIG_MASK = project_path("qc", "framework", "figures", "mask_status_counts.png")
QC_FRAMEWORK_FIG_METRIC = project_path("qc", "framework", "figures", "metric_scatter.png")
PREFILTER_BINSIZE = int(PREFILTER_CFG.get("binsize", 100000))
PREFILTER_MAX_ITER = int(PREFILTER_CFG.get("max_iterations", 3))
ATOMIC_BINSIZE = int(
    REFERENCE_ASSETS_CFG.get(
        "atomic_binsize",
        min(PREFILTER_BINSIZE, int(WISE_CFG.get("binsize", PREFILTER_BINSIZE)), BASELINE_QC_BIN_SIZE),
    )
)
ANALYSIS_BINSIZE_DEFAULT = int(REFERENCE_ASSETS_CFG.get("analysis_binsize", int(WISE_CFG.get("binsize", BASELINE_QC_BIN_SIZE))))
QC_BINSIZE = int(REFERENCE_ASSETS_CFG.get("qc_binsize", BASELINE_QC_BIN_SIZE))
HARD_MASK_N_FRACTION = float(REFERENCE_ASSETS_CFG.get("hard_n_fraction", 0.2))
SOFT_MASK_GC_LOW = float(REFERENCE_ASSETS_CFG.get("soft_gc_low", 0.2))
SOFT_MASK_GC_HIGH = float(REFERENCE_ASSETS_CFG.get("soft_gc_high", 0.8))
DYNAMIC_MASK_Z_FRAC = float(REFERENCE_ASSETS_CFG.get("dynamic_z_frac_threshold", 0.05))
DYNAMIC_MASK_MEDIAN_ABS_Z = float(REFERENCE_ASSETS_CFG.get("dynamic_median_abs_z_threshold", 1.5))

TUNING_ENABLED = bool(TUNE_CFG.get("enable", False))
TUNING_BIN_SIZES = [int(item) for item in TUNE_CFG.get("bin_sizes", [int(WISE_CFG["binsize"])])]
TUNING_WORKDIR = project_path("wisecondorx", "tuning")
TUNING_SUMMARY = project_path("wisecondorx", "tuning", "bin_pca_grid.tsv")
TUNING_BINSIZE_SUMMARY = project_path("wisecondorx", "tuning", "binsize_ranking.tsv")
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

REF_SET_CFG_RAW = BUILD_REF_CFG.get("reference_sets", {})
if not REF_SET_CFG_RAW:
    REF_SET_CFG_RAW = {
        "pass_only": {"decisions": ["PASS"], "publish_as_default": True},
        "pass_warn": {"decisions": ["PASS", "WARN"], "publish_as_default": False},
    }
ACTIVE_REF_SET_NAMES_RAW = BUILD_REF_CFG.get("active_reference_sets", [])
ACTIVE_REF_SET_NAMES = [str(name).strip() for name in ACTIVE_REF_SET_NAMES_RAW if str(name).strip()]
if ACTIVE_REF_SET_NAMES:
    unknown_active_ref_sets = sorted(set(ACTIVE_REF_SET_NAMES) - set(REF_SET_CFG_RAW))
    if unknown_active_ref_sets:
        raise ValueError(
            "build_reference.active_reference_sets contains unknown entries: "
            + ",".join(unknown_active_ref_sets)
        )

REF_SET_ORDER = [name for name in ("pass_only", "pass_warn") if name in REF_SET_CFG_RAW] + [
    name for name in sorted(REF_SET_CFG_RAW) if name not in {"pass_only", "pass_warn"}
]
if ACTIVE_REF_SET_NAMES:
    active_ref_set_lookup = set(ACTIVE_REF_SET_NAMES)
    REF_SET_ORDER = [name for name in REF_SET_ORDER if name in active_ref_set_lookup]
if not REF_SET_ORDER:
    raise ValueError("No active reference_sets configured for build_reference.")
REF_SET_CFG = {}
for ref_set in REF_SET_ORDER:
    raw_cfg = dict(REF_SET_CFG_RAW.get(ref_set) or {})
    decisions = [str(item).strip().upper() for item in raw_cfg.get("decisions", []) if str(item).strip()]
    if not decisions:
        raise ValueError(f"build_reference.reference_sets.{ref_set}.decisions must not be empty.")
    min_ref_override_raw = raw_cfg.get("min_reference_samples_override")
    min_ref_override = None if min_ref_override_raw in (None, "") else int(min_ref_override_raw)
    REF_SET_CFG[ref_set] = {
        "decisions": decisions,
        "publish_as_default": bool(raw_cfg.get("publish_as_default", ref_set == "pass_only")),
        "min_reference_samples_override": min_ref_override,
    }
DEFAULT_REF_SET_CANDIDATES = [name for name, cfg in REF_SET_CFG.items() if cfg["publish_as_default"]]
if len(DEFAULT_REF_SET_CANDIDATES) > 1:
    raise ValueError("Only one build_reference.reference_sets entry may set publish_as_default=true.")
DEFAULT_REF_SET = DEFAULT_REF_SET_CANDIDATES[0] if DEFAULT_REF_SET_CANDIDATES else REF_SET_ORDER[0]
REF_SET_SAMPLE_LISTS = {ref_set: project_path("reference", "cohorts", ref_set, "selected_samples.txt") for ref_set in REF_SET_ORDER}
REF_SET_PREFILTER_DIR_BY_SEX = {
    ref_set: {sex: project_path("reference", "cohorts", ref_set, sex, "prefilter") for sex in REF_SEXES}
    for ref_set in REF_SET_ORDER
}
REF_SET_PREFILTER_MERGED_INLIERS = {
    ref_set: project_path("reference", "cohorts", ref_set, "prefilter", "reference_inlier_samples.txt")
    for ref_set in REF_SET_ORDER
}
REF_SET_TUNING_WORKDIR = {ref_set: project_path("wisecondorx", "tuning", ref_set) for ref_set in REF_SET_ORDER}
REF_SET_TUNING_SUMMARY = {ref_set: project_path("wisecondorx", "tuning", ref_set, "bin_pca_grid.tsv") for ref_set in REF_SET_ORDER}
REF_SET_TUNING_BINSIZE_SUMMARY = {
    ref_set: project_path("wisecondorx", "tuning", ref_set, "binsize_ranking.tsv")
    for ref_set in REF_SET_ORDER
}
REF_SET_TUNING_BEST = {ref_set: project_path("wisecondorx", "tuning", ref_set, "best_params.yaml") for ref_set in REF_SET_ORDER}
REF_SET_TUNING_QC = {ref_set: project_path("wisecondorx", "tuning", ref_set, "reference_sample_qc.tsv") for ref_set in REF_SET_ORDER}
REF_SET_TUNING_PLOT = {
    ref_set: project_path("wisecondorx", "tuning", ref_set, "best_bin_pca_elbow.svg")
    for ref_set in REF_SET_ORDER
}
REF_SET_TUNING_QC_STATS_PLOT = {
    ref_set: project_path("wisecondorx", "tuning", ref_set, "reference_qc_metrics.svg")
    for ref_set in REF_SET_ORDER
}
REF_SET_TUNING_INLIERS = {
    ref_set: project_path("wisecondorx", "tuning", ref_set, "reference_inlier_samples.txt")
    for ref_set in REF_SET_ORDER
}
REF_SET_OUTPUTS_BY_SEX = {
    ref_set: {
        sex: project_path("reference", "cohorts", ref_set, sex, "result", f"ref_{sex}_best.npz")
        for sex in REF_SEXES
    }
    for ref_set in REF_SET_ORDER
}
REF_SET_GENDER_OUTPUT = {
    ref_set: project_path("reference", "cohorts", ref_set, "gender", "result", "ref_gender_best.npz")
    for ref_set in REF_SET_ORDER
}
REF_SET_COMMON_BINSIZE = {
    ref_set: project_path("reference", "cohorts", ref_set, "gender", "common_best_binsize.txt")
    for ref_set in REF_SET_ORDER
}
if BUILD_REF_BY_SEX_ENABLED and len(REF_SET_ORDER) > 1 and not TUNING_ENABLED:
    raise ValueError("Multiple reference_sets currently require core.wisecondorx.tuning.enable=true.")

REF_TARGET_FILES = []
if BUILD_REF_BY_SEX_ENABLED:
    if TUNING_ENABLED:
        REF_TARGET_FILES = [
            REF_SET_OUTPUTS_BY_SEX[ref_set][sex]
            for ref_set in REF_SET_ORDER
            for sex in REF_SEXES
        ] + [
            REF_SET_GENDER_OUTPUT[ref_set]
            for ref_set in REF_SET_ORDER
        ] + [
            REF_SET_COMMON_BINSIZE[ref_set]
            for ref_set in REF_SET_ORDER
        ] + [
            REF_OUTPUTS_BY_SEX[sex]
            for sex in REF_SEXES
        ] + [
            GENDER_REF_OUTPUT,
            COMMON_REF_BINSIZE,
        ]
    else:
        REF_TARGET_FILES = [REF_OUTPUTS_BY_SEX[sex] for sex in REF_SEXES] + [GENDER_REF_OUTPUT, COMMON_REF_BINSIZE]
elif REF_OUTPUT:
    REF_TARGET_FILES = [REF_OUTPUT]
