CNV_ENABLED = bool(CNV_CFG.get("enable", False))
CNV_DIR = resolve_path(CNV_CFG.get("output_dir", "wisecondorx/cnv"))
CNV_QC_DIR = resolve_path(CNV_CFG.get("qc_dir", "wisecondorx/cnv/qc"))
CNV_PREDICT_DIR = resolve_path(CNV_CFG.get("predict_dir", "wisecondorx/cnv/predict"))
CNV_GENDER_DIR = resolve_path(CNV_CFG.get("gender_dir", "wisecondorx/cnv/gender"))
CNV_CONVERT_BINSIZE = int(CNV_CFG.get("convert_binsize", 100000))
CNV_ZSCORE = float(CNV_CFG.get("zscore", 8))
CNV_ALPHA = float(CNV_CFG.get("alpha", 0.001))
CNV_MASKREPEATS = int(CNV_CFG.get("maskrepeats", 5))
CNV_MINREFBINS = int(CNV_CFG.get("minrefbins", 150))
CNV_PREDICT_SEED = int(CNV_CFG.get("seed", 1))
CNV_SEX_CALL_METHOD = str(CNV_SEX_CFG.get("method", "wisecondorx_plus_bam_depth"))
CNV_SEX_XX_MIN_X_REL = float(CNV_SEX_CFG.get("bam_xx_min_x_relative", 0.95))
CNV_SEX_XY_MAX_X_REL = float(CNV_SEX_CFG.get("bam_xy_max_x_relative", 0.80))
CNV_SEX_XY_MIN_Y_REL = float(CNV_SEX_CFG.get("bam_xy_min_y_relative", 0.20))
CNV_SEX_XX_MAX_Y_REL = float(CNV_SEX_CFG.get("bam_xx_max_y_relative", 0.15))
CNV_QC_MIN_TOTAL = float(CNV_QC_CFG.get("min_total_counts", 1000000))
CNV_QC_MIN_NONZERO = float(CNV_QC_CFG.get("min_nonzero_fraction", 0.4))
CNV_QC_MAX_MAD = float(CNV_QC_CFG.get("max_mad_log1p", 1.2))
CNV_POSTPROCESS_ENABLE_BRANCH_B = bool(CNV_POSTPROCESS_CFG.get("enable_branch_b", True))
CNV_POSTPROCESS_PRESERVE_BRANCH_A = bool(CNV_POSTPROCESS_CFG.get("preserve_branch_a", True))
CNV_POSTPROCESS_CORRECTION_MODEL = str(CNV_POSTPROCESS_CFG.get("default_correction", "2d_loess_gc_mappability"))
CNV_POSTPROCESS_GENOME_BUILD = str(CNV_POSTPROCESS_ANNOTATION_CFG.get("genome_build", "hg19"))
CNV_CORRECTION_LOESS_FRAC = float(CNV_CORRECTION_CFG.get("loess_frac", 0.2))
CNV_CORRECTION_MIN_VALID_BINS = int(CNV_CORRECTION_CFG.get("min_valid_bins", 200))
CNV_CORRECTION_ROBUST_ITERS = int(CNV_CORRECTION_CFG.get("robust_iters", 2))
CNV_CORRECTION_INCLUDE_MASK_LABELS = [str(item) for item in CNV_CORRECTION_CFG.get("include_mask_labels", ["pass", "soft"])]
CNV_POSTPROCESS_PAR_REGIONS = [
    str(item)
    for item in CNV_POSTPROCESS_ANNOTATION_CFG.get(
        "par_regions",
        [
            "chrX:60001-2699520",
            "chrX:154931044-155260560",
        ],
    )
]
CNV_POSTPROCESS_XTR_REGIONS = [str(item) for item in CNV_POSTPROCESS_ANNOTATION_CFG.get("xtr_regions", [])]
CNV_POSTPROCESS_SEGMENTAL_DUPLICATION_BED = (
    resolve_path(CNV_POSTPROCESS_ANNOTATION_CFG.get("segmental_duplication_bed", ""))
    if CNV_POSTPROCESS_ANNOTATION_CFG.get("segmental_duplication_bed", "")
    else ""
)
CNV_POSTPROCESS_LOW_MAPPABILITY_BED = (
    resolve_path(CNV_POSTPROCESS_ANNOTATION_CFG.get("low_mappability_bed", ""))
    if CNV_POSTPROCESS_ANNOTATION_CFG.get("low_mappability_bed", "")
    else ""
)
CNV_POSTPROCESS_GAP_BED = (
    resolve_path(CNV_POSTPROCESS_ANNOTATION_CFG.get("gap_centromere_telomere_bed", ""))
    if CNV_POSTPROCESS_ANNOTATION_CFG.get("gap_centromere_telomere_bed", "")
    else ""
)
CNV_POSTPROCESS_REPEAT_BED = (
    resolve_path(CNV_POSTPROCESS_ANNOTATION_CFG.get("repeat_rich_bed", ""))
    if CNV_POSTPROCESS_ANNOTATION_CFG.get("repeat_rich_bed", "")
    else ""
)
CNV_POSTPROCESS_BLACKLIST_BED = (
    resolve_path(CNV_POSTPROCESS_ANNOTATION_CFG.get("blacklist_bed", ""))
    if CNV_POSTPROCESS_ANNOTATION_CFG.get("blacklist_bed", "")
    else ""
)
CNV_POSTPROCESS_SEX_HOMOLOGY_BED = (
    resolve_path(CNV_POSTPROCESS_ANNOTATION_CFG.get("sex_homology_bed", ""))
    if CNV_POSTPROCESS_ANNOTATION_CFG.get("sex_homology_bed", "")
    else ""
)
CNV_POSTPROCESS_AMBIGUOUS_ALIGNMENT_BED = (
    resolve_path(CNV_POSTPROCESS_ANNOTATION_CFG.get("ambiguous_alignment_bed", ""))
    if CNV_POSTPROCESS_ANNOTATION_CFG.get("ambiguous_alignment_bed", "")
    else ""
)
CNV_POSTPROCESS_DIR = str(Path(CNV_DIR) / "postprocess")
CNV_B_CORRECTED_BINS = str(Path(CNV_POSTPROCESS_DIR) / "correction" / "{sample}.bins.tsv")
CNV_B_CORRECTION_SUMMARY = str(Path(CNV_POSTPROCESS_DIR) / "correction" / "{sample}.summary.json")
CNV_B_BINS = str(Path(CNV_POSTPROCESS_DIR) / "calling" / "{sample}.bins.tsv")
CNV_B_CANDIDATES = str(Path(CNV_POSTPROCESS_DIR) / "calling" / "{sample}.candidate_events.tsv")
CNV_B_CALLING_SUMMARY = str(Path(CNV_POSTPROCESS_DIR) / "calling" / "{sample}.summary.json")
CNV_B_CALIBRATED_BINS = str(Path(CNV_POSTPROCESS_DIR) / "calibration" / "{sample}.bins.tsv")
CNV_B_CALIBRATED_CANDIDATES = str(Path(CNV_POSTPROCESS_DIR) / "calibration" / "{sample}.candidate_events.tsv")
CNV_B_CALIBRATION_SUMMARY = str(Path(CNV_POSTPROCESS_DIR) / "calibration" / "{sample}.summary.json")
CNV_B_MOSAIC_CANDIDATES = str(Path(CNV_POSTPROCESS_DIR) / "mosaic_fraction" / "{sample}.candidate_events.tsv")
CNV_B_MOSAIC_SUMMARY = str(Path(CNV_POSTPROCESS_DIR) / "mosaic_fraction" / "{sample}.summary.json")
CNV_B_FINAL_EVENTS = str(Path(CNV_POSTPROCESS_DIR) / "artifact_rules" / "{sample}.candidate_events.tsv")
CNV_B_ARTIFACT_SUMMARY = str(Path(CNV_POSTPROCESS_DIR) / "artifact_rules" / "{sample}.summary.tsv")
CNV_B_FINAL_JSON = str(Path(CNV_POSTPROCESS_DIR) / "artifact_rules" / "{sample}.candidate_events.json")
CNV_CALLING_MIN_BINS = int(CNV_CALLING_CFG.get("min_bins", 5))
CNV_CALLING_MAX_SEGMENTS = int(CNV_CALLING_CFG.get("max_segments_per_chrom", 12))
CNV_CALLING_SPLIT_THRESHOLD = float(CNV_CALLING_CFG.get("split_threshold", 2.5))
CNV_CALLING_HMM_SHIFT = float(CNV_CALLING_CFG.get("hmm_state_shift", 2.5))
CNV_CALLING_HMM_STAY_PROB = float(CNV_CALLING_CFG.get("hmm_stay_prob", 0.995))
CNV_CALLING_MIN_EVENT_BINS = int(CNV_CALLING_CFG.get("min_event_bins", 3))
CNV_CALLING_MIN_EVENT_Z = float(CNV_CALLING_CFG.get("min_event_z", 1.5))
CNV_CAL_NULL_LOW = float(CNV_CALIBRATION_CFG.get("null_quantile_low", 0.1))
CNV_CAL_NULL_HIGH = float(CNV_CALIBRATION_CFG.get("null_quantile_high", 0.9))
CNV_CAL_MIN_NULL_BINS = int(CNV_CALIBRATION_CFG.get("min_null_bins", 200))
CNV_CAL_EVENT_Z_THRESHOLD = float(CNV_CALIBRATION_CFG.get("event_z_threshold", 2.5))
CNV_MOSAIC_MIN_EFFECTIVE_BINS = float(CNV_MOSAIC_CFG.get("min_effective_bins", 5))
CNV_MOSAIC_MIN_CLEAN_FRACTION = float(CNV_MOSAIC_CFG.get("min_clean_fraction", 0.5))
CNV_MOSAIC_MAX_HIGH_RISK_FRACTION = float(CNV_MOSAIC_CFG.get("max_high_risk_fraction", 0.25))
CNV_MOSAIC_MIN_ABS_LOG2_RATIO = float(CNV_MOSAIC_CFG.get("min_abs_log2_ratio", 0.03))
CNV_MOSAIC_LOW_FRACTION_THRESHOLD = float(CNV_MOSAIC_CFG.get("low_fraction_threshold", 0.15))
CNV_MOSAIC_BASELINE_MIN_BINS = int(CNV_MOSAIC_CFG.get("baseline_min_bins", 200))
CNV_MOSAIC_CI_ZSCORE = float(CNV_MOSAIC_CFG.get("ci_zscore", 1.96))
CNV_ARTIFACT_MIN_BINS = int(CNV_ARTIFACT_CFG.get("min_event_bins", 3))
CNV_ARTIFACT_MIN_ABS_Z = float(CNV_ARTIFACT_CFG.get("min_abs_calibrated_z", 2.0))
CNV_ARTIFACT_MAX_CHROM_FRAC = float(CNV_ARTIFACT_CFG.get("max_chrom_fraction", 0.35))
CNV_ARTIFACT_EDGE_WINDOW = int(CNV_ARTIFACT_CFG.get("edge_bin_window", 2))
CNV_ARTIFACT_MAX_QVALUE = float(CNV_ARTIFACT_CFG.get("max_qvalue", 0.25))
CNV_ARTIFACT_KEEP_REVIEW = int(bool(CNV_ARTIFACT_CFG.get("keep_review", True)))
CNV_ARTIFACT_HIGH_CONF_Z = float(CNV_ARTIFACT_CFG.get("high_confidence_z", 4.0))
CNV_ARTIFACT_HIGH_CONF_QVALUE = float(CNV_ARTIFACT_CFG.get("high_confidence_qvalue", 0.05))
CNV_ML_BACKEND = str(CNV_ML_CFG.get("backend", "auto"))
CNV_ML_CV_FOLDS = int(CNV_ML_CFG.get("cv_folds", 5))
CNV_ML_LABELS_TSV = resolve_path(CNV_ML_CFG.get("labels_tsv", "")) if str(CNV_ML_CFG.get("labels_tsv", "")).strip() else ""
CNV_EVAL_TRUTH_TSV = resolve_path(CNV_EVALUATION_CFG.get("truth_tsv", "")) if str(CNV_EVALUATION_CFG.get("truth_tsv", "")).strip() else ""
CNV_BENCHMARK_TRUTH_TSV = (
    resolve_path(CNV_BENCHMARK_CFG.get("truth_tsv", ""))
    if str(CNV_BENCHMARK_CFG.get("truth_tsv", "")).strip()
    else CNV_EVAL_TRUTH_TSV
)
CNV_MOSAIC_FRACTION_TRUTH_TSV = (
    resolve_path(CNV_BENCHMARK_CFG.get("mosaic_fraction_truth_tsv", ""))
    if str(CNV_BENCHMARK_CFG.get("mosaic_fraction_truth_tsv", "")).strip()
    else ""
)
CNV_BENCHMARK_ADMIXTURE_LEVELS = [float(item) for item in CNV_BENCHMARK_CFG.get("admixture_levels", [1.0, 0.75, 0.5, 0.3, 0.2, 0.1])]
CNV_BENCHMARK_LOW_FRACTION_THRESHOLDS = [
    float(item) for item in CNV_BENCHMARK_CFG.get("low_fraction_thresholds", [0.05, 0.10, 0.15, 0.20, 0.30])
]
CNV_BENCHMARK_BRANCH_B_Z_THRESHOLD = float(CNV_BENCHMARK_CFG.get("branch_b_z_threshold", CNV_CAL_EVENT_Z_THRESHOLD))
CNV_BENCHMARK_BRANCH_A_Z_THRESHOLD = float(CNV_BENCHMARK_CFG.get("branch_a_z_threshold", CNV_ZSCORE))
CNV_EVAL_DIR = str(Path(CNV_DIR) / "evaluation")
CNV_ML_DIR = str(Path(CNV_DIR) / "ml")
CNV_BENCHMARK_DIR = str(Path(CNV_DIR) / "benchmarking")
CNV_REPORT_DIR = str(Path(CNV_DIR) / "report")
CNV_EVAL_SAMPLE_METRICS = str(Path(CNV_EVAL_DIR) / "sample_metrics.tsv")
CNV_EVAL_EVENT_METRICS = str(Path(CNV_EVAL_DIR) / "event_metrics.tsv")
CNV_EVAL_CALIBRATION = str(Path(CNV_EVAL_DIR) / "calibration.tsv")
CNV_EVAL_SUMMARY = str(Path(CNV_EVAL_DIR) / "summary.json")
CNV_ML_FEATURES = str(Path(CNV_ML_DIR) / "candidate_event_features.tsv")
CNV_ML_CV_METRICS = str(Path(CNV_ML_DIR) / "cv_metrics.tsv")
CNV_ML_CALIBRATION = str(Path(CNV_ML_DIR) / "calibration.tsv")
CNV_ML_IMPORTANCE = str(Path(CNV_ML_DIR) / "feature_importance.tsv")
CNV_ML_PREDICTIONS = str(Path(CNV_ML_DIR) / "predictions.tsv")
CNV_ML_SUMMARY = str(Path(CNV_ML_DIR) / "summary.json")
CNV_BENCHMARK_SIMULATION = str(Path(CNV_BENCHMARK_DIR) / "simulation_truth_hits.tsv")
CNV_BENCHMARK_ADMIXTURE = str(Path(CNV_BENCHMARK_DIR) / "admixture_sensitivity.tsv")
CNV_BENCHMARK_SUMMARY = str(Path(CNV_BENCHMARK_DIR) / "summary.json")
CNV_MOSAIC_TRUTH_VALIDATION_SUMMARY = str(Path(CNV_BENCHMARK_DIR) / "mosaic_truth_validation.json")
CNV_REPORT_TSV = str(Path(CNV_REPORT_DIR) / "cnv_summary.tsv")
CNV_REPORT_JSON = str(Path(CNV_REPORT_DIR) / "cnv_summary.json")
CNV_REPORT_MD = str(Path(CNV_REPORT_DIR) / "cnv_summary.md")
CNV_REPORT_HTML = str(Path(CNV_REPORT_DIR) / "cnv_summary.html")
CNV_NPZ = str(Path(CNV_DIR) / "npz" / "{sample}.npz")
CNV_GENDER_TSV = str(Path(CNV_GENDER_DIR) / "{sample}.gender.tsv")
CNV_QC_TSV = str(Path(CNV_QC_DIR) / "{sample}.qc.tsv")
CNV_QC_PLOT = str(Path(CNV_QC_DIR) / "{sample}.qc.svg")
CNV_QC_PASS = str(Path(CNV_QC_DIR) / "{sample}.pass")
CNV_DONE = str(Path(CNV_PREDICT_DIR) / "{sample}.done")
CNV_A_ABERRATIONS_BED = str(Path(CNV_PREDICT_DIR) / "{sample}_aberrations.bed")
