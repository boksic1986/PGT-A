CNV_ENABLED = bool(CNV_CFG.get("enable", False))
CNV_DIR = resolve_path(CNV_CFG.get("output_dir", "wisecondorx/cnv"))
CNV_PREPROCESS_CFG = CNV_CFG.get("preprocess", {})
CNV_PREPROCESS_STRATEGY = str(CNV_PREPROCESS_CFG.get("strategy", "mask_only"))
CNV_QC_DIR = resolve_path(CNV_CFG.get("qc_dir", "wisecondorx/cnv/qc"))
CNV_PREDICT_DIR = resolve_path(CNV_CFG.get("predict_dir", "wisecondorx/cnv/predict"))
CNV_GENDER_DIR = resolve_path(CNV_CFG.get("gender_dir", "wisecondorx/cnv/gender"))
CNV_CONVERT_BINSIZE = int(CNV_CFG.get("convert_binsize", 100000))
CNV_ZSCORE = float(CNV_CFG.get("zscore", 8))
CNV_ALPHA = float(CNV_CFG.get("alpha", 0.001))
CNV_MASKREPEATS = int(CNV_CFG.get("maskrepeats", 5))
CNV_MINREFBINS = int(CNV_CFG.get("minrefbins", 150))
CNV_PREDICT_SEED = int(CNV_CFG.get("seed", 1))
CNV_BRANCH_A_CFG = CNV_CFG.get("branch_a", {})
CNV_BRANCH_A_MERGE_GAP_BP = int(CNV_BRANCH_A_CFG.get("merge_gap_bp", 0))
CNV_BRANCH_A_STRONG_Z = float(CNV_BRANCH_A_CFG.get("strong_z", 10.0))
CNV_BRANCH_A_OUTPUT_DIR = resolve_path(
    CNV_BRANCH_A_CFG.get("output_dir", str(Path(CNV_DIR) / "a_branch"))
)
CNV_BRANCH_A_VALIDATION_DIR = resolve_path(
    CNV_BRANCH_A_CFG.get("validation_dir", str(Path(CNV_DIR) / "branch_a_validation"))
)
CNV_BRANCH_A_LOG_DIR = resolve_path(
    CNV_BRANCH_A_CFG.get("log_dir", str(Path("logs") / "cnv"))
)
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
CNV_POSTPROCESS_DIR = resolve_path(
    CNV_POSTPROCESS_CFG.get("output_dir", str(Path(CNV_DIR) / "postprocess"))
)
CNV_PLOTS_DIR = str(Path(CNV_DIR) / "plots")
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
CNV_B_EVIDENCE_LEDGER = str(Path(CNV_POSTPROCESS_DIR) / "evidence_ledger" / "{sample}.candidate_evidence.tsv")
CNV_B_EVIDENCE_SUMMARY = str(Path(CNV_POSTPROCESS_DIR) / "evidence_ledger" / "{sample}.summary.json")
CNV_B_MATCHED_NEGATIVE = str(Path(CNV_POSTPROCESS_DIR) / "matched_negative" / "{sample}.candidate_evidence.tsv")
CNV_B_MATCHED_NEGATIVE_SUMMARY = str(Path(CNV_POSTPROCESS_DIR) / "matched_negative" / "{sample}.summary.json")
CNV_LOWRES_EVIDENCE_CFG = CNV_CFG.get("lowres_evidence", {})
CNV_LOWRES_EVIDENCE_ENABLE = bool(CNV_LOWRES_EVIDENCE_CFG.get("enable", False))
CNV_LOWRES_2MB_EVENTS_DIR = (
    resolve_path(CNV_LOWRES_EVIDENCE_CFG.get("events_2mb_dir", ""))
    if str(CNV_LOWRES_EVIDENCE_CFG.get("events_2mb_dir", "")).strip()
    else ""
)
CNV_LOWRES_3MB_EVENTS_DIR = (
    resolve_path(CNV_LOWRES_EVIDENCE_CFG.get("events_3mb_dir", ""))
    if str(CNV_LOWRES_EVIDENCE_CFG.get("events_3mb_dir", "")).strip()
    else ""
)
CNV_LOWRES_2MB_EVENTS = (
    str(Path(CNV_LOWRES_2MB_EVENTS_DIR) / "{sample}.candidate_events.tsv")
    if CNV_LOWRES_2MB_EVENTS_DIR
    else ""
)
CNV_LOWRES_3MB_EVENTS = (
    str(Path(CNV_LOWRES_3MB_EVENTS_DIR) / "{sample}.candidate_events.tsv")
    if CNV_LOWRES_3MB_EVENTS_DIR
    else ""
)
CNV_LOWRES_REF_NPZ_PATHS = [
    resolve_path(item)
    for item in CNV_LOWRES_EVIDENCE_CFG.get("reference_npz", [])
    if str(item).strip()
]
for glob_pattern in CNV_LOWRES_EVIDENCE_CFG.get("reference_npz_glob", []):
    if str(glob_pattern).strip():
        CNV_LOWRES_REF_NPZ_PATHS.extend(
            sorted(resolve_path(path) for path in __import__("glob").glob(resolve_path(glob_pattern)))
        )
CNV_LOWRES_REF_NPZ_PATHS = list(dict.fromkeys(CNV_LOWRES_REF_NPZ_PATHS))
CNV_LOWRES_REF_SAMPLE_IDS = [
    str(item)
    for item in CNV_LOWRES_EVIDENCE_CFG.get("reference_sample_ids", [])
    if str(item).strip()
]
if CNV_LOWRES_EVIDENCE_ENABLE and not CNV_LOWRES_REF_NPZ_PATHS:
    raise ValueError(
        "core.wisecondorx.cnv.lowres_evidence.enable=true requires "
        "lowres_evidence.reference_npz or lowres_evidence.reference_npz_glob"
    )
if (
    CNV_LOWRES_EVIDENCE_ENABLE
    and CNV_LOWRES_REF_SAMPLE_IDS
    and len(CNV_LOWRES_REF_SAMPLE_IDS) != len(CNV_LOWRES_REF_NPZ_PATHS)
):
    raise ValueError(
        "core.wisecondorx.cnv.lowres_evidence.reference_sample_ids must match lowres_evidence.reference_npz length"
    )
CNV_LOWRES_REF_MODERATE_MAD_Z = float(CNV_LOWRES_EVIDENCE_CFG.get("moderate_mad_z", 2.0))
CNV_LOWRES_REF_HIGH_MAD_Z = float(CNV_LOWRES_EVIDENCE_CFG.get("high_mad_z", 4.0))
CNV_B_REF_STABILITY_DIR = str(Path(CNV_POSTPROCESS_DIR) / "ref_stability")
CNV_B_REF_STABILITY_BINS = str(Path(CNV_B_REF_STABILITY_DIR) / "ref_bin_stability.tsv")
CNV_B_REF_STABILITY_SUMMARY = str(Path(CNV_B_REF_STABILITY_DIR) / "summary.json")
CNV_B_REF_STABILITY_EVENTS = str(Path(CNV_B_REF_STABILITY_DIR) / "{sample}.candidate_ref_stability.tsv")
CNV_B_REF_STABILITY_EVENTS_SUMMARY = str(Path(CNV_B_REF_STABILITY_DIR) / "{sample}.summary.json")
CNV_B_LOWRES_EVIDENCE = str(Path(CNV_POSTPROCESS_DIR) / "lowres_evidence" / "{sample}.candidate_evidence.tsv")
CNV_B_LOWRES_EVIDENCE_SUMMARY = str(Path(CNV_POSTPROCESS_DIR) / "lowres_evidence" / "{sample}.summary.json")
CNV_B_V2_CLASSIFIER = str(Path(CNV_POSTPROCESS_DIR) / "v2_classifier" / "{sample}.candidate_classification.tsv")
CNV_B_V2_CLASSIFIER_SUMMARY = str(Path(CNV_POSTPROCESS_DIR) / "v2_classifier" / "{sample}.summary.json")
CNV_B_V2_BENCHMARK_DIR = str(Path(CNV_POSTPROCESS_DIR) / "v2_benchmark")
CNV_B_V2_BENCHMARK_TRUTH_METRICS = str(Path(CNV_B_V2_BENCHMARK_DIR) / "truth_metrics.tsv")
CNV_B_V2_BENCHMARK_SAMPLE_SUMMARY = str(Path(CNV_B_V2_BENCHMARK_DIR) / "sample_summary.tsv")
CNV_B_V2_BENCHMARK_FILTERED_EVENTS = str(Path(CNV_B_V2_BENCHMARK_DIR) / "filtered_events.tsv")
CNV_B_V2_BENCHMARK_FILTERED_EVENTS_JSON = str(Path(CNV_B_V2_BENCHMARK_DIR) / "filtered_events.json")
CNV_B_V2_BENCHMARK_REPORT_EVENTS = str(Path(CNV_B_V2_BENCHMARK_DIR) / "report_events.tsv")
CNV_B_V2_BENCHMARK_REPORT_EVENTS_JSON = str(Path(CNV_B_V2_BENCHMARK_DIR) / "report_events.json")
CNV_B_V2_BENCHMARK_SUMMARY = str(Path(CNV_B_V2_BENCHMARK_DIR) / "summary.json")
CNV_BRANCH_S_EVIDENCE = str(Path(CNV_POSTPROCESS_DIR) / "branch_s" / "{sample}.sex_chrom_evidence.tsv")
CNV_BRANCH_S_SCORES = str(Path(CNV_POSTPROCESS_DIR) / "branch_s" / "{sample}.sca_state_scores.tsv")
CNV_BRANCH_S_SUMMARY = str(Path(CNV_POSTPROCESS_DIR) / "branch_s" / "{sample}.summary.json")
CNV_B_PLOT_SVG = str(Path(CNV_PLOTS_DIR) / "{sample}.final_cnv.svg")
CNV_NEGATIVE_BANK_CFG = CNV_CFG.get("negative_bank", {})
CNV_NEGATIVE_BANK_SAMPLES_TSV = (
    resolve_path(CNV_NEGATIVE_BANK_CFG.get("samples_tsv", ""))
    if str(CNV_NEGATIVE_BANK_CFG.get("samples_tsv", "")).strip()
    else ""
)
CNV_NEGATIVE_BANK_VERSION = str(CNV_NEGATIVE_BANK_CFG.get("version", "branch_ab_v2"))
CNV_NEGATIVE_BANK_BACKGROUND_LEDGERS = [
    resolve_path(item)
    for item in CNV_NEGATIVE_BANK_CFG.get("background_ledgers", [])
    if str(item).strip()
]
CNV_NEGATIVE_BANK_MIN_BACKGROUND = int(CNV_NEGATIVE_BANK_CFG.get("min_background", 5))
CNV_NEGATIVE_BANK_SIMILAR_LENGTH_FOLD = float(CNV_NEGATIVE_BANK_CFG.get("similar_length_fold", 2.0))
CNV_NEGATIVE_BANK_FEATURE_COLUMN = str(CNV_NEGATIVE_BANK_CFG.get("feature_column", "corrected_amplitude"))
CNV_NEGATIVE_BANK_SHADOW_BACKGROUND_LABELS = [
    str(item).strip().upper()
    for item in CNV_NEGATIVE_BANK_CFG.get("shadow_background_labels", [])
    if str(item).strip()
]
CNV_NEGATIVE_BANK_LABELS = str(Path(CNV_POSTPROCESS_DIR) / "negative_bank" / "negative_bank_labels.tsv")
CNV_NEGATIVE_BANK_SUMMARY = str(Path(CNV_POSTPROCESS_DIR) / "negative_bank" / "negative_bank_summary.json")
CNV_REFERENCE_ID = str(CNV_CFG.get("reference_id", "UNKNOWN_REFERENCE"))
CNV_CALLING_MIN_BINS = int(CNV_CALLING_CFG.get("min_bins", 5))
CNV_CALLING_MAX_SEGMENTS = int(CNV_CALLING_CFG.get("max_segments_per_chrom", 12))
CNV_CALLING_SPLIT_THRESHOLD = float(CNV_CALLING_CFG.get("split_threshold", 2.5))
CNV_CALLING_HMM_SHIFT = float(CNV_CALLING_CFG.get("hmm_state_shift", 2.5))
CNV_CALLING_HMM_STAY_PROB = float(CNV_CALLING_CFG.get("hmm_stay_prob", 0.995))
CNV_CALLING_HMM_ROLE = str(CNV_CALLING_CFG.get("hmm_role", "sidecar"))
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
CNV_ARTIFACT_A_BRANCH_REVIEW_MIN_ABS_Z = float(CNV_ARTIFACT_CFG.get("a_branch_review_min_abs_z", 15.0))
CNV_ARTIFACT_A_BRANCH_SENSITIVE_REVIEW_MIN_ABS_Z = float(
    CNV_ARTIFACT_CFG.get("a_branch_sensitive_review_min_abs_z", 7.0)
)
CNV_ARTIFACT_A_BRANCH_SENSITIVE_REVIEW_MIN_BINS = float(
    CNV_ARTIFACT_CFG.get("a_branch_sensitive_review_min_bins", 10.0)
)
CNV_ARTIFACT_A_BRANCH_SENSITIVE_REVIEW_MAX_HIGH_RISK_FRAC = float(
    CNV_ARTIFACT_CFG.get("a_branch_sensitive_review_max_high_risk_fraction", 0.20)
)
CNV_ARTIFACT_A_BRANCH_SENSITIVE_REVIEW_MAX_REGION_RISK = float(
    CNV_ARTIFACT_CFG.get("a_branch_sensitive_review_max_region_risk", 0.20)
)
CNV_ARTIFACT_A_BRANCH_SENSITIVE_REVIEW_MIN_SAME_DIRECTION_Z = float(
    CNV_ARTIFACT_CFG.get("a_branch_sensitive_review_min_same_direction_z", 0.25)
)
CNV_ARTIFACT_A_BRANCH_BOUNDARY_PROTECT_MIN_ABS_Z = float(
    CNV_ARTIFACT_CFG.get("a_branch_boundary_protect_min_abs_z", 30.0)
)
CNV_ARTIFACT_A_BRANCH_DISCORDANT_PROTECT_MIN_ABS_Z = float(
    CNV_ARTIFACT_CFG.get("a_branch_discordant_protect_min_abs_z", 50.0)
)
CNV_ARTIFACT_BRANCH_B_DIRECTION_MIN_ABS_Z = float(CNV_ARTIFACT_CFG.get("branch_b_direction_min_abs_z", 0.25))
CNV_ARTIFACT_NARROW_BOUNDARY_MAX_BINS = int(CNV_ARTIFACT_CFG.get("narrow_boundary_artifact_max_bins", 15))
CNV_ARTIFACT_NARROW_BOUNDARY_MAX_AVAILABLE_CHROM_FRACTION = float(
    CNV_ARTIFACT_CFG.get("narrow_boundary_artifact_max_available_chrom_fraction", 0.08)
)
CNV_ARTIFACT_NARROW_BOUNDARY_PROTECT_MIN_A_ABS_Z = float(
    CNV_ARTIFACT_CFG.get("narrow_boundary_artifact_protect_min_a_abs_z", 30.0)
)
CNV_ARTIFACT_SCA_XY_XGAIN_MAX_BAM_X_REL = float(CNV_ARTIFACT_CFG.get("sca_xy_xgain_max_bam_x_relative", 0.80))
CNV_ARTIFACT_SCA_XY_XGAIN_FOCAL_EDGE_MAX_BINS = int(
    CNV_ARTIFACT_CFG.get("sca_xy_xgain_focal_edge_max_bins", 20)
)
CNV_ARTIFACT_CNVSEQ_REPORTABLE_MIN_BP = int(CNV_ARTIFACT_CFG.get("cnvseq_reportable_min_bp", 2000000))
CNV_ARTIFACT_CNVSEQ_REVIEW_MIN_BP = int(CNV_ARTIFACT_CFG.get("cnvseq_review_min_bp", 1000000))
CNV_ARTIFACT_CNVSEQ_LARGE_EVENT_MIN_BP = int(CNV_ARTIFACT_CFG.get("cnvseq_large_event_min_bp", 10000000))
CNV_ARTIFACT_CNVSEQ_BOUNDARY_MAX_ABS_Z = float(CNV_ARTIFACT_CFG.get("cnvseq_boundary_max_abs_z", 4.0))
CNV_ARTIFACT_CNVSEQ_WHOLE_CHROM_AVAILABLE_FRACTION = float(
    CNV_ARTIFACT_CFG.get("cnvseq_whole_chrom_available_fraction", 0.90)
)
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
CNV_REPORT_BRANCH_A_VALIDATION_SUMMARIES = [
    resolve_path(item)
    for item in CNV_REPORT_CFG.get("branch_a_validation_summaries", [])
    if str(item).strip()
]
CNV_REFERENCE_AUDIT_CFG = CNV_CFG.get("reference_audit", {})
CNV_REFERENCE_AUDIT_SAMPLE_IDS = [
    str(item)
    for item in CNV_REFERENCE_AUDIT_CFG.get("samples", SAMPLES)
]
CNV_REFERENCE_AUDIT_SAMPLE_METADATA = (
    resolve_path(CNV_REFERENCE_AUDIT_CFG.get("sample_metadata", ""))
    if str(CNV_REFERENCE_AUDIT_CFG.get("sample_metadata", "")).strip()
    else ""
)
CNV_REFERENCE_AUDIT_REFERENCE_SAMPLES_FILE = (
    resolve_path(CNV_REFERENCE_AUDIT_CFG.get("reference_samples_file", ""))
    if str(CNV_REFERENCE_AUDIT_CFG.get("reference_samples_file", "")).strip()
    else ""
)
CNV_REFERENCE_AUDIT_BIN_ANNOTATIONS = (
    resolve_path(CNV_REFERENCE_AUDIT_CFG.get("bin_annotations", ""))
    if str(CNV_REFERENCE_AUDIT_CFG.get("bin_annotations", "")).strip()
    else ""
)
CNV_REFERENCE_AUDIT_EVIDENCE_LEDGER_DIR = (
    resolve_path(CNV_REFERENCE_AUDIT_CFG.get("evidence_ledger_dir", ""))
    if str(CNV_REFERENCE_AUDIT_CFG.get("evidence_ledger_dir", "")).strip()
    else str(Path(CNV_POSTPROCESS_DIR) / "evidence_ledger")
)
CNV_REFERENCE_AUDIT_STRONG_Z = float(CNV_REFERENCE_AUDIT_CFG.get("strong_z", 10.0))
CNV_REFERENCE_AUDIT_BROAD_EVENT_MIN_BP = int(CNV_REFERENCE_AUDIT_CFG.get("broad_event_min_bp", 10000000))
CNV_REFERENCE_AUDIT_SAMPLE_SPECIFIC_FRACTION_THRESHOLD = float(
    CNV_REFERENCE_AUDIT_CFG.get("sample_specific_fraction_threshold", 0.05)
)
CNV_REFERENCE_AUDIT_HIGH_RISK_BURDEN_THRESHOLD = float(
    CNV_REFERENCE_AUDIT_CFG.get("high_risk_burden_threshold", 0.20)
)
CNV_REFERENCE_AUDIT_SHARED_SIGNAL_MIN_SAMPLES = int(
    CNV_REFERENCE_AUDIT_CFG.get("shared_signal_min_samples", 3)
)
CNV_REFERENCE_AUDIT_DIR = str(Path(CNV_DIR) / "reference_audit")
CNV_REFERENCE_AUDIT_TSV = str(Path(CNV_REFERENCE_AUDIT_DIR) / "reference_candidate_audit.tsv")
CNV_REFERENCE_AUDIT_SUMMARY = str(Path(CNV_REFERENCE_AUDIT_DIR) / "reference_candidate_audit.summary.json")
CNV_NPZ = str(Path(CNV_DIR) / "npz" / "{sample}.npz")
CNV_MASKED_NPZ = str(Path(CNV_DIR) / "npz_mask_only" / "{sample}.npz")
CNV_MASK_SUMMARY = str(Path(CNV_DIR) / "npz_mask_only" / "{sample}.mask_summary.json")
CNV_PREDICT_INPUT_NPZ = CNV_MASKED_NPZ if CNV_PREPROCESS_STRATEGY == "mask_only" else CNV_NPZ
CNV_GENDER_TSV = str(Path(CNV_GENDER_DIR) / "{sample}.gender.tsv")
CNV_QC_TSV = str(Path(CNV_QC_DIR) / "{sample}.qc.tsv")
CNV_QC_PLOT = str(Path(CNV_QC_DIR) / "{sample}.qc.svg")
CNV_QC_PASS = str(Path(CNV_QC_DIR) / "{sample}.pass")
CNV_DONE = str(Path(CNV_PREDICT_DIR) / "{sample}.done")
CNV_A_ABERRATIONS_BED = str(Path(CNV_PREDICT_DIR) / "{sample}_aberrations.bed")
CNV_A_CANDIDATES = str(Path(CNV_BRANCH_A_OUTPUT_DIR) / "{sample}.candidate_events.tsv")
CNV_A_CANDIDATE_SUMMARY = str(Path(CNV_BRANCH_A_OUTPUT_DIR) / "{sample}.summary.json")
CNV_BRANCH_A_VALIDATION_SAMPLE_SUMMARY = str(Path(CNV_BRANCH_A_VALIDATION_DIR) / "sample_summary.tsv")
CNV_BRANCH_A_VALIDATION_TRUTH_METRICS = str(Path(CNV_BRANCH_A_VALIDATION_DIR) / "truth_metrics.tsv")
CNV_BRANCH_A_VALIDATION_SUMMARY = str(Path(CNV_BRANCH_A_VALIDATION_DIR) / "summary.json")
