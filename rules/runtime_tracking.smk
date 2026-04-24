from pathlib import Path


def _sample_bucket(prefix, sample_count):
    return f"{prefix}_n{sample_count}"


def _runtime_record(
    *,
    module_name,
    rule_name,
    snakemake_target,
    benchmark_path,
    log_path,
    output_paths,
    sample_id=None,
    cohort_type=None,
    sample_bucket=None,
    binsize=None,
    pca_components=None,
    threads=None,
    parameter_path=None,
):
    return {
        "module": module_name,
        "rule_name": rule_name,
        "snakemake_target": snakemake_target,
        "sample_id": sample_id,
        "cohort_type": cohort_type,
        "sample_bucket": sample_bucket,
        "binsize": binsize,
        "pca_components": pca_components,
        "threads": threads,
        "benchmark_path": benchmark_path,
        "log_path": log_path,
        "output_paths": list(output_paths),
        "parameter_path": parameter_path,
    }


RUNTIME_TRACKING_RECORDS = []
BASELINE_BUCKET = _sample_bucket("baseline_cohort", len(BASELINE_SAMPLE_IDS))
REFERENCE_BUCKET = _sample_bucket("reference_cohort", len(REF_SAMPLE_IDS))

if RUNTIME_TRACKING_ENABLED and RUNTIME_BENCHMARK_FILES:
    for sample_id in SAMPLES:
        RUNTIME_TRACKING_RECORDS.append(
            _runtime_record(
                module_name="mapping",
                rule_name="fastp_bwa",
                snakemake_target="mapping",
                sample_id=sample_id,
                cohort_type="mixed",
                sample_bucket="single_sample",
                benchmark_path=BENCH_FASTP_BWA.format(sample=sample_id),
                log_path=project_path("logs", "bwa", f"{sample_id}.log"),
                output_paths=[
                    FASTP_R1.format(sample=sample_id),
                    FASTP_R2.format(sample=sample_id),
                    FASTP_HTML.format(sample=sample_id),
                    FASTP_JSON.format(sample=sample_id),
                    SORTED_BAM.format(sample=sample_id),
                    SORTED_BAI.format(sample=sample_id),
                ],
                threads=16,
            )
        )

    if "baseline_qc" in REQUESTED_TARGETS:
        RUNTIME_TRACKING_RECORDS.extend(
            [
                _runtime_record(
                    module_name="qc",
                    rule_name="baseline_bam_uniformity_qc",
                    snakemake_target="baseline_qc",
                    cohort_type="baseline",
                    sample_bucket=BASELINE_BUCKET,
                    benchmark_path=BENCH_BASELINE_BAM_UNIFORMITY_QC,
                    log_path=project_path("logs", "qc", "baseline", "batch_qc.log"),
                    output_paths=(
                        expand(BASELINE_QC_TSV, sample=BASELINE_SAMPLE_IDS)
                        + expand(BASELINE_QC_PROFILE_TSV, sample=BASELINE_SAMPLE_IDS)
                        + expand(BASELINE_QC_REF_SUMMARY_TSV, sample=BASELINE_SAMPLE_IDS)
                        + expand(BASELINE_QC_PLOT, sample=BASELINE_SAMPLE_IDS)
                        + expand(BASELINE_QC_GC_PLOT, sample=BASELINE_SAMPLE_IDS)
                    ),
                    binsize=BASELINE_QC_BIN_SIZE,
                    threads=BASELINE_QC_THREADS,
                ),
                _runtime_record(
                    module_name="qc",
                    rule_name="aggregate_baseline_qc",
                    snakemake_target="baseline_qc",
                    cohort_type="baseline",
                    sample_bucket=BASELINE_BUCKET,
                    benchmark_path=BENCH_AGGREGATE_BASELINE_QC,
                    log_path=project_path("logs", "qc", "baseline", "aggregate.log"),
                    output_paths=[
                        BASELINE_QC_SUMMARY,
                        BASELINE_QC_PASS_SAMPLES,
                        BASELINE_QC_RETAINED_SAMPLES,
                        BASELINE_QC_OUTLIER_SAMPLES,
                        BASELINE_QC_REPORT_MD,
                        BASELINE_QC_FIG_DECISION,
                        BASELINE_QC_FIG_MAPPED,
                        BASELINE_QC_FIG_METRICS,
                        BASELINE_QC_FIG_BATCH,
                        BASELINE_QC_FIG_CORR,
                    ],
                    binsize=BASELINE_QC_BIN_SIZE,
                    threads=1,
                ),
                _runtime_record(
                    module_name="qc",
                    rule_name="baseline_multiscale_bam_uniformity_qc",
                    snakemake_target="baseline_qc",
                    cohort_type="baseline",
                    sample_bucket=BASELINE_BUCKET,
                    benchmark_path=BENCH_BASELINE_MULTISCALE_BAM_UNIFORMITY_QC,
                    log_path=project_path("logs", "qc", "baseline_multiscale", "batch_qc.log"),
                    output_paths=(
                        expand(BASELINE_MULTI_TSV, bin_label=BASELINE_QC_BIN_LABELS, sample=BASELINE_SAMPLE_IDS)
                        + expand(BASELINE_MULTI_PROFILE_TSV, bin_label=BASELINE_QC_BIN_LABELS, sample=BASELINE_SAMPLE_IDS)
                        + expand(BASELINE_MULTI_REF_SUMMARY_TSV, bin_label=BASELINE_QC_BIN_LABELS, sample=BASELINE_SAMPLE_IDS)
                        + expand(BASELINE_MULTI_PLOT, bin_label=BASELINE_QC_BIN_LABELS, sample=BASELINE_SAMPLE_IDS)
                        + expand(BASELINE_MULTI_GC_PLOT, bin_label=BASELINE_QC_BIN_LABELS, sample=BASELINE_SAMPLE_IDS)
                    ),
                    threads=BASELINE_QC_THREADS,
                ),
                _runtime_record(
                    module_name="qc",
                    rule_name="baseline_qc_review",
                    snakemake_target="baseline_qc",
                    cohort_type="baseline",
                    sample_bucket=BASELINE_BUCKET,
                    benchmark_path=BENCH_BASELINE_QC_REVIEW,
                    log_path=project_path("logs", "qc", "baseline_multiscale", "report.log"),
                    output_paths=[
                        BASELINE_MULTI_BINSIZE_SUMMARY,
                        BASELINE_MULTI_SAMPLE_MATRIX,
                        BASELINE_MULTI_STABLE_RETAINED,
                        BASELINE_MULTI_FAILED_DIAG,
                        BASELINE_REPORT_MD,
                        BASELINE_MULTI_FIG_DECISION_TREND,
                        BASELINE_MULTI_FIG_SAMPLE_HEATMAP,
                        BASELINE_MULTI_FIG_LIBRARY,
                        BASELINE_MULTI_FIG_FAIL_CHROM,
                    ],
                    threads=1,
                ),
            ]
        )

        for bin_label, bin_size in zip(BASELINE_QC_BIN_LABELS, BASELINE_QC_BIN_SIZES):
            RUNTIME_TRACKING_RECORDS.append(
                _runtime_record(
                    module_name="qc",
                    rule_name="aggregate_baseline_qc_multiscale",
                    snakemake_target="baseline_qc",
                    cohort_type="baseline",
                    sample_bucket=BASELINE_BUCKET,
                    benchmark_path=BENCH_AGGREGATE_BASELINE_QC_MULTISCALE.format(bin_label=bin_label),
                    log_path=project_path("logs", "qc", "baseline_multiscale", bin_label, "aggregate.log"),
                    output_paths=[
                        BASELINE_MULTI_SUMMARY.format(bin_label=bin_label),
                        BASELINE_MULTI_PASS_SAMPLES.format(bin_label=bin_label),
                        BASELINE_MULTI_RETAINED_SAMPLES.format(bin_label=bin_label),
                        BASELINE_MULTI_OUTLIER_SAMPLES.format(bin_label=bin_label),
                        BASELINE_MULTI_REPORT_MD.format(bin_label=bin_label),
                        BASELINE_MULTI_FIG_DECISION.format(bin_label=bin_label),
                        BASELINE_MULTI_FIG_MAPPED.format(bin_label=bin_label),
                        BASELINE_MULTI_FIG_METRICS.format(bin_label=bin_label),
                        BASELINE_MULTI_FIG_BATCH.format(bin_label=bin_label),
                        BASELINE_MULTI_FIG_CORR.format(bin_label=bin_label),
                    ],
                    binsize=bin_size,
                    threads=1,
                )
            )

        for sample_id in BASELINE_SAMPLE_IDS:
            RUNTIME_TRACKING_RECORDS.append(
                _runtime_record(
                    module_name="qc",
                    rule_name="baseline_sample_samtools_diagnostics",
                    snakemake_target="baseline_qc",
                    sample_id=sample_id,
                    cohort_type="baseline",
                    sample_bucket="single_sample",
                    benchmark_path=BENCH_BASELINE_SAMPLE_SAMTOOLS_DIAGNOSTICS.format(sample=sample_id),
                    log_path=project_path("logs", "qc", "baseline_diagnostics", f"{sample_id}.log"),
                    output_paths=[
                        BASELINE_FLAGSTAT.format(sample=sample_id),
                        BASELINE_IDXSTATS.format(sample=sample_id),
                    ],
                    threads=BASELINE_QC_SAMTOOLS_THREADS,
                )
            )

    if "reference_qc" in REQUESTED_TARGETS and TUNING_ENABLED:
        if BUILD_REF_BY_SEX_ENABLED:
            for sex in ("XX", "XY"):
                sex_ids = REF_SAMPLE_IDS_BY_SEX[sex]
                RUNTIME_TRACKING_RECORDS.append(
                    _runtime_record(
                        module_name="reference",
                        rule_name=f"reference_prefilter_{sex.lower()}",
                        snakemake_target="reference_qc",
                        cohort_type=f"reference_{sex}",
                        sample_bucket=_sample_bucket(f"reference_{sex}", len(sex_ids)),
                        benchmark_path=BENCH_REFERENCE_PREFILTER_BY_SEX[sex],
                        log_path=project_path("logs", "wisecondorx", f"reference_prefilter_{sex}.log"),
                        output_paths=[
                            str(Path(REF_PREFILTER_DIR_BY_SEX[sex]) / "reference_sample_qc.tsv"),
                            str(Path(REF_PREFILTER_DIR_BY_SEX[sex]) / "reference_sample_qc.svg"),
                            str(Path(REF_PREFILTER_DIR_BY_SEX[sex]) / "reference_inlier_samples.txt"),
                            str(Path(REF_PREFILTER_DIR_BY_SEX[sex]) / "prefilter_summary.yaml"),
                        ],
                        binsize=PREFILTER_BINSIZE,
                        threads=4,
                    )
                )
        else:
            RUNTIME_TRACKING_RECORDS.append(
                _runtime_record(
                    module_name="reference",
                    rule_name="reference_prefilter",
                    snakemake_target="reference_qc",
                    cohort_type="reference",
                    sample_bucket=REFERENCE_BUCKET,
                    benchmark_path=BENCH_REFERENCE_PREFILTER,
                    log_path=project_path("logs", "wisecondorx", "reference_prefilter.log"),
                    output_paths=[
                        str(Path(REF_PREFILTER_DIR) / "reference_sample_qc.tsv"),
                        str(Path(REF_PREFILTER_DIR) / "reference_sample_qc.svg"),
                        str(Path(REF_PREFILTER_DIR) / "reference_inlier_samples.txt"),
                        str(Path(REF_PREFILTER_DIR) / "prefilter_summary.yaml"),
                    ],
                    binsize=PREFILTER_BINSIZE,
                    threads=4,
                )
            )

        RUNTIME_TRACKING_RECORDS.append(
            _runtime_record(
                module_name="reference",
                rule_name="tune_wisecondorx_reference_qc",
                snakemake_target="reference_qc",
                cohort_type="reference",
                sample_bucket=REFERENCE_BUCKET,
                benchmark_path=BENCH_TUNE_WISECONDORX_REFERENCE_QC,
                log_path=project_path("logs", "wisecondorx", "tuning.log"),
                output_paths=[TUNING_SUMMARY, TUNING_BEST, TUNING_QC, TUNING_PLOT, TUNING_QC_STATS_PLOT, TUNING_INLIERS],
                threads=4,
                parameter_path=TUNING_BEST,
            )
        )

    if "reference" in REQUESTED_TARGETS:
        if BUILD_REF_BY_SEX_ENABLED:
            for sex in ("XX", "XY"):
                sex_ids = REF_SAMPLE_IDS_BY_SEX[sex]
                parameter_path = TUNING_BEST if TUNING_ENABLED else None
                binsize = None if TUNING_ENABLED else int(WISE_CFG["binsize"])
                RUNTIME_TRACKING_RECORDS.append(
                    _runtime_record(
                        module_name="reference",
                        rule_name=f"build_wisecondorx_reference_{sex.lower()}",
                        snakemake_target="reference",
                        cohort_type=f"reference_{sex}",
                        sample_bucket=_sample_bucket(f"reference_{sex}", len(sex_ids)),
                        benchmark_path=BENCH_BUILD_REFERENCE_BY_SEX[sex],
                        log_path=project_path("logs", "wisecondorx", f"build_reference_{sex}.log"),
                        output_paths=[REF_OUTPUTS_BY_SEX[sex]],
                        binsize=binsize,
                        threads=4,
                        parameter_path=parameter_path,
                    )
                )

            RUNTIME_TRACKING_RECORDS.append(
                _runtime_record(
                    module_name="reference",
                    rule_name="build_wisecondorx_gender_reference",
                    snakemake_target="reference",
                    cohort_type="reference",
                    sample_bucket=REFERENCE_BUCKET,
                    benchmark_path=BENCH_BUILD_GENDER_REFERENCE,
                    log_path=project_path("logs", "wisecondorx", "build_reference_gender.log"),
                    output_paths=[GENDER_REF_OUTPUT],
                    binsize=None if TUNING_ENABLED else int(WISE_CFG["binsize"]),
                    threads=4,
                    parameter_path=TUNING_BEST if TUNING_ENABLED else COMMON_REF_BINSIZE,
                )
            )
        else:
            RUNTIME_TRACKING_RECORDS.append(
                _runtime_record(
                    module_name="reference",
                    rule_name="build_wisecondorx_reference",
                    snakemake_target="reference",
                    cohort_type="reference",
                    sample_bucket=REFERENCE_BUCKET,
                    benchmark_path=BENCH_BUILD_REFERENCE,
                    log_path=project_path("logs", "wisecondorx", "build_reference.log"),
                    output_paths=[REF_OUTPUT],
                    binsize=None if TUNING_ENABLED else int(WISE_CFG["binsize"]),
                    threads=4,
                    parameter_path=TUNING_BEST if TUNING_ENABLED else None,
                )
            )

    if "cnv_qc" in REQUESTED_TARGETS and CNV_ENABLED:
        for sample_id in SAMPLES:
            RUNTIME_TRACKING_RECORDS.append(
                _runtime_record(
                    module_name="predict",
                    rule_name="wisecondorx_convert_for_cnv",
                    snakemake_target="cnv_qc",
                    sample_id=sample_id,
                    cohort_type="abnormal",
                    sample_bucket="single_sample",
                    benchmark_path=BENCH_WISECONDORX_CONVERT_FOR_CNV.format(sample=sample_id),
                    log_path=project_path("logs", "cnv", f"{sample_id}.convert.log"),
                    output_paths=[CNV_NPZ.format(sample=sample_id)],
                    binsize=None if PREDICT_BY_SEX_ENABLED else CNV_CONVERT_BINSIZE,
                    threads=4,
                    parameter_path=COMMON_REF_BINSIZE if PREDICT_BY_SEX_ENABLED else None,
                )
            )
            if PREDICT_BY_SEX_ENABLED:
                RUNTIME_TRACKING_RECORDS.append(
                    _runtime_record(
                        module_name="predict",
                        rule_name="wisecondorx_gender_for_predict",
                        snakemake_target="cnv_qc",
                        sample_id=sample_id,
                        cohort_type="abnormal",
                        sample_bucket="single_sample",
                        benchmark_path=BENCH_WISECONDORX_GENDER_FOR_PREDICT.format(sample=sample_id),
                        log_path=project_path("logs", "cnv", f"{sample_id}.gender.log"),
                        output_paths=[CNV_GENDER_TSV.format(sample=sample_id)],
                        threads=1,
                    )
                )
            RUNTIME_TRACKING_RECORDS.append(
                _runtime_record(
                    module_name="predict",
                    rule_name="wisecondorx_qc_for_predict",
                    snakemake_target="cnv_qc",
                    sample_id=sample_id,
                    cohort_type="abnormal",
                    sample_bucket="single_sample",
                    benchmark_path=BENCH_WISECONDORX_QC_FOR_PREDICT.format(sample=sample_id),
                    log_path=project_path("logs", "cnv", f"{sample_id}.qc.log"),
                    output_paths=[
                        CNV_QC_TSV.format(sample=sample_id),
                        CNV_QC_PLOT.format(sample=sample_id),
                        CNV_QC_PASS.format(sample=sample_id),
                    ],
                    threads=1,
                )
            )

    if "cnv" in REQUESTED_TARGETS and CNV_ENABLED:
        for sample_id in SAMPLES:
            RUNTIME_TRACKING_RECORDS.append(
                _runtime_record(
                    module_name="predict",
                    rule_name="wisecondorx_predict_cnv",
                    snakemake_target="cnv",
                    sample_id=sample_id,
                    cohort_type="abnormal",
                    sample_bucket="single_sample",
                    benchmark_path=BENCH_WISECONDORX_PREDICT_CNV.format(sample=sample_id),
                    log_path=project_path("logs", "cnv", f"{sample_id}.predict.log"),
                    output_paths=[CNV_DONE.format(sample=sample_id)],
                    threads=2,
                    parameter_path=COMMON_REF_BINSIZE if PREDICT_BY_SEX_ENABLED else None,
                )
            )

    rule collect_runtime_tracking:
        input:
            upstream=RUNTIME_TRACKING_INPUT_FILES
        output:
            db=RUNTIME_DB,
            done=RUNTIME_TRACKING_SENTINEL
        params:
            records=RUNTIME_TRACKING_RECORDS
        threads: 1
        script:
            SCRIPT_COLLECT_BENCHMARKS_TO_SQLITE
