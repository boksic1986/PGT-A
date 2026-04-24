RUN_METADATA = project_path("logs", "run_metadata.tsv")
MONITOR_DIR = project_path("monitor")
RUNTIME_TRACKING_CFG = config.get("runtime_tracking", {})
RUNTIME_TRACKING_ENABLED = bool(RUNTIME_TRACKING_CFG.get("enable", True))
RUNTIME_DB = str(Path(MONITOR_DIR) / "runtime.db")
RUNTIME_TRACKING_SENTINEL = str(Path(MONITOR_DIR) / "runtime_tracking.done")
BENCHMARK_ROOT = str(Path(MONITOR_DIR) / "benchmarks")
BENCH_FASTP_BWA = str(Path(BENCHMARK_ROOT) / "mapping" / "fastp_bwa" / "{sample}.tsv")
BENCH_BASELINE_BAM_UNIFORMITY_QC = str(Path(BENCHMARK_ROOT) / "qc" / "baseline_bam_uniformity_qc.tsv")
BENCH_AGGREGATE_BASELINE_QC = str(Path(BENCHMARK_ROOT) / "qc" / "aggregate_baseline_qc.tsv")
BENCH_BASELINE_MULTISCALE_BAM_UNIFORMITY_QC = str(Path(BENCHMARK_ROOT) / "qc" / "baseline_multiscale_bam_uniformity_qc.tsv")
BENCH_AGGREGATE_BASELINE_QC_MULTISCALE = str(Path(BENCHMARK_ROOT) / "qc" / "aggregate_baseline_qc_multiscale" / "{bin_label}.tsv")
BENCH_BASELINE_SAMPLE_SAMTOOLS_DIAGNOSTICS = str(Path(BENCHMARK_ROOT) / "qc" / "baseline_sample_samtools_diagnostics" / "{sample}.tsv")
BENCH_BASELINE_QC_REVIEW = str(Path(BENCHMARK_ROOT) / "qc" / "baseline_qc_review.tsv")
BENCH_REFERENCE_PREFILTER = str(Path(BENCHMARK_ROOT) / "reference" / "reference_prefilter.tsv")
BENCH_REFERENCE_PREFILTER_BY_SEX = {
    sex: str(Path(BENCHMARK_ROOT) / "reference" / "reference_prefilter" / f"{sex}.tsv")
    for sex in ("XX", "XY")
}
BENCH_TUNE_WISECONDORX_REFERENCE_QC = str(Path(BENCHMARK_ROOT) / "reference" / "tune_wisecondorx_reference_qc.tsv")
BENCH_BUILD_REFERENCE = str(Path(BENCHMARK_ROOT) / "reference" / "build_reference.tsv")
BENCH_BUILD_REFERENCE_BY_SEX = {
    sex: str(Path(BENCHMARK_ROOT) / "reference" / "build_reference" / f"{sex}.tsv")
    for sex in ("XX", "XY")
}
BENCH_BUILD_GENDER_REFERENCE = str(Path(BENCHMARK_ROOT) / "reference" / "build_gender_reference.tsv")
BENCH_SELECT_REFERENCE_COHORT = {
    ref_set: str(Path(BENCHMARK_ROOT) / "reference" / "select_cohort" / f"{ref_set}.tsv")
    for ref_set in REF_SET_ORDER
}
BENCH_REFERENCE_PREFILTER_BY_SET_SEX = {
    ref_set: {
        sex: str(Path(BENCHMARK_ROOT) / "reference" / "reference_prefilter" / ref_set / f"{sex}.tsv")
        for sex in REF_SEXES
    }
    for ref_set in REF_SET_ORDER
}
BENCH_TUNE_WISECONDORX_REFERENCE_QC_BY_SET = {
    ref_set: str(Path(BENCHMARK_ROOT) / "reference" / "tune_wisecondorx_reference_qc" / f"{ref_set}.tsv")
    for ref_set in REF_SET_ORDER
}
BENCH_BUILD_REFERENCE_BY_SET_SEX = {
    ref_set: {
        sex: str(Path(BENCHMARK_ROOT) / "reference" / "build_reference" / ref_set / f"{sex}.tsv")
        for sex in REF_SEXES
    }
    for ref_set in REF_SET_ORDER
}
BENCH_BUILD_GENDER_REFERENCE_BY_SET = {
    ref_set: str(Path(BENCHMARK_ROOT) / "reference" / "build_gender_reference" / f"{ref_set}.tsv")
    for ref_set in REF_SET_ORDER
}
BENCH_WISECONDORX_CONVERT_FOR_CNV = str(Path(BENCHMARK_ROOT) / "predict" / "convert" / "{sample}.tsv")
BENCH_WISECONDORX_GENDER_FOR_PREDICT = str(Path(BENCHMARK_ROOT) / "predict" / "gender" / "{sample}.tsv")
BENCH_WISECONDORX_QC_FOR_PREDICT = str(Path(BENCHMARK_ROOT) / "predict" / "qc" / "{sample}.tsv")
BENCH_WISECONDORX_PREDICT_CNV = str(Path(BENCHMARK_ROOT) / "predict" / "predict" / "{sample}.tsv")
