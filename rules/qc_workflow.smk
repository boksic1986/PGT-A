if "baseline_qc" in REQUESTED_TARGETS:
    rule baseline_bam_uniformity_qc:
        input:
            bams=expand(SORTED_BAM, sample=BASELINE_SAMPLE_IDS),
            metadata=RUN_METADATA
        output:
            qc_tsvs=expand(BASELINE_QC_TSV, sample=BASELINE_SAMPLE_IDS),
            profiles=expand(BASELINE_QC_PROFILE_TSV, sample=BASELINE_SAMPLE_IDS),
            ref_summaries=expand(BASELINE_QC_REF_SUMMARY_TSV, sample=BASELINE_SAMPLE_IDS),
            plots=expand(BASELINE_QC_PLOT, sample=BASELINE_SAMPLE_IDS),
            gc_plots=expand(BASELINE_QC_GC_PLOT, sample=BASELINE_SAMPLE_IDS)
        log:
            project_path("logs", "qc", "baseline", "batch_qc.log")
        benchmark:
            BENCH_BASELINE_BAM_UNIFORMITY_QC
        params:
            python_bin=config["biosoft"]["python"],
            script=SCRIPT_BATCH_BAM_UNIFORMITY_QC,
            script_action=SCRIPT_BATCH_BAM_UNIFORMITY_QC_ACTION,
            outdir=BASELINE_QC_DIR,
            reference_fasta=config["core"]["reference_genome"],
            bin_size=BASELINE_QC_BIN_SIZE,
            min_mapped_warn=BASELINE_QC_THRESHOLDS["min_mapped_warn"],
            min_mapped_fail=BASELINE_QC_THRESHOLDS["min_mapped_fail"],
            max_zero_frac_warn=BASELINE_QC_THRESHOLDS["max_zero_frac_warn"],
            max_zero_frac_fail=BASELINE_QC_THRESHOLDS["max_zero_frac_fail"],
            max_bin_cv_warn=BASELINE_QC_THRESHOLDS["max_bin_cv_warn"],
            max_bin_cv_fail=BASELINE_QC_THRESHOLDS["max_bin_cv_fail"],
            max_adj_mad_warn=BASELINE_QC_THRESHOLDS["max_adj_mad_warn"],
            max_adj_mad_fail=BASELINE_QC_THRESHOLDS["max_adj_mad_fail"],
            max_gini_warn=BASELINE_QC_THRESHOLDS["max_gini_warn"],
            max_gini_fail=BASELINE_QC_THRESHOLDS["max_gini_fail"],
            min_pearson_warn=BASELINE_QC_THRESHOLDS["min_pearson_warn"],
            min_pearson_fail=BASELINE_QC_THRESHOLDS["min_pearson_fail"],
            min_spearman_warn=BASELINE_QC_THRESHOLDS["min_spearman_warn"],
            min_spearman_fail=BASELINE_QC_THRESHOLDS["min_spearman_fail"],
            max_median_abs_z_warn=BASELINE_QC_THRESHOLDS["max_median_abs_z_warn"],
            max_median_abs_z_fail=BASELINE_QC_THRESHOLDS["max_median_abs_z_fail"],
            max_outlier3_warn=BASELINE_QC_THRESHOLDS["max_outlier3_warn"],
            max_outlier3_fail=BASELINE_QC_THRESHOLDS["max_outlier3_fail"],
            gc_correction_method=BASELINE_QC_GC_CORRECTION_METHOD,
            gc_correction_frac=BASELINE_QC_GC_CORRECTION_FRAC,
            gc_correction_poly_degree=BASELINE_QC_GC_CORRECTION_POLY_DEGREE,
            gc_correction_min_valid_bins=BASELINE_QC_GC_CORRECTION_MIN_VALID_BINS,
            gc_correction_robust_iters=BASELINE_QC_GC_CORRECTION_ROBUST_ITERS
        threads: BASELINE_QC_THREADS
        shell:
            r"""
            mkdir -p "{params.outdir}" "$(dirname {log})"
            (
                echo "=== PIPELINE AUDIT ==="
                cat {input.metadata:q}
                echo "=== COMMAND ==="
            ) > {log:q}
            export NUMEXPR_MAX_THREADS={threads}
            export NUMEXPR_NUM_THREADS={threads}
            {params.python_bin:q} {params.script:q} {params.script_action:q} \
                --bams {input.bams:q} \
                --reference-fasta {params.reference_fasta:q} \
                --bin-size {params.bin_size:q} \
                --threads {threads} \
                --min-mapped-warn {params.min_mapped_warn:q} \
                --min-mapped-fail {params.min_mapped_fail:q} \
                --max-zero-frac-warn {params.max_zero_frac_warn:q} \
                --max-zero-frac-fail {params.max_zero_frac_fail:q} \
                --max-bin-cv-warn {params.max_bin_cv_warn:q} \
                --max-bin-cv-fail {params.max_bin_cv_fail:q} \
                --max-adj-mad-warn {params.max_adj_mad_warn:q} \
                --max-adj-mad-fail {params.max_adj_mad_fail:q} \
                --max-gini-warn {params.max_gini_warn:q} \
                --max-gini-fail {params.max_gini_fail:q} \
                --min-pearson-warn {params.min_pearson_warn:q} \
                --min-pearson-fail {params.min_pearson_fail:q} \
                --min-spearman-warn {params.min_spearman_warn:q} \
                --min-spearman-fail {params.min_spearman_fail:q} \
                --max-median-abs-z-warn {params.max_median_abs_z_warn:q} \
                --max-median-abs-z-fail {params.max_median_abs_z_fail:q} \
                --max-outlier3-warn {params.max_outlier3_warn:q} \
                --max-outlier3-fail {params.max_outlier3_fail:q} \
                --gc-correction-method {params.gc_correction_method:q} \
                --gc-correction-frac {params.gc_correction_frac:q} \
                --gc-correction-poly-degree {params.gc_correction_poly_degree:q} \
                --gc-correction-min-valid-bins {params.gc_correction_min_valid_bins:q} \
                --gc-correction-robust-iters {params.gc_correction_robust_iters:q} \
                --outdir {params.outdir:q} \
                --log {log:q} \
                >> {log:q} 2>&1
            """

    rule aggregate_baseline_qc:
        input:
            qc_tsvs=expand(BASELINE_QC_TSV, sample=BASELINE_SAMPLE_IDS),
            profile_tsvs=expand(BASELINE_QC_PROFILE_TSV, sample=BASELINE_SAMPLE_IDS),
            metadata=RUN_METADATA
        output:
            summary=BASELINE_QC_SUMMARY,
            passed=BASELINE_QC_PASS_SAMPLES,
            retained=BASELINE_QC_RETAINED_SAMPLES,
            outliers=BASELINE_QC_OUTLIER_SAMPLES,
            report=BASELINE_QC_REPORT_MD,
            fig_decision=BASELINE_QC_FIG_DECISION,
            fig_mapped=BASELINE_QC_FIG_MAPPED,
            fig_metrics=BASELINE_QC_FIG_METRICS,
            fig_batch=BASELINE_QC_FIG_BATCH,
            fig_corr=BASELINE_QC_FIG_CORR
        log:
            project_path("logs", "qc", "baseline", "aggregate.log")
        benchmark:
            BENCH_AGGREGATE_BASELINE_QC
        params:
            python_bin=config["biosoft"]["python"],
            script=SCRIPT_AGGREGATE_BASELINE_QC,
            script_action=SCRIPT_AGGREGATE_BASELINE_QC_ACTION,
            figures_dir=BASELINE_QC_FIGURES_DIR
        threads: 1
        shell:
            r"""
            mkdir -p "$(dirname {output.summary})" "$(dirname {log})" "{params.figures_dir}"
            (
                echo "=== PIPELINE AUDIT ==="
                cat {input.metadata:q}
                echo "=== COMMAND ==="
            ) > {log:q}
            {params.python_bin:q} {params.script:q} {params.script_action:q} \
                --qc-tsvs {input.qc_tsvs:q} \
                --profile-tsvs {input.profile_tsvs:q} \
                --summary-output {output.summary:q} \
                --pass-samples-output {output.passed:q} \
                --retained-samples-output {output.retained:q} \
                --outlier-samples-output {output.outliers:q} \
                --report-output {output.report:q} \
                --figures-dir {params.figures_dir:q} \
                --log {log:q} \
                >> {log:q} 2>&1
            """

    rule baseline_multiscale_bam_uniformity_qc:
        input:
            bams=expand(SORTED_BAM, sample=BASELINE_SAMPLE_IDS),
            metadata=RUN_METADATA
        output:
            qc_tsvs=expand(BASELINE_MULTI_TSV, bin_label=BASELINE_QC_BIN_LABELS, sample=BASELINE_SAMPLE_IDS),
            profiles=expand(BASELINE_MULTI_PROFILE_TSV, bin_label=BASELINE_QC_BIN_LABELS, sample=BASELINE_SAMPLE_IDS),
            ref_summaries=expand(BASELINE_MULTI_REF_SUMMARY_TSV, bin_label=BASELINE_QC_BIN_LABELS, sample=BASELINE_SAMPLE_IDS),
            plots=expand(BASELINE_MULTI_PLOT, bin_label=BASELINE_QC_BIN_LABELS, sample=BASELINE_SAMPLE_IDS),
            gc_plots=expand(BASELINE_MULTI_GC_PLOT, bin_label=BASELINE_QC_BIN_LABELS, sample=BASELINE_SAMPLE_IDS)
        log:
            project_path("logs", "qc", "baseline_multiscale", "batch_qc.log")
        benchmark:
            BENCH_BASELINE_MULTISCALE_BAM_UNIFORMITY_QC
        params:
            python_bin=config["biosoft"]["python"],
            script=SCRIPT_BATCH_BAM_UNIFORMITY_QC,
            script_action=SCRIPT_BATCH_BAM_UNIFORMITY_QC_ACTION,
            outdir=BASELINE_MULTI_DIR,
            reference_fasta=config["core"]["reference_genome"],
            bin_sizes=" ".join(str(item) for item in BASELINE_QC_BIN_SIZES),
            min_mapped_warn=BASELINE_QC_THRESHOLDS["min_mapped_warn"],
            min_mapped_fail=BASELINE_QC_THRESHOLDS["min_mapped_fail"],
            max_zero_frac_warn=BASELINE_QC_THRESHOLDS["max_zero_frac_warn"],
            max_zero_frac_fail=BASELINE_QC_THRESHOLDS["max_zero_frac_fail"],
            max_bin_cv_warn=BASELINE_QC_THRESHOLDS["max_bin_cv_warn"],
            max_bin_cv_fail=BASELINE_QC_THRESHOLDS["max_bin_cv_fail"],
            max_adj_mad_warn=BASELINE_QC_THRESHOLDS["max_adj_mad_warn"],
            max_adj_mad_fail=BASELINE_QC_THRESHOLDS["max_adj_mad_fail"],
            max_gini_warn=BASELINE_QC_THRESHOLDS["max_gini_warn"],
            max_gini_fail=BASELINE_QC_THRESHOLDS["max_gini_fail"],
            min_pearson_warn=BASELINE_QC_THRESHOLDS["min_pearson_warn"],
            min_pearson_fail=BASELINE_QC_THRESHOLDS["min_pearson_fail"],
            min_spearman_warn=BASELINE_QC_THRESHOLDS["min_spearman_warn"],
            min_spearman_fail=BASELINE_QC_THRESHOLDS["min_spearman_fail"],
            max_median_abs_z_warn=BASELINE_QC_THRESHOLDS["max_median_abs_z_warn"],
            max_median_abs_z_fail=BASELINE_QC_THRESHOLDS["max_median_abs_z_fail"],
            max_outlier3_warn=BASELINE_QC_THRESHOLDS["max_outlier3_warn"],
            max_outlier3_fail=BASELINE_QC_THRESHOLDS["max_outlier3_fail"],
            gc_correction_method=BASELINE_QC_GC_CORRECTION_METHOD,
            gc_correction_frac=BASELINE_QC_GC_CORRECTION_FRAC,
            gc_correction_poly_degree=BASELINE_QC_GC_CORRECTION_POLY_DEGREE,
            gc_correction_min_valid_bins=BASELINE_QC_GC_CORRECTION_MIN_VALID_BINS,
            gc_correction_robust_iters=BASELINE_QC_GC_CORRECTION_ROBUST_ITERS
        threads: BASELINE_QC_THREADS
        shell:
            r"""
            mkdir -p "{params.outdir}" "$(dirname {log})"
            (
                echo "=== PIPELINE AUDIT ==="
                cat {input.metadata:q}
                echo "=== COMMAND ==="
            ) > {log:q}
            export NUMEXPR_MAX_THREADS={threads}
            export NUMEXPR_NUM_THREADS={threads}
            {params.python_bin:q} {params.script:q} {params.script_action:q} \
                --bams {input.bams:q} \
                --reference-fasta {params.reference_fasta:q} \
                --bin-sizes {params.bin_sizes} \
                --threads {threads} \
                --min-mapped-warn {params.min_mapped_warn:q} \
                --min-mapped-fail {params.min_mapped_fail:q} \
                --max-zero-frac-warn {params.max_zero_frac_warn:q} \
                --max-zero-frac-fail {params.max_zero_frac_fail:q} \
                --max-bin-cv-warn {params.max_bin_cv_warn:q} \
                --max-bin-cv-fail {params.max_bin_cv_fail:q} \
                --max-adj-mad-warn {params.max_adj_mad_warn:q} \
                --max-adj-mad-fail {params.max_adj_mad_fail:q} \
                --max-gini-warn {params.max_gini_warn:q} \
                --max-gini-fail {params.max_gini_fail:q} \
                --min-pearson-warn {params.min_pearson_warn:q} \
                --min-pearson-fail {params.min_pearson_fail:q} \
                --min-spearman-warn {params.min_spearman_warn:q} \
                --min-spearman-fail {params.min_spearman_fail:q} \
                --max-median-abs-z-warn {params.max_median_abs_z_warn:q} \
                --max-median-abs-z-fail {params.max_median_abs_z_fail:q} \
                --max-outlier3-warn {params.max_outlier3_warn:q} \
                --max-outlier3-fail {params.max_outlier3_fail:q} \
                --gc-correction-method {params.gc_correction_method:q} \
                --gc-correction-frac {params.gc_correction_frac:q} \
                --gc-correction-poly-degree {params.gc_correction_poly_degree:q} \
                --gc-correction-min-valid-bins {params.gc_correction_min_valid_bins:q} \
                --gc-correction-robust-iters {params.gc_correction_robust_iters:q} \
                --outdir {params.outdir:q} \
                --log {log:q} \
                >> {log:q} 2>&1
            """

    rule aggregate_baseline_qc_multiscale:
        input:
            qc_tsvs=lambda wildcards: expand(BASELINE_MULTI_TSV, bin_label=wildcards.bin_label, sample=BASELINE_SAMPLE_IDS),
            profile_tsvs=lambda wildcards: expand(BASELINE_MULTI_PROFILE_TSV, bin_label=wildcards.bin_label, sample=BASELINE_SAMPLE_IDS),
            metadata=RUN_METADATA
        output:
            summary=BASELINE_MULTI_SUMMARY,
            passed=BASELINE_MULTI_PASS_SAMPLES,
            retained=BASELINE_MULTI_RETAINED_SAMPLES,
            outliers=BASELINE_MULTI_OUTLIER_SAMPLES,
            report=BASELINE_MULTI_REPORT_MD,
            fig_decision=BASELINE_MULTI_FIG_DECISION,
            fig_mapped=BASELINE_MULTI_FIG_MAPPED,
            fig_metrics=BASELINE_MULTI_FIG_METRICS,
            fig_batch=BASELINE_MULTI_FIG_BATCH,
            fig_corr=BASELINE_MULTI_FIG_CORR
        log:
            project_path("logs", "qc", "baseline_multiscale", "{bin_label}", "aggregate.log")
        benchmark:
            BENCH_AGGREGATE_BASELINE_QC_MULTISCALE
        params:
            python_bin=config["biosoft"]["python"],
            script=SCRIPT_AGGREGATE_BASELINE_QC,
            script_action=SCRIPT_AGGREGATE_BASELINE_QC_ACTION,
            figures_dir=lambda wildcards: project_path("qc", "baseline_multiscale", wildcards.bin_label, "figures")
        threads: 1
        shell:
            r"""
            mkdir -p "$(dirname {output.summary})" "$(dirname {log})" "{params.figures_dir}"
            (
                echo "=== PIPELINE AUDIT ==="
                cat {input.metadata:q}
                echo "=== COMMAND ==="
            ) > {log:q}
            {params.python_bin:q} {params.script:q} {params.script_action:q} \
                --qc-tsvs {input.qc_tsvs:q} \
                --profile-tsvs {input.profile_tsvs:q} \
                --summary-output {output.summary:q} \
                --pass-samples-output {output.passed:q} \
                --retained-samples-output {output.retained:q} \
                --outlier-samples-output {output.outliers:q} \
                --report-output {output.report:q} \
                --figures-dir {params.figures_dir:q} \
                --log {log:q} \
                >> {log:q} 2>&1
            """

    rule baseline_sample_samtools_diagnostics:
        input:
            bam=SORTED_BAM,
            bai=SORTED_BAI,
            metadata=RUN_METADATA
        output:
            flagstat=BASELINE_FLAGSTAT,
            idxstats=BASELINE_IDXSTATS
        log:
            project_path("logs", "qc", "baseline_diagnostics", "{sample}.log")
        benchmark:
            BENCH_BASELINE_SAMPLE_SAMTOOLS_DIAGNOSTICS
        params:
            samtools=config["biosoft"]["samtools"]
        threads: BASELINE_QC_SAMTOOLS_THREADS
        shell:
            r"""
            mkdir -p "$(dirname {output.flagstat})" "$(dirname {log})"
            (
                echo "=== PIPELINE AUDIT ==="
                cat {input.metadata:q}
                echo "=== COMMAND ==="
            ) > {log:q}
            {params.samtools:q} flagstat -@ {threads} {input.bam:q} > {output.flagstat:q} 2>> {log:q}
            {params.samtools:q} idxstats {input.bam:q} > {output.idxstats:q} 2>> {log:q}
            """

    rule baseline_qc_review:
        input:
            summaries=expand(BASELINE_MULTI_SUMMARY, bin_label=BASELINE_QC_BIN_LABELS),
            profiles=expand(BASELINE_MULTI_PROFILE_TSV, bin_label=BASELINE_QC_BIN_LABELS, sample=BASELINE_SAMPLE_IDS),
            fastp_jsons=expand(FASTP_JSON, sample=BASELINE_SAMPLE_IDS),
            flagstats=expand(BASELINE_FLAGSTAT, sample=BASELINE_SAMPLE_IDS),
            idxstats=expand(BASELINE_IDXSTATS, sample=BASELINE_SAMPLE_IDS)
        output:
            binsize_summary=BASELINE_MULTI_BINSIZE_SUMMARY,
            sample_matrix=BASELINE_MULTI_SAMPLE_MATRIX,
            stable_retained=BASELINE_MULTI_STABLE_RETAINED,
            failed_diag=BASELINE_MULTI_FAILED_DIAG,
            report=BASELINE_REPORT_MD,
            fig_decision_trend=BASELINE_MULTI_FIG_DECISION_TREND,
            fig_sample_heatmap=BASELINE_MULTI_FIG_SAMPLE_HEATMAP,
            fig_library=BASELINE_MULTI_FIG_LIBRARY,
            fig_fail_chrom=BASELINE_MULTI_FIG_FAIL_CHROM
        log:
            project_path("logs", "qc", "baseline_multiscale", "report.log")
        benchmark:
            BENCH_BASELINE_QC_REVIEW
        params:
            python_bin=config["biosoft"]["python"],
            script=SCRIPT_BASELINE_QC_REPORT,
            script_action=SCRIPT_BASELINE_QC_REPORT_ACTION,
            figures_dir=project_path("qc", "baseline_multiscale", "figures")
        threads: 1
        shell:
            r"""
            mkdir -p "$(dirname {output.binsize_summary})" "$(dirname {output.report})" "$(dirname {log})" "{params.figures_dir}"
            (
                echo "=== COMMAND ==="
            ) > {log:q}
            {params.python_bin:q} {params.script:q} {params.script_action:q} \
                --summary-tsvs {input.summaries:q} \
                --profile-tsvs {input.profiles:q} \
                --fastp-jsons {input.fastp_jsons:q} \
                --flagstats {input.flagstats:q} \
                --idxstats {input.idxstats:q} \
                --binsize-summary-output {output.binsize_summary:q} \
                --sample-matrix-output {output.sample_matrix:q} \
                --stable-retained-output {output.stable_retained:q} \
                --failed-diagnostics-output {output.failed_diag:q} \
                --report-output {output.report:q} \
                --figures-dir {params.figures_dir:q} \
                --log {log:q} \
                >> {log:q} 2>&1
            """
