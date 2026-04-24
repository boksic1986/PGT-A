if QC_FRAMEWORK_TARGET_FILES:
    rule run_structured_qc_framework:
        input:
            summary=BASELINE_QC_SUMMARY,
            profiles=expand(BASELINE_QC_PROFILE_TSV, sample=BASELINE_SAMPLE_IDS),
            annotations=REFERENCE_QC_BIN_ANNOTATIONS,
            combined_mask=REFERENCE_COMBINED_MASK_TSV,
            metadata=RUN_METADATA
        output:
            run_qc_tsv=RUN_QC_TSV,
            run_qc_json=RUN_QC_JSON,
            sample_qc_tsv=SAMPLE_QC_TSV,
            sample_qc_json=SAMPLE_QC_JSON,
            bin_qc_tsv=BIN_QC_TSV,
            bin_qc_json=BIN_QC_JSON,
            event_qc_tsv=EVENT_QC_TSV,
            event_qc_json=EVENT_QC_JSON,
            fig_sample=QC_FRAMEWORK_FIG_SAMPLE,
            fig_event=QC_FRAMEWORK_FIG_EVENT,
            fig_mask=QC_FRAMEWORK_FIG_MASK,
            fig_metric=QC_FRAMEWORK_FIG_METRIC
        log:
            project_path("logs", "qc", "framework", "structured_qc.log")
        params:
            python_bin=config["biosoft"]["python"],
            script=SCRIPT_RUN_QC_FRAMEWORK,
            script_action=SCRIPT_RUN_QC_FRAMEWORK_ACTION
        threads: 1
        shell:
            r"""
            mkdir -p "$(dirname {output.run_qc_tsv})" "$(dirname {log})"
            {params.python_bin:q} {params.script:q} {params.script_action:q} \
                --baseline-summary {input.summary:q} \
                --profile-tsvs {input.profiles:q} \
                --bin-annotations {input.annotations:q} \
                --combined-mask {input.combined_mask:q} \
                --run-metadata {input.metadata:q} \
                --run-qc-tsv {output.run_qc_tsv:q} \
                --run-qc-json {output.run_qc_json:q} \
                --sample-qc-tsv {output.sample_qc_tsv:q} \
                --sample-qc-json {output.sample_qc_json:q} \
                --bin-qc-tsv {output.bin_qc_tsv:q} \
                --bin-qc-json {output.bin_qc_json:q} \
                --event-qc-tsv {output.event_qc_tsv:q} \
                --event-qc-json {output.event_qc_json:q} \
                --sample-status-figure {output.fig_sample:q} \
                --event-status-figure {output.fig_event:q} \
                --mask-status-figure {output.fig_mask:q} \
                --metric-scatter-figure {output.fig_metric:q} \
                >> {log:q} 2>&1
            """
