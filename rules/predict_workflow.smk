if CNV_ENABLED:
    pass
    rule wisecondorx_convert_for_cnv:
        input:
            bam=SORTED_BAM,
            metadata=RUN_METADATA
        output:
            npz=CNV_NPZ
        log:
            project_path("logs", "cnv", "{sample}.convert.log")
        params:
            wise=config["biosoft"]["WisecondorX"],
            binsize=CNV_CONVERT_BINSIZE
        threads: 4
        shell:
            r"""
            mkdir -p "$(dirname {output.npz})" "$(dirname {log})"
            (
                echo "=== PIPELINE AUDIT ==="
                cat {input.metadata:q}
                echo "=== COMMAND ==="
            ) > {log:q}
            {params.wise:q} convert {input.bam:q} {output.npz:q} --binsize {params.binsize} >> {log:q} 2>&1
            """

    rule wisecondorx_qc_for_predict:
        input:
            npz=CNV_NPZ,
            metadata=RUN_METADATA
        output:
            tsv=CNV_QC_TSV,
            plot=CNV_QC_PLOT,
            passed=CNV_QC_PASS
        log:
            project_path("logs", "cnv", "{sample}.qc.log")
        params:
            python_bin=config["biosoft"]["python"],
            min_total=CNV_QC_MIN_TOTAL,
            min_nonzero=CNV_QC_MIN_NONZERO,
            max_mad=CNV_QC_MAX_MAD
        threads: 1
        shell:
            r"""
            mkdir -p "$(dirname {output.tsv})" "$(dirname {log})"
            (
                echo "=== PIPELINE AUDIT ==="
                cat {input.metadata:q}
                echo "=== COMMAND ==="
            ) > {log:q}
            {params.python_bin:q} {SCRIPT_CNV_QC:q} \
                --sample-id {wildcards.sample:q} \
                --npz {input.npz:q} \
                --output-tsv {output.tsv:q} \
                --output-plot {output.plot:q} \
                --pass-marker {output.passed:q} \
                --min-total-counts {params.min_total} \
                --min-nonzero-fraction {params.min_nonzero} \
                --max-mad-log1p {params.max_mad} \
                --log {log:q} \
                >> {log:q} 2>&1
            """

    rule wisecondorx_predict_cnv:
        input:
            npz=CNV_NPZ,
            ref=REF_OUTPUT,
            qc_report=CNV_QC_TSV,
            qc_pass=CNV_QC_PASS,
            metadata=RUN_METADATA
        output:
            done=CNV_DONE
        log:
            project_path("logs", "cnv", "{sample}.predict.log")
        params:
            wise=config["biosoft"]["WisecondorX"],
            zscore=CNV_ZSCORE,
            alpha=CNV_ALPHA,
            maskrepeats=CNV_MASKREPEATS,
            minrefbins=CNV_MINREFBINS,
            output_prefix=lambda wildcards: str(Path(CNV_PREDICT_DIR) / wildcards.sample)
        threads: 2
        shell:
            r"""
            mkdir -p "$(dirname {output.done})" "$(dirname {log})"
            (
                echo "=== PIPELINE AUDIT ==="
                cat {input.metadata:q}
                echo "=== QC REPORT ==="
                cat {input.qc_report:q}
                echo "=== COMMAND ==="
            ) > {log:q}
            {params.wise:q} predict {input.npz:q} {input.ref:q} {params.output_prefix:q} \
                --zscore {params.zscore} \
                --alpha {params.alpha} \
                --maskrepeats {params.maskrepeats} \
                --minrefbins {params.minrefbins} \
                --cpus {threads} >> {log:q} 2>&1
            touch {output.done:q}
            """
