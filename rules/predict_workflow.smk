if CNV_ENABLED:
    if PREDICT_BY_SEX_ENABLED:
        rule wisecondorx_convert_for_cnv:
            input:
                bam=SORTED_BAM,
                common_binsize=COMMON_REF_BINSIZE,
                metadata=RUN_METADATA
            output:
                npz=CNV_NPZ
            log:
                project_path("logs", "cnv", "{sample}.convert.log")
            benchmark:
                BENCH_WISECONDORX_CONVERT_FOR_CNV
            params:
                wise=config["biosoft"]["WisecondorX"],
                binsize=lambda wildcards, input: read_int_from_file(input.common_binsize)
            threads: 4
            shell:
                r"""
                mkdir -p "$(dirname {output.npz})" "$(dirname {log})"
                (
                    echo "=== PIPELINE AUDIT ==="
                    cat {input.metadata:q}
                    echo "=== COMMON BINSIZE ==="
                    cat {input.common_binsize:q}
                    echo "=== COMMAND ==="
                ) > {log:q}
                {params.wise:q} convert {input.bam:q} {output.npz:q} --binsize {params.binsize} >> {log:q} 2>&1
                """

        rule wisecondorx_gender_for_predict:
            input:
                npz=CNV_NPZ,
                bam=SORTED_BAM,
                ref=GENDER_REF_OUTPUT,
                metadata=RUN_METADATA
            output:
                tsv=CNV_GENDER_TSV
            log:
                project_path("logs", "cnv", "{sample}.gender.log")
            benchmark:
                BENCH_WISECONDORX_GENDER_FOR_PREDICT
            params:
                wise=config["biosoft"]["WisecondorX"],
                samtools=config["biosoft"]["samtools"],
                method=CNV_SEX_CALL_METHOD,
                xx_min_x_relative=CNV_SEX_XX_MIN_X_REL,
                xy_max_x_relative=CNV_SEX_XY_MAX_X_REL,
                xy_min_y_relative=CNV_SEX_XY_MIN_Y_REL,
                xx_max_y_relative=CNV_SEX_XX_MAX_Y_REL
            threads: 1
            run:
                from pgta.core.logging import setup_logger, write_rule_audit_log
                from pgta.predict.sex_routing import run_wisecondorx_gender

                write_rule_audit_log(log[0], input.metadata)
                logger = setup_logger("wisecondorx_gender", log[0])
                run_wisecondorx_gender(
                    wisecondorx=params.wise,
                    sample_npz=input.npz,
                    gender_reference=input.ref,
                    output_tsv=output.tsv,
                    sample_id=wildcards.sample,
                    logger=logger,
                    bam_path=input.bam,
                    samtools=params.samtools,
                    method=params.method,
                    xx_min_x_relative=params.xx_min_x_relative,
                    xy_max_x_relative=params.xy_max_x_relative,
                    xy_min_y_relative=params.xy_min_y_relative,
                    xx_max_y_relative=params.xx_max_y_relative,
                )

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
            benchmark:
                BENCH_WISECONDORX_QC_FOR_PREDICT
            params:
                min_total=CNV_QC_MIN_TOTAL,
                min_nonzero=CNV_QC_MIN_NONZERO,
                max_mad=CNV_QC_MAX_MAD
            threads: 1
            run:
                from pgta.core.logging import setup_logger, write_rule_audit_log
                from pgta.predict.cnv_qc import run_cnv_qc

                write_rule_audit_log(log[0], input.metadata)
                logger = setup_logger("cnv_qc", log[0])
                run_cnv_qc(
                    sample_id=wildcards.sample,
                    npz=input.npz,
                    output_tsv=output.tsv,
                    output_plot=output.plot,
                    pass_marker=output.passed,
                    min_total_counts=params.min_total,
                    min_nonzero_fraction=params.min_nonzero,
                    max_mad_log1p=params.max_mad,
                    logger=logger,
                )

        rule wisecondorx_predict_cnv:
            input:
                npz=CNV_NPZ,
                gender_tsv=CNV_GENDER_TSV,
                qc_report=CNV_QC_TSV,
                qc_pass=CNV_QC_PASS,
                metadata=RUN_METADATA
            output:
                a_branch_bed=CNV_A_ABERRATIONS_BED,
                done=CNV_DONE
            log:
                project_path("logs", "cnv", "{sample}.predict.log")
            benchmark:
                BENCH_WISECONDORX_PREDICT_CNV
            params:
                wise=config["biosoft"]["WisecondorX"],
                ref=lambda wildcards, input: select_predict_reference(input.gender_tsv),
                gender=lambda wildcards, input: select_predict_gender(input.gender_tsv),
                zscore=CNV_ZSCORE,
                alpha=CNV_ALPHA,
                maskrepeats=CNV_MASKREPEATS,
                minrefbins=CNV_MINREFBINS,
                seed=CNV_PREDICT_SEED,
                output_prefix=lambda wildcards: str(Path(CNV_PREDICT_DIR) / wildcards.sample)
            threads: 2
            shell:
                r"""
                mkdir -p "$(dirname {output.done})" "$(dirname {log})"
                export PATH="$(dirname {params.wise}):$PATH"
                (
                    echo "=== PIPELINE AUDIT ==="
                    cat {input.metadata:q}
                    echo "=== GENDER CALL ==="
                    cat {input.gender_tsv:q}
                    echo "=== QC REPORT ==="
                    cat {input.qc_report:q}
                    echo "=== COMMAND ==="
                ) > {log:q}
                {params.wise:q} predict {input.npz:q} {params.ref:q} {params.output_prefix:q} \
                    --gender {params.gender} \
                    --bed \
                    --plot \
                    --zscore {params.zscore} \
                    --alpha {params.alpha} \
                    --maskrepeats {params.maskrepeats} \
                    --minrefbins {params.minrefbins} \
                    --seed {params.seed} >> {log:q} 2>&1
                touch {output.done:q}
                """

        if CNV_POSTPROCESS_ENABLE_BRANCH_B:
            rule cnv_correction_branch_b:
                input:
                    npz=CNV_NPZ,
                    qc_pass=CNV_QC_PASS,
                    annotations=REFERENCE_ANALYSIS_BIN_ANNOTATIONS,
                    combined_mask=REFERENCE_COMBINED_MASK_TSV,
                    metadata=RUN_METADATA
                output:
                    bins=CNV_B_CORRECTED_BINS,
                    summary=CNV_B_CORRECTION_SUMMARY
                log:
                    project_path("logs", "cnv", "{sample}.branch_b.correction.log")
                threads: 1
                run:
                    from pgta.core.logging import write_rule_audit_log
                    import subprocess

                    write_rule_audit_log(log[0], input.metadata, extra_sections=[("QC PASS", open(input.qc_pass, "r", encoding="utf-8").read())])
                    command = [
                        config["biosoft"]["python"],
                        SCRIPT_CNV_CORRECTION,
                        SCRIPT_CNV_CORRECTION_ACTION,
                        "--sample-id", wildcards.sample,
                        "--npz", input.npz,
                        "--annotations", input.annotations,
                        "--combined-mask", input.combined_mask,
                        "--output-bins", output.bins,
                        "--output-summary", output.summary,
                        "--correction-model", CNV_POSTPROCESS_CORRECTION_MODEL,
                        "--loess-frac", str(CNV_CORRECTION_LOESS_FRAC),
                        "--min-valid-bins", str(CNV_CORRECTION_MIN_VALID_BINS),
                        "--robust-iters", str(CNV_CORRECTION_ROBUST_ITERS),
                        "--log", log[0],
                    ]
                    for label in CNV_CORRECTION_INCLUDE_MASK_LABELS:
                        command.extend(["--include-mask-label", label])
                    subprocess.run(command, check=True)

            rule cnv_calling_branch_b:
                input:
                    bins=CNV_B_CORRECTED_BINS,
                    qc_pass=CNV_QC_PASS,
                    metadata=RUN_METADATA
                output:
                    bins=CNV_B_BINS,
                    candidates=CNV_B_CANDIDATES,
                    summary=CNV_B_CALLING_SUMMARY
                log:
                    project_path("logs", "cnv", "{sample}.branch_b.calling.log")
                threads: 1
                run:
                    from pgta.core.logging import write_rule_audit_log
                    import subprocess

                    write_rule_audit_log(log[0], input.metadata, extra_sections=[("QC PASS", open(input.qc_pass, "r", encoding="utf-8").read())])
                    subprocess.run(
                        [
                            config["biosoft"]["python"],
                            SCRIPT_CNV_CALLING,
                            SCRIPT_CNV_CALLING_ACTION,
                            "--sample-id", wildcards.sample,
                            "--input-bins", input.bins,
                            "--output-bins", output.bins,
                            "--output-candidates", output.candidates,
                            "--output-summary", output.summary,
                            "--branch", "B",
                            "--correction-model", CNV_POSTPROCESS_CORRECTION_MODEL,
                            "--min-bins", str(CNV_CALLING_MIN_BINS),
                            "--max-segments-per-chrom", str(CNV_CALLING_MAX_SEGMENTS),
                            "--split-threshold", str(CNV_CALLING_SPLIT_THRESHOLD),
                            "--hmm-state-shift", str(CNV_CALLING_HMM_SHIFT),
                            "--hmm-stay-prob", str(CNV_CALLING_HMM_STAY_PROB),
                            "--min-event-bins", str(CNV_CALLING_MIN_EVENT_BINS),
                            "--min-event-z", str(CNV_CALLING_MIN_EVENT_Z),
                            "--log", log[0],
                        ],
                        check=True,
                    )

            rule cnv_calibration_branch_b:
                input:
                    bins=CNV_B_BINS,
                    candidates=CNV_B_CANDIDATES,
                    metadata=RUN_METADATA
                output:
                    bins=CNV_B_CALIBRATED_BINS,
                    candidates=CNV_B_CALIBRATED_CANDIDATES,
                    summary=CNV_B_CALIBRATION_SUMMARY
                log:
                    project_path("logs", "cnv", "{sample}.branch_b.calibration.log")
                threads: 1
                run:
                    from pgta.core.logging import write_rule_audit_log
                    import subprocess

                    write_rule_audit_log(log[0], input.metadata)
                    subprocess.run(
                        [
                            config["biosoft"]["python"],
                            SCRIPT_CNV_CALIBRATION,
                            SCRIPT_CNV_CALIBRATION_ACTION,
                            "--sample-id", wildcards.sample,
                            "--input-bins", input.bins,
                            "--input-candidates", input.candidates,
                            "--output-bins", output.bins,
                            "--output-candidates", output.candidates,
                            "--output-summary", output.summary,
                            "--null-quantile-low", str(CNV_CAL_NULL_LOW),
                            "--null-quantile-high", str(CNV_CAL_NULL_HIGH),
                            "--min-null-bins", str(CNV_CAL_MIN_NULL_BINS),
                            "--event-z-threshold", str(CNV_CAL_EVENT_Z_THRESHOLD),
                            "--log", log[0],
                        ],
                        check=True,
                    )

            rule cnv_artifact_rules_branch_b:
                input:
                    bins=CNV_B_CALIBRATED_BINS,
                    candidates=CNV_B_CALIBRATED_CANDIDATES,
                    gender_tsv=CNV_GENDER_TSV,
                    metadata=RUN_METADATA
                output:
                    mosaic_candidates=CNV_B_MOSAIC_CANDIDATES,
                    mosaic_summary=CNV_B_MOSAIC_SUMMARY,
                    events=CNV_B_FINAL_EVENTS,
                    summary=CNV_B_ARTIFACT_SUMMARY,
                    json=CNV_B_FINAL_JSON
                log:
                    project_path("logs", "cnv", "{sample}.branch_b.artifact_rules.log")
                threads: 1
                run:
                    from pgta.core.logging import write_rule_audit_log
                    import subprocess

                    write_rule_audit_log(log[0], input.metadata)
                    subprocess.run(
                        [
                            config["biosoft"]["python"],
                            SCRIPT_CNV_MOSAIC_FRACTION,
                            SCRIPT_CNV_MOSAIC_FRACTION_ACTION,
                            "--sample-id", wildcards.sample,
                            "--input-bins", input.bins,
                            "--input-candidates", input.candidates,
                            "--output-candidates", output.mosaic_candidates,
                            "--output-summary", output.mosaic_summary,
                            "--min-effective-bins", str(CNV_MOSAIC_MIN_EFFECTIVE_BINS),
                            "--min-clean-fraction", str(CNV_MOSAIC_MIN_CLEAN_FRACTION),
                            "--max-high-risk-fraction", str(CNV_MOSAIC_MAX_HIGH_RISK_FRACTION),
                            "--min-abs-log2-ratio", str(CNV_MOSAIC_MIN_ABS_LOG2_RATIO),
                            "--low-fraction-threshold", str(CNV_MOSAIC_LOW_FRACTION_THRESHOLD),
                            "--baseline-min-bins", str(CNV_MOSAIC_BASELINE_MIN_BINS),
                            "--ci-zscore", str(CNV_MOSAIC_CI_ZSCORE),
                            "--log", log[0],
                        ],
                        check=True,
                    )
                    command = [
                        config["biosoft"]["python"],
                        SCRIPT_CNV_ARTIFACT_RULES,
                        SCRIPT_CNV_ARTIFACT_RULES_ACTION,
                        "--sample-id", wildcards.sample,
                        "--input-bins", input.bins,
                        "--input-candidates", output.mosaic_candidates,
                        "--gender-tsv", input.gender_tsv,
                        "--output-events", output.events,
                        "--output-summary", output.summary,
                        "--output-json", output.json,
                        "--genome-build", CNV_POSTPROCESS_GENOME_BUILD,
                        "--min-event-bins", str(CNV_ARTIFACT_MIN_BINS),
                        "--min-abs-calibrated-z", str(CNV_ARTIFACT_MIN_ABS_Z),
                        "--max-chrom-fraction", str(CNV_ARTIFACT_MAX_CHROM_FRAC),
                        "--edge-bin-window", str(CNV_ARTIFACT_EDGE_WINDOW),
                        "--max-qvalue", str(CNV_ARTIFACT_MAX_QVALUE),
                        "--keep-review", str(CNV_ARTIFACT_KEEP_REVIEW),
                        "--high-confidence-z", str(CNV_ARTIFACT_HIGH_CONF_Z),
                        "--high-confidence-qvalue", str(CNV_ARTIFACT_HIGH_CONF_QVALUE),
                        "--log", log[0],
                    ]
                    for region in CNV_POSTPROCESS_PAR_REGIONS:
                        command.extend(["--par-region", region])
                    subprocess.run(command, check=True)

    else:
        rule wisecondorx_convert_for_cnv:
            input:
                bam=SORTED_BAM,
                metadata=RUN_METADATA
            output:
                npz=CNV_NPZ
            log:
                project_path("logs", "cnv", "{sample}.convert.log")
            benchmark:
                BENCH_WISECONDORX_CONVERT_FOR_CNV
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
            benchmark:
                BENCH_WISECONDORX_QC_FOR_PREDICT
            params:
                min_total=CNV_QC_MIN_TOTAL,
                min_nonzero=CNV_QC_MIN_NONZERO,
                max_mad=CNV_QC_MAX_MAD
            threads: 1
            run:
                from pgta.core.logging import setup_logger, write_rule_audit_log
                from pgta.predict.cnv_qc import run_cnv_qc

                write_rule_audit_log(log[0], input.metadata)
                logger = setup_logger("cnv_qc", log[0])
                run_cnv_qc(
                    sample_id=wildcards.sample,
                    npz=input.npz,
                    output_tsv=output.tsv,
                    output_plot=output.plot,
                    pass_marker=output.passed,
                    min_total_counts=params.min_total,
                    min_nonzero_fraction=params.min_nonzero,
                    max_mad_log1p=params.max_mad,
                    logger=logger,
                )

        rule wisecondorx_predict_cnv:
            input:
                npz=CNV_NPZ,
                ref=REF_OUTPUT,
                qc_report=CNV_QC_TSV,
                qc_pass=CNV_QC_PASS,
                metadata=RUN_METADATA
            output:
                a_branch_bed=CNV_A_ABERRATIONS_BED,
                done=CNV_DONE
            log:
                project_path("logs", "cnv", "{sample}.predict.log")
            benchmark:
                BENCH_WISECONDORX_PREDICT_CNV
            params:
                wise=config["biosoft"]["WisecondorX"],
                zscore=CNV_ZSCORE,
                alpha=CNV_ALPHA,
                maskrepeats=CNV_MASKREPEATS,
                minrefbins=CNV_MINREFBINS,
                seed=CNV_PREDICT_SEED,
                output_prefix=lambda wildcards: str(Path(CNV_PREDICT_DIR) / wildcards.sample)
            threads: 2
            shell:
                r"""
                mkdir -p "$(dirname {output.done})" "$(dirname {log})"
                export PATH="$(dirname {params.wise}):$PATH"
                (
                    echo "=== PIPELINE AUDIT ==="
                    cat {input.metadata:q}
                    echo "=== QC REPORT ==="
                    cat {input.qc_report:q}
                    echo "=== COMMAND ==="
                ) > {log:q}
                {params.wise:q} predict {input.npz:q} {input.ref:q} {params.output_prefix:q} \
                    --bed \
                    --plot \
                    --zscore {params.zscore} \
                    --alpha {params.alpha} \
                    --maskrepeats {params.maskrepeats} \
                    --minrefbins {params.minrefbins} \
                    --seed {params.seed} >> {log:q} 2>&1
                touch {output.done:q}
                """

        if CNV_POSTPROCESS_ENABLE_BRANCH_B:
            rule cnv_correction_branch_b:
                input:
                    npz=CNV_NPZ,
                    qc_pass=CNV_QC_PASS,
                    annotations=REFERENCE_ANALYSIS_BIN_ANNOTATIONS,
                    combined_mask=REFERENCE_COMBINED_MASK_TSV,
                    metadata=RUN_METADATA
                output:
                    bins=CNV_B_CORRECTED_BINS,
                    summary=CNV_B_CORRECTION_SUMMARY
                log:
                    project_path("logs", "cnv", "{sample}.branch_b.correction.log")
                threads: 1
                run:
                    from pgta.core.logging import write_rule_audit_log
                    import subprocess

                    write_rule_audit_log(log[0], input.metadata, extra_sections=[("QC PASS", open(input.qc_pass, "r", encoding="utf-8").read())])
                    command = [
                        config["biosoft"]["python"],
                        SCRIPT_CNV_CORRECTION,
                        SCRIPT_CNV_CORRECTION_ACTION,
                        "--sample-id", wildcards.sample,
                        "--npz", input.npz,
                        "--annotations", input.annotations,
                        "--combined-mask", input.combined_mask,
                        "--output-bins", output.bins,
                        "--output-summary", output.summary,
                        "--correction-model", CNV_POSTPROCESS_CORRECTION_MODEL,
                        "--loess-frac", str(CNV_CORRECTION_LOESS_FRAC),
                        "--min-valid-bins", str(CNV_CORRECTION_MIN_VALID_BINS),
                        "--robust-iters", str(CNV_CORRECTION_ROBUST_ITERS),
                        "--log", log[0],
                    ]
                    for label in CNV_CORRECTION_INCLUDE_MASK_LABELS:
                        command.extend(["--include-mask-label", label])
                    subprocess.run(command, check=True)

            rule cnv_calling_branch_b:
                input:
                    bins=CNV_B_CORRECTED_BINS,
                    qc_pass=CNV_QC_PASS,
                    metadata=RUN_METADATA
                output:
                    bins=CNV_B_BINS,
                    candidates=CNV_B_CANDIDATES,
                    summary=CNV_B_CALLING_SUMMARY
                log:
                    project_path("logs", "cnv", "{sample}.branch_b.calling.log")
                threads: 1
                run:
                    from pgta.core.logging import write_rule_audit_log
                    import subprocess

                    write_rule_audit_log(log[0], input.metadata, extra_sections=[("QC PASS", open(input.qc_pass, "r", encoding="utf-8").read())])
                    subprocess.run(
                        [
                            config["biosoft"]["python"],
                            SCRIPT_CNV_CALLING,
                            SCRIPT_CNV_CALLING_ACTION,
                            "--sample-id", wildcards.sample,
                            "--input-bins", input.bins,
                            "--output-bins", output.bins,
                            "--output-candidates", output.candidates,
                            "--output-summary", output.summary,
                            "--branch", "B",
                            "--correction-model", CNV_POSTPROCESS_CORRECTION_MODEL,
                            "--min-bins", str(CNV_CALLING_MIN_BINS),
                            "--max-segments-per-chrom", str(CNV_CALLING_MAX_SEGMENTS),
                            "--split-threshold", str(CNV_CALLING_SPLIT_THRESHOLD),
                            "--hmm-state-shift", str(CNV_CALLING_HMM_SHIFT),
                            "--hmm-stay-prob", str(CNV_CALLING_HMM_STAY_PROB),
                            "--min-event-bins", str(CNV_CALLING_MIN_EVENT_BINS),
                            "--min-event-z", str(CNV_CALLING_MIN_EVENT_Z),
                            "--log", log[0],
                        ],
                        check=True,
                    )

            rule cnv_calibration_branch_b:
                input:
                    bins=CNV_B_BINS,
                    candidates=CNV_B_CANDIDATES,
                    metadata=RUN_METADATA
                output:
                    bins=CNV_B_CALIBRATED_BINS,
                    candidates=CNV_B_CALIBRATED_CANDIDATES,
                    summary=CNV_B_CALIBRATION_SUMMARY
                log:
                    project_path("logs", "cnv", "{sample}.branch_b.calibration.log")
                threads: 1
                run:
                    from pgta.core.logging import write_rule_audit_log
                    import subprocess

                    write_rule_audit_log(log[0], input.metadata)
                    subprocess.run(
                        [
                            config["biosoft"]["python"],
                            SCRIPT_CNV_CALIBRATION,
                            SCRIPT_CNV_CALIBRATION_ACTION,
                            "--sample-id", wildcards.sample,
                            "--input-bins", input.bins,
                            "--input-candidates", input.candidates,
                            "--output-bins", output.bins,
                            "--output-candidates", output.candidates,
                            "--output-summary", output.summary,
                            "--null-quantile-low", str(CNV_CAL_NULL_LOW),
                            "--null-quantile-high", str(CNV_CAL_NULL_HIGH),
                            "--min-null-bins", str(CNV_CAL_MIN_NULL_BINS),
                            "--event-z-threshold", str(CNV_CAL_EVENT_Z_THRESHOLD),
                            "--log", log[0],
                        ],
                        check=True,
                    )

            rule cnv_artifact_rules_branch_b:
                input:
                    bins=CNV_B_CALIBRATED_BINS,
                    candidates=CNV_B_CALIBRATED_CANDIDATES,
                    metadata=RUN_METADATA
                output:
                    mosaic_candidates=CNV_B_MOSAIC_CANDIDATES,
                    mosaic_summary=CNV_B_MOSAIC_SUMMARY,
                    events=CNV_B_FINAL_EVENTS,
                    summary=CNV_B_ARTIFACT_SUMMARY,
                    json=CNV_B_FINAL_JSON
                log:
                    project_path("logs", "cnv", "{sample}.branch_b.artifact_rules.log")
                threads: 1
                run:
                    from pgta.core.logging import write_rule_audit_log
                    import subprocess

                    write_rule_audit_log(log[0], input.metadata)
                    subprocess.run(
                        [
                            config["biosoft"]["python"],
                            SCRIPT_CNV_MOSAIC_FRACTION,
                            SCRIPT_CNV_MOSAIC_FRACTION_ACTION,
                            "--sample-id", wildcards.sample,
                            "--input-bins", input.bins,
                            "--input-candidates", input.candidates,
                            "--output-candidates", output.mosaic_candidates,
                            "--output-summary", output.mosaic_summary,
                            "--min-effective-bins", str(CNV_MOSAIC_MIN_EFFECTIVE_BINS),
                            "--min-clean-fraction", str(CNV_MOSAIC_MIN_CLEAN_FRACTION),
                            "--max-high-risk-fraction", str(CNV_MOSAIC_MAX_HIGH_RISK_FRACTION),
                            "--min-abs-log2-ratio", str(CNV_MOSAIC_MIN_ABS_LOG2_RATIO),
                            "--low-fraction-threshold", str(CNV_MOSAIC_LOW_FRACTION_THRESHOLD),
                            "--baseline-min-bins", str(CNV_MOSAIC_BASELINE_MIN_BINS),
                            "--ci-zscore", str(CNV_MOSAIC_CI_ZSCORE),
                            "--log", log[0],
                        ],
                        check=True,
                    )
                    command = [
                        config["biosoft"]["python"],
                        SCRIPT_CNV_ARTIFACT_RULES,
                        SCRIPT_CNV_ARTIFACT_RULES_ACTION,
                        "--sample-id", wildcards.sample,
                        "--input-bins", input.bins,
                        "--input-candidates", output.mosaic_candidates,
                        "--output-events", output.events,
                        "--output-summary", output.summary,
                        "--output-json", output.json,
                        "--genome-build", CNV_POSTPROCESS_GENOME_BUILD,
                        "--min-event-bins", str(CNV_ARTIFACT_MIN_BINS),
                        "--min-abs-calibrated-z", str(CNV_ARTIFACT_MIN_ABS_Z),
                        "--max-chrom-fraction", str(CNV_ARTIFACT_MAX_CHROM_FRAC),
                        "--edge-bin-window", str(CNV_ARTIFACT_EDGE_WINDOW),
                        "--max-qvalue", str(CNV_ARTIFACT_MAX_QVALUE),
                        "--keep-review", str(CNV_ARTIFACT_KEEP_REVIEW),
                        "--high-confidence-z", str(CNV_ARTIFACT_HIGH_CONF_Z),
                        "--high-confidence-qvalue", str(CNV_ARTIFACT_HIGH_CONF_QVALUE),
                        "--log", log[0],
                    ]
                    for region in CNV_POSTPROCESS_PAR_REGIONS:
                        command.extend(["--par-region", region])
                    subprocess.run(command, check=True)

    if CNV_MOSAIC_FRACTION_TRUTH_TSV and ("cnv_benchmark" in AVAILABLE_TARGETS or "cnv_report" in AVAILABLE_TARGETS):
        rule cnv_mosaic_truth_validation:
            input:
                metadata=RUN_METADATA
            output:
                summary=CNV_MOSAIC_TRUTH_VALIDATION_SUMMARY
            log:
                project_path("logs", "cnv", "mosaic_truth_validation.log")
            threads: 1
            run:
                from pgta.core.logging import write_rule_audit_log
                import subprocess

                write_rule_audit_log(log[0], input.metadata)
                subprocess.run(
                    [
                        config["biosoft"]["python"],
                        SCRIPT_VALIDATE_FRACTION_TRUTH,
                        SCRIPT_VALIDATE_FRACTION_TRUTH_ACTION,
                        "--input-tsv", CNV_MOSAIC_FRACTION_TRUTH_TSV,
                        "--output-json", output.summary,
                    ],
                    check=True,
                )

    if "cnv_eval" in AVAILABLE_TARGETS:
        rule cnv_evaluation:
            input:
                metadata=RUN_METADATA,
                events=(expand(CNV_B_FINAL_EVENTS, sample=SAMPLES) if CNV_POSTPROCESS_ENABLE_BRANCH_B else []),
                genders=(expand(CNV_GENDER_TSV, sample=SAMPLES) if PREDICT_BY_SEX_ENABLED else []),
                qcs=expand(CNV_QC_TSV, sample=SAMPLES)
            output:
                sample_metrics=CNV_EVAL_SAMPLE_METRICS,
                event_metrics=CNV_EVAL_EVENT_METRICS,
                calibration=CNV_EVAL_CALIBRATION,
                summary=CNV_EVAL_SUMMARY
            log:
                project_path("logs", "cnv", "evaluation.log")
            threads: 1
            run:
                from pgta.core.logging import write_rule_audit_log
                import subprocess

                write_rule_audit_log(log[0], input.metadata)
                command = [
                    config["biosoft"]["python"],
                    SCRIPT_CNV_EVALUATION,
                    SCRIPT_CNV_EVALUATION_ACTION,
                    "--output-sample-metrics", output.sample_metrics,
                    "--output-event-metrics", output.event_metrics,
                    "--output-calibration", output.calibration,
                    "--output-summary", output.summary,
                    "--truth-tsv", CNV_EVAL_TRUTH_TSV,
                    "--log", log[0],
                ]
                for path_value in input.events:
                    command.extend(["--event-tsv", path_value])
                for path_value in input.genders:
                    command.extend(["--gender-tsv", path_value])
                for path_value in input.qcs:
                    command.extend(["--qc-tsv", path_value])
                subprocess.run(command, check=True)

    if "cnv_ml" in AVAILABLE_TARGETS:
        rule cnv_ml_candidate_classifier:
            input:
                metadata=RUN_METADATA,
                events=(expand(CNV_B_FINAL_EVENTS, sample=SAMPLES) if CNV_POSTPROCESS_ENABLE_BRANCH_B else [])
            output:
                features=CNV_ML_FEATURES,
                cv_metrics=CNV_ML_CV_METRICS,
                calibration=CNV_ML_CALIBRATION,
                importance=CNV_ML_IMPORTANCE,
                predictions=CNV_ML_PREDICTIONS,
                summary=CNV_ML_SUMMARY
            log:
                project_path("logs", "cnv", "ml.log")
            threads: 1
            run:
                from pgta.core.logging import write_rule_audit_log
                import subprocess

                write_rule_audit_log(log[0], input.metadata)
                command = [
                    config["biosoft"]["python"],
                    SCRIPT_CNV_ML,
                    SCRIPT_CNV_ML_ACTION,
                    "--output-features", output.features,
                    "--output-cv-metrics", output.cv_metrics,
                    "--output-calibration", output.calibration,
                    "--output-importance", output.importance,
                    "--output-predictions", output.predictions,
                    "--output-summary", output.summary,
                    "--backend", CNV_ML_BACKEND,
                    "--cv-folds", str(CNV_ML_CV_FOLDS),
                    "--labels-tsv", CNV_ML_LABELS_TSV,
                    "--log", log[0],
                ]
                for path_value in input.events:
                    command.extend(["--event-tsv", path_value])
                subprocess.run(command, check=True)

    if "cnv_benchmark" in AVAILABLE_TARGETS:
        rule cnv_benchmark_framework:
            input:
                metadata=RUN_METADATA,
                events=(expand(CNV_B_FINAL_EVENTS, sample=SAMPLES) if CNV_POSTPROCESS_ENABLE_BRANCH_B else []),
                a_branch=(expand(CNV_A_ABERRATIONS_BED, sample=SAMPLES) if CNV_POSTPROCESS_PRESERVE_BRANCH_A else []),
                mosaic_truth_validation=([CNV_MOSAIC_TRUTH_VALIDATION_SUMMARY] if CNV_MOSAIC_FRACTION_TRUTH_TSV else [])
            output:
                simulation=CNV_BENCHMARK_SIMULATION,
                admixture=CNV_BENCHMARK_ADMIXTURE,
                summary=CNV_BENCHMARK_SUMMARY
            log:
                project_path("logs", "cnv", "benchmark.log")
            threads: 1
            run:
                from pgta.core.logging import write_rule_audit_log
                import subprocess

                write_rule_audit_log(log[0], input.metadata)
                command = [
                    config["biosoft"]["python"],
                    SCRIPT_CNV_BENCHMARK,
                    SCRIPT_CNV_BENCHMARK_ACTION,
                    "--output-simulation", output.simulation,
                    "--output-admixture", output.admixture,
                    "--output-summary", output.summary,
                    "--truth-tsv", CNV_BENCHMARK_TRUTH_TSV,
                    "--branch-b-z-threshold", str(CNV_BENCHMARK_BRANCH_B_Z_THRESHOLD),
                    "--branch-a-z-threshold", str(CNV_BENCHMARK_BRANCH_A_Z_THRESHOLD),
                    "--log", log[0],
                ]
                for path_value in input.events:
                    command.extend(["--event-tsv", path_value])
                for path_value in input.a_branch:
                    command.extend(["--a-branch-bed", path_value])
                for admixture_level in CNV_BENCHMARK_ADMIXTURE_LEVELS:
                    command.extend(["--admixture-level", str(admixture_level)])
                for threshold in CNV_BENCHMARK_LOW_FRACTION_THRESHOLDS:
                    command.extend(["--low-fraction-threshold", str(threshold)])
                subprocess.run(command, check=True)

    if "cnv_report" in AVAILABLE_TARGETS:
        rule cnv_report_summary:
            input:
                metadata=RUN_METADATA,
                events=(expand(CNV_B_FINAL_EVENTS, sample=SAMPLES) if CNV_POSTPROCESS_ENABLE_BRANCH_B else []),
                genders=(expand(CNV_GENDER_TSV, sample=SAMPLES) if PREDICT_BY_SEX_ENABLED else []),
                qcs=expand(CNV_QC_TSV, sample=SAMPLES),
                a_branch=(expand(CNV_A_ABERRATIONS_BED, sample=SAMPLES) if CNV_POSTPROCESS_PRESERVE_BRANCH_A else []),
                evaluation_summary=([CNV_EVAL_SUMMARY] if "cnv_eval" in REQUESTED_TARGETS else []),
                ml_summary=([CNV_ML_SUMMARY] if "cnv_ml" in REQUESTED_TARGETS else []),
                benchmark_summary=([CNV_BENCHMARK_SUMMARY] if "cnv_benchmark" in REQUESTED_TARGETS else []),
                mosaic_truth_validation=([CNV_MOSAIC_TRUTH_VALIDATION_SUMMARY] if CNV_MOSAIC_FRACTION_TRUTH_TSV else [])
            output:
                tsv=CNV_REPORT_TSV,
                json=CNV_REPORT_JSON,
                md=CNV_REPORT_MD,
                html=CNV_REPORT_HTML
            log:
                project_path("logs", "cnv", "report.log")
            threads: 1
            run:
                from pgta.core.logging import write_rule_audit_log
                import subprocess

                write_rule_audit_log(log[0], input.metadata)
                command = [
                    config["biosoft"]["python"],
                    SCRIPT_CNV_REPORT,
                    SCRIPT_CNV_REPORT_ACTION,
                    "--output-tsv", output.tsv,
                    "--output-json", output.json,
                    "--output-md", output.md,
                    "--output-html", output.html,
                    "--evaluation-summary", CNV_EVAL_SUMMARY,
                    "--ml-summary", CNV_ML_SUMMARY,
                    "--benchmark-summary", CNV_BENCHMARK_SUMMARY,
                    "--log", log[0],
                ]
                if input.mosaic_truth_validation:
                    command.extend(["--truth-validation-summary", input.mosaic_truth_validation[0]])
                for path_value in input.events:
                    command.extend(["--event-tsv", path_value])
                for path_value in input.genders:
                    command.extend(["--gender-tsv", path_value])
                for path_value in input.qcs:
                    command.extend(["--qc-tsv", path_value])
                for path_value in input.a_branch:
                    command.extend(["--a-branch-bed", path_value])
                subprocess.run(command, check=True)
