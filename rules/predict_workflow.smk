if CNV_ENABLED:
    if PREDICT_BY_SEX_ENABLED:
        rule wisecondorx_convert_for_cnv:
            input:
                bam=lambda wildcards: sample_bam_path(wildcards.sample),
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
                bam=lambda wildcards: sample_bam_path(wildcards.sample),
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
                npz=CNV_PREDICT_INPUT_NPZ,
                gender_tsv=CNV_GENDER_TSV,
                qc_report=CNV_QC_TSV,
                qc_pass=CNV_QC_PASS,
                blacklist=([CNV_POSTPROCESS_BLACKLIST_BED] if CNV_POSTPROCESS_BLACKLIST_BED else []),
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
                blacklist=CNV_POSTPROCESS_BLACKLIST_BED,
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
                blacklist_args=""
                if [ -n {params.blacklist:q} ]; then
                    if [ -s {params.blacklist:q} ]; then
                        blacklist_args="--blacklist {params.blacklist}"
                    else
                        echo "blacklist not passed because file missing or empty: {params.blacklist}" >> {log:q}
                    fi
                else
                    echo "blacklist not passed because path is empty" >> {log:q}
                fi
                {params.wise:q} predict {input.npz:q} {params.ref:q} {params.output_prefix:q} \
                    --gender {params.gender} \
                    --bed \
                    --plot \
                    --zscore {params.zscore} \
                    --alpha {params.alpha} \
                    --maskrepeats {params.maskrepeats} \
                    --minrefbins {params.minrefbins} \
                    --seed {params.seed} $blacklist_args >> {log:q} 2>&1
                touch {output.done:q}
                """

    else:
        rule wisecondorx_convert_for_cnv:
            input:
                bam=lambda wildcards: sample_bam_path(wildcards.sample),
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
                npz=CNV_PREDICT_INPUT_NPZ,
                ref=REF_OUTPUT,
                qc_report=CNV_QC_TSV,
                qc_pass=CNV_QC_PASS,
                blacklist=([CNV_POSTPROCESS_BLACKLIST_BED] if CNV_POSTPROCESS_BLACKLIST_BED else []),
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
                blacklist=CNV_POSTPROCESS_BLACKLIST_BED,
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
                blacklist_args=""
                if [ -n {params.blacklist:q} ]; then
                    if [ -s {params.blacklist:q} ]; then
                        blacklist_args="--blacklist {params.blacklist}"
                    else
                        echo "blacklist not passed because file missing or empty: {params.blacklist}" >> {log:q}
                    fi
                else
                    echo "blacklist not passed because path is empty" >> {log:q}
                fi
                {params.wise:q} predict {input.npz:q} {input.ref:q} {params.output_prefix:q} \
                    --bed \
                    --plot \
                    --zscore {params.zscore} \
                    --alpha {params.alpha} \
                    --maskrepeats {params.maskrepeats} \
                    --minrefbins {params.minrefbins} \
                    --seed {params.seed} $blacklist_args >> {log:q} 2>&1
                touch {output.done:q}
                """

    if CNV_PREPROCESS_STRATEGY == "mask_only":
        rule cnv_mask_npz_for_predict:
            input:
                npz=CNV_NPZ,
                combined_mask=REFERENCE_COMBINED_MASK_TSV,
                metadata=RUN_METADATA
            output:
                npz=CNV_MASKED_NPZ,
                summary=CNV_MASK_SUMMARY
            log:
                project_path("logs", "cnv", "{sample}.mask_npz.log")
            params:
                python_bin=config["biosoft"]["python"]
            threads: 1
            shell:
                r"""
                mkdir -p "$(dirname {output.npz})" "$(dirname {log})"
                (
                    echo "=== PIPELINE AUDIT ==="
                    cat {input.metadata:q}
                    echo "=== COMMAND ==="
                ) > {log:q}
                {params.python_bin:q} -m pgta.reference.preprocess mask_npz \
                    --input-npz {input.npz:q} \
                    --combined-mask {input.combined_mask:q} \
                    --output-npz {output.npz:q} \
                    --output-summary {output.summary:q} \
                    >> {log:q} 2>&1
                """

    rule a_branch_candidate_assembly:
        input:
            bed=CNV_A_ABERRATIONS_BED,
            done=CNV_DONE,
            metadata=RUN_METADATA
        output:
            candidates=CNV_A_CANDIDATES,
            summary=CNV_A_CANDIDATE_SUMMARY
        log:
            str(Path(CNV_BRANCH_A_LOG_DIR) / "{sample}.a_branch_candidates.log")
        params:
            merge_gap_bp=CNV_BRANCH_A_MERGE_GAP_BP,
            strong_z=CNV_BRANCH_A_STRONG_Z
        threads: 1
        run:
            from pgta.core.logging import write_rule_audit_log
            import subprocess

            write_rule_audit_log(log[0], input.metadata)
            subprocess.run(
                [
                    config["biosoft"]["python"],
                    SCRIPT_A_BRANCH_CANDIDATES,
                    SCRIPT_A_BRANCH_CANDIDATES_ACTION,
                    "--sample-id", wildcards.sample,
                    "--input-bed", input.bed,
                    "--output-candidates", output.candidates,
                    "--output-summary", output.summary,
                    "--merge-gap-bp", str(params.merge_gap_bp),
                    "--strong-z", str(params.strong_z),
                    "--log", log[0],
                ],
                check=True,
            )

    rule cnv_branch_a_validation:
        input:
            a_candidates=expand(CNV_A_CANDIDATES, sample=SAMPLES),
            qcs=expand(CNV_QC_TSV, sample=SAMPLES),
            genders=(expand(CNV_GENDER_TSV, sample=SAMPLES) if PREDICT_BY_SEX_ENABLED else []),
            truth=([CNV_EVAL_TRUTH_TSV] if CNV_EVAL_TRUTH_TSV else []),
            metadata=RUN_METADATA
        output:
            sample_summary=CNV_BRANCH_A_VALIDATION_SAMPLE_SUMMARY,
            truth_metrics=CNV_BRANCH_A_VALIDATION_TRUTH_METRICS,
            summary=CNV_BRANCH_A_VALIDATION_SUMMARY
        log:
            str(Path(CNV_BRANCH_A_LOG_DIR) / "branch_a_validation.log")
        threads: 1
        run:
            from pgta.core.logging import write_rule_audit_log
            import subprocess

            write_rule_audit_log(log[0], input.metadata)
            command = [
                config["biosoft"]["python"],
                SCRIPT_BRANCH_A_VALIDATION,
                SCRIPT_BRANCH_A_VALIDATION_ACTION,
                "--reference-id", CNV_REFERENCE_ID,
                "--binsize", str(CNV_CONVERT_BINSIZE),
                "--preprocess-mask-version", CNV_PREPROCESS_STRATEGY,
                "--wisecondorx-predict-command",
                (
                    "WisecondorX predict <npz> <sex-routed-ref> <output-prefix> "
                    f"--zscore {CNV_ZSCORE} --alpha {CNV_ALPHA} "
                    f"--maskrepeats {CNV_MASKREPEATS} --minrefbins {CNV_MINREFBINS}"
                ),
                "--blacklist-bed", CNV_POSTPROCESS_BLACKLIST_BED,
                "--minrefbins", str(CNV_MINREFBINS),
                "--zscore", str(CNV_ZSCORE),
                "--alpha", str(CNV_ALPHA),
                "--maskrepeats", str(CNV_MASKREPEATS),
                "--output-sample-summary", output.sample_summary,
                "--output-truth-metrics", output.truth_metrics,
                "--output-summary", output.summary,
            ]
            for sample_id in SAMPLES:
                command.extend(["--sample-id", str(sample_id)])
            for candidate_tsv in input.a_candidates:
                command.extend(["--candidate-tsv", str(candidate_tsv)])
            for qc_tsv in input.qcs:
                command.extend(["--qc-tsv", str(qc_tsv)])
            for gender_tsv in input.genders:
                command.extend(["--gender-tsv", str(gender_tsv)])
            if input.truth:
                command.extend(["--truth-tsv", str(input.truth[0])])
            subprocess.run(command, check=True)

    if "reference_audit" in AVAILABLE_TARGETS:
        rule cnv_reference_candidate_audit:
            input:
                a_candidates=expand(CNV_A_CANDIDATES, sample=CNV_REFERENCE_AUDIT_SAMPLE_IDS),
                qcs=expand(CNV_QC_TSV, sample=CNV_REFERENCE_AUDIT_SAMPLE_IDS),
                genders=(expand(CNV_GENDER_TSV, sample=CNV_REFERENCE_AUDIT_SAMPLE_IDS) if PREDICT_BY_SEX_ENABLED else []),
                sample_metadata=([CNV_REFERENCE_AUDIT_SAMPLE_METADATA] if CNV_REFERENCE_AUDIT_SAMPLE_METADATA else []),
                reference_samples=([CNV_REFERENCE_AUDIT_REFERENCE_SAMPLES_FILE] if CNV_REFERENCE_AUDIT_REFERENCE_SAMPLES_FILE else []),
                bin_annotations=([CNV_REFERENCE_AUDIT_BIN_ANNOTATIONS] if CNV_REFERENCE_AUDIT_BIN_ANNOTATIONS else []),
                metadata=RUN_METADATA
            output:
                audit=CNV_REFERENCE_AUDIT_TSV,
                summary=CNV_REFERENCE_AUDIT_SUMMARY
            log:
                project_path("logs", "cnv", "reference_candidate_audit.log")
            threads: 1
            run:
                from pgta.core.logging import write_rule_audit_log
                import subprocess

                write_rule_audit_log(log[0], input.metadata)
                command = [
                    config["biosoft"]["python"],
                    SCRIPT_REFERENCE_CANDIDATE_AUDIT,
                    SCRIPT_REFERENCE_CANDIDATE_AUDIT_ACTION,
                    "--a-candidates-dir", CNV_BRANCH_A_OUTPUT_DIR,
                    "--qc-dir", CNV_QC_DIR,
                    "--reference-id", CNV_REFERENCE_ID,
                    "--output-audit", output.audit,
                    "--output-summary", output.summary,
                    "--strong-z", str(CNV_REFERENCE_AUDIT_STRONG_Z),
                    "--broad-event-min-bp", str(CNV_REFERENCE_AUDIT_BROAD_EVENT_MIN_BP),
                    "--sample-specific-fraction-threshold", str(CNV_REFERENCE_AUDIT_SAMPLE_SPECIFIC_FRACTION_THRESHOLD),
                    "--high-risk-burden-threshold", str(CNV_REFERENCE_AUDIT_HIGH_RISK_BURDEN_THRESHOLD),
                    "--shared-signal-min-samples", str(CNV_REFERENCE_AUDIT_SHARED_SIGNAL_MIN_SAMPLES),
                ]
                if PREDICT_BY_SEX_ENABLED:
                    command.extend(["--gender-dir", CNV_GENDER_DIR])
                if CNV_REFERENCE_AUDIT_EVIDENCE_LEDGER_DIR:
                    command.extend(["--evidence-ledger-dir", CNV_REFERENCE_AUDIT_EVIDENCE_LEDGER_DIR])
                if input.sample_metadata:
                    command.extend(["--sample-metadata", input.sample_metadata[0]])
                else:
                    for sample_id in CNV_REFERENCE_AUDIT_SAMPLE_IDS:
                        command.extend(["--sample-id", str(sample_id)])
                if input.reference_samples:
                    command.extend(["--reference-samples-file", input.reference_samples[0]])
                if input.bin_annotations:
                    command.extend(["--bin-annotations", input.bin_annotations[0]])
                subprocess.run(command, check=True)

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
                a_candidates=CNV_A_CANDIDATES,
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
                        "--input-a-candidates", input.a_candidates,
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
                        "--hmm-role", CNV_CALLING_HMM_ROLE,
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
                gender_tsv=([CNV_GENDER_TSV] if PREDICT_BY_SEX_ENABLED else []),
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
                    "--a-branch-review-min-abs-z", str(CNV_ARTIFACT_A_BRANCH_REVIEW_MIN_ABS_Z),
                    "--a-branch-sensitive-review-min-abs-z",
                    str(CNV_ARTIFACT_A_BRANCH_SENSITIVE_REVIEW_MIN_ABS_Z),
                    "--a-branch-sensitive-review-min-bins",
                    str(CNV_ARTIFACT_A_BRANCH_SENSITIVE_REVIEW_MIN_BINS),
                    "--a-branch-sensitive-review-max-high-risk-fraction",
                    str(CNV_ARTIFACT_A_BRANCH_SENSITIVE_REVIEW_MAX_HIGH_RISK_FRAC),
                    "--a-branch-sensitive-review-max-region-risk",
                    str(CNV_ARTIFACT_A_BRANCH_SENSITIVE_REVIEW_MAX_REGION_RISK),
                    "--a-branch-sensitive-review-min-same-direction-z",
                    str(CNV_ARTIFACT_A_BRANCH_SENSITIVE_REVIEW_MIN_SAME_DIRECTION_Z),
                    "--a-branch-boundary-protect-min-abs-z",
                    str(CNV_ARTIFACT_A_BRANCH_BOUNDARY_PROTECT_MIN_ABS_Z),
                    "--a-branch-discordant-protect-min-abs-z",
                    str(CNV_ARTIFACT_A_BRANCH_DISCORDANT_PROTECT_MIN_ABS_Z),
                    "--branch-b-direction-min-abs-z", str(CNV_ARTIFACT_BRANCH_B_DIRECTION_MIN_ABS_Z),
                    "--narrow-boundary-artifact-max-bins", str(CNV_ARTIFACT_NARROW_BOUNDARY_MAX_BINS),
                    "--narrow-boundary-artifact-max-available-chrom-fraction",
                    str(CNV_ARTIFACT_NARROW_BOUNDARY_MAX_AVAILABLE_CHROM_FRACTION),
                    "--narrow-boundary-artifact-protect-min-a-abs-z",
                    str(CNV_ARTIFACT_NARROW_BOUNDARY_PROTECT_MIN_A_ABS_Z),
                    "--sca-xy-xgain-max-bam-x-relative", str(CNV_ARTIFACT_SCA_XY_XGAIN_MAX_BAM_X_REL),
                    "--sca-xy-xgain-focal-edge-max-bins", str(CNV_ARTIFACT_SCA_XY_XGAIN_FOCAL_EDGE_MAX_BINS),
                    "--cnvseq-reportable-min-bp", str(CNV_ARTIFACT_CNVSEQ_REPORTABLE_MIN_BP),
                    "--cnvseq-review-min-bp", str(CNV_ARTIFACT_CNVSEQ_REVIEW_MIN_BP),
                    "--cnvseq-large-event-min-bp", str(CNV_ARTIFACT_CNVSEQ_LARGE_EVENT_MIN_BP),
                    "--cnvseq-boundary-max-abs-z", str(CNV_ARTIFACT_CNVSEQ_BOUNDARY_MAX_ABS_Z),
                    "--cnvseq-whole-chrom-available-fraction", str(CNV_ARTIFACT_CNVSEQ_WHOLE_CHROM_AVAILABLE_FRACTION),
                    "--log", log[0],
                ]
                if input.gender_tsv:
                    command.extend(["--gender-tsv", input.gender_tsv[0]])
                for region in CNV_POSTPROCESS_PAR_REGIONS:
                    command.extend(["--par-region", region])
                subprocess.run(command, check=True)

        rule cnv_branch_b_evidence_ledger:
            input:
                a_candidates=CNV_A_CANDIDATES,
                bins=CNV_B_CALIBRATED_BINS,
                events=CNV_B_FINAL_EVENTS,
                metadata=RUN_METADATA
            output:
                ledger=CNV_B_EVIDENCE_LEDGER,
                summary=CNV_B_EVIDENCE_SUMMARY
            log:
                project_path("logs", "cnv", "{sample}.branch_b.evidence_ledger.log")
            threads: 1
            run:
                from pgta.core.logging import write_rule_audit_log
                import subprocess

                write_rule_audit_log(log[0], input.metadata)
                subprocess.run(
                    [
                        config["biosoft"]["python"],
                        SCRIPT_CNV_EVIDENCE_LEDGER,
                        SCRIPT_CNV_EVIDENCE_LEDGER_ACTION,
                        "--sample-id", wildcards.sample,
                        "--a-candidates", input.a_candidates,
                        "--branch-b-events", input.events,
                        "--input-bins", input.bins,
                        "--output-ledger", output.ledger,
                        "--output-summary", output.summary,
                        "--reference-id", CNV_REFERENCE_ID,
                        "--negative-bank-version", CNV_NEGATIVE_BANK_VERSION,
                        "--log", log[0],
                    ],
                    check=True,
                )

        rule cnv_branch_ab_plot:
            input:
                bins=CNV_B_CALIBRATED_BINS,
                events=CNV_B_FINAL_EVENTS,
                a_branch=([CNV_A_ABERRATIONS_BED] if CNV_POSTPROCESS_PRESERVE_BRANCH_A else []),
                metadata=RUN_METADATA
            output:
                svg=CNV_B_PLOT_SVG
            log:
                project_path("logs", "cnv", "{sample}.branch_ab_plot.log")
            threads: 1
            run:
                from pgta.core.logging import write_rule_audit_log
                import subprocess

                write_rule_audit_log(log[0], input.metadata)
                command = [
                    config["biosoft"]["python"],
                    SCRIPT_CNV_PLOT,
                    SCRIPT_CNV_PLOT_ACTION,
                    "--sample-id", wildcards.sample,
                    "--input-bins", input.bins,
                    "--input-events", input.events,
                    "--output-svg", output.svg,
                    "--log", log[0],
                ]
                if input.a_branch:
                    command.extend(["--input-a-branch", input.a_branch[0]])
                subprocess.run(command, check=True)

    if CNV_NEGATIVE_BANK_SAMPLES_TSV:
        rule cnv_negative_bank_labels:
            input:
                samples=CNV_NEGATIVE_BANK_SAMPLES_TSV,
                metadata=RUN_METADATA
            output:
                labels=CNV_NEGATIVE_BANK_LABELS,
                summary=CNV_NEGATIVE_BANK_SUMMARY
            log:
                project_path("logs", "cnv", "negative_bank_labels.log")
            threads: 1
            run:
                from pgta.core.logging import write_rule_audit_log
                import subprocess

                write_rule_audit_log(log[0], input.metadata)
                subprocess.run(
                    [
                        config["biosoft"]["python"],
                        SCRIPT_NEGATIVE_BANK_LABELS,
                        SCRIPT_NEGATIVE_BANK_LABELS_ACTION,
                        "--input-samples", input.samples,
                        "--output-labels", output.labels,
                        "--output-summary", output.summary,
                        "--version", CNV_NEGATIVE_BANK_VERSION,
                        "--log", log[0],
                    ],
                    check=True,
                )

        if CNV_POSTPROCESS_ENABLE_BRANCH_B:
            rule cnv_branch_b_matched_negative:
                input:
                    ledger=CNV_B_EVIDENCE_LEDGER,
                    labels=CNV_NEGATIVE_BANK_LABELS,
                    background_ledgers=CNV_NEGATIVE_BANK_BACKGROUND_LEDGERS + expand(CNV_B_EVIDENCE_LEDGER, sample=SAMPLES),
                    metadata=RUN_METADATA
                output:
                    ledger=CNV_B_MATCHED_NEGATIVE,
                    summary=CNV_B_MATCHED_NEGATIVE_SUMMARY
                log:
                    project_path("logs", "cnv", "{sample}.branch_b.matched_negative.log")
                threads: 1
                run:
                    from pgta.core.logging import write_rule_audit_log
                    import subprocess

                    write_rule_audit_log(log[0], input.metadata)
                    command = [
                        config["biosoft"]["python"],
                        SCRIPT_MATCHED_NEGATIVE,
                        SCRIPT_MATCHED_NEGATIVE_ACTION,
                        "--sample-id", wildcards.sample,
                        "--input-ledger", input.ledger,
                        "--negative-bank-labels", input.labels,
                        "--output-ledger", output.ledger,
                        "--output-summary", output.summary,
                        "--version", CNV_NEGATIVE_BANK_VERSION,
                        "--feature-column", CNV_NEGATIVE_BANK_FEATURE_COLUMN,
                        "--min-background", str(CNV_NEGATIVE_BANK_MIN_BACKGROUND),
                        "--similar-length-fold", str(CNV_NEGATIVE_BANK_SIMILAR_LENGTH_FOLD),
                        "--log", log[0],
                    ]
                    for label in CNV_NEGATIVE_BANK_SHADOW_BACKGROUND_LABELS:
                        command.extend(["--shadow-background-label", label])
                    for background_ledger in input.background_ledgers:
                        command.extend(["--background-ledger", background_ledger])
                    subprocess.run(command, check=True)

    if CNV_POSTPROCESS_ENABLE_BRANCH_B:
        rule cnv_branch_b_v2_classifier:
            input:
                ledger=(CNV_B_MATCHED_NEGATIVE if CNV_NEGATIVE_BANK_SAMPLES_TSV else CNV_B_EVIDENCE_LEDGER),
                metadata=RUN_METADATA
            output:
                classification=CNV_B_V2_CLASSIFIER,
                summary=CNV_B_V2_CLASSIFIER_SUMMARY
            log:
                project_path("logs", "cnv", "{sample}.branch_b.v2_classifier.log")
            threads: 1
            run:
                from pgta.core.logging import write_rule_audit_log
                import subprocess

                write_rule_audit_log(log[0], input.metadata)
                subprocess.run(
                    [
                        config["biosoft"]["python"],
                        SCRIPT_BRANCH_B_V2_CLASSIFIER,
                        SCRIPT_BRANCH_B_V2_CLASSIFIER_ACTION,
                        "--sample-id", wildcards.sample,
                        "--input-ledger", input.ledger,
                        "--output-classification", output.classification,
                        "--output-summary", output.summary,
                        "--version", CNV_NEGATIVE_BANK_VERSION,
                        "--log", log[0],
                    ],
                    check=True,
                )

        if "branch_b_v2_benchmark" in AVAILABLE_TARGETS:
            rule cnv_branch_b_v2_benchmark:
                input:
                    classifications=expand(CNV_B_V2_CLASSIFIER, sample=SAMPLES),
                    truth=([CNV_BENCHMARK_TRUTH_TSV] if CNV_BENCHMARK_TRUTH_TSV else []),
                    metadata=RUN_METADATA
                output:
                    truth_metrics=CNV_B_V2_BENCHMARK_TRUTH_METRICS,
                    sample_summary=CNV_B_V2_BENCHMARK_SAMPLE_SUMMARY,
                    summary=CNV_B_V2_BENCHMARK_SUMMARY
                log:
                    project_path("logs", "cnv", "branch_b_v2_benchmark.log")
                threads: 1
                run:
                    from pgta.core.logging import write_rule_audit_log
                    import subprocess

                    write_rule_audit_log(log[0], input.metadata)
                    command = [
                        config["biosoft"]["python"],
                        SCRIPT_BRANCH_B_V2_BENCHMARK,
                        SCRIPT_BRANCH_B_V2_BENCHMARK_ACTION,
                        "--reference-id", CNV_REFERENCE_ID,
                        "--output-truth-metrics", output.truth_metrics,
                        "--output-sample-summary", output.sample_summary,
                        "--output-summary", output.summary,
                    ]
                    if input.truth:
                        command.extend(["--truth-tsv", str(input.truth[0])])
                    for sample_id in SAMPLES:
                        command.extend(["--sample-id", str(sample_id)])
                    for path_value in input.classifications:
                        command.extend(["--classification-tsv", path_value])
                    subprocess.run(command, check=True)

        rule cnv_branch_s_shadow:
            input:
                bins=CNV_B_CALIBRATED_BINS,
                a_candidates=CNV_A_CANDIDATES,
                gender_tsv=([CNV_GENDER_TSV] if PREDICT_BY_SEX_ENABLED else []),
                metadata=RUN_METADATA
            output:
                evidence=CNV_BRANCH_S_EVIDENCE,
                scores=CNV_BRANCH_S_SCORES,
                summary=CNV_BRANCH_S_SUMMARY
            log:
                project_path("logs", "cnv", "{sample}.branch_s.shadow.log")
            threads: 1
            run:
                from pgta.core.logging import write_rule_audit_log
                import subprocess

                write_rule_audit_log(log[0], input.metadata)
                command = [
                    config["biosoft"]["python"],
                    SCRIPT_BRANCH_S_SHADOW,
                    SCRIPT_BRANCH_S_SHADOW_ACTION,
                    "--sample-id", wildcards.sample,
                    "--input-bins", input.bins,
                    "--a-candidates", input.a_candidates,
                    "--output-evidence", output.evidence,
                    "--output-scores", output.scores,
                    "--output-summary", output.summary,
                    "--log", log[0],
                ]
                if input.gender_tsv:
                    command.extend(["--gender-tsv", input.gender_tsv[0]])
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
                plots=(expand(CNV_B_PLOT_SVG, sample=SAMPLES) if CNV_POSTPROCESS_ENABLE_BRANCH_B else []),
                genders=(expand(CNV_GENDER_TSV, sample=SAMPLES) if PREDICT_BY_SEX_ENABLED else []),
                qcs=expand(CNV_QC_TSV, sample=SAMPLES),
                a_branch=(expand(CNV_A_ABERRATIONS_BED, sample=SAMPLES) if CNV_POSTPROCESS_PRESERVE_BRANCH_A else []),
                branch_a_validation_summaries=CNV_REPORT_BRANCH_A_VALIDATION_SUMMARIES,
                branch_b_evidence_summaries=(expand(CNV_B_EVIDENCE_SUMMARY, sample=SAMPLES) if CNV_POSTPROCESS_ENABLE_BRANCH_B else []),
                branch_s_summaries=(expand(CNV_BRANCH_S_SUMMARY, sample=SAMPLES) if CNV_POSTPROCESS_ENABLE_BRANCH_B else []),
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
            params:
                reference_id=CNV_REFERENCE_ID,
                wisecondorx_predict_command=(
                    "WisecondorX predict <npz> <ref> <output-prefix> "
                    f"--zscore {CNV_ZSCORE} --alpha {CNV_ALPHA} "
                    f"--maskrepeats {CNV_MASKREPEATS} --minrefbins {CNV_MINREFBINS}"
                )
            threads: 1
            run:
                from pgta.core.logging import write_rule_audit_log
                import subprocess

                write_rule_audit_log(log[0], input.metadata)
                command = [
                    config["biosoft"]["python"],
                    SCRIPT_CNV_REPORT,
                    SCRIPT_CNV_REPORT_ACTION,
                    "--reference-id", params.reference_id,
                    "--wisecondorx-predict-command", params.wisecondorx_predict_command,
                    "--output-tsv", output.tsv,
                    "--output-json", output.json,
                    "--output-md", output.md,
                    "--output-html", output.html,
                    "--log", log[0],
                ]
                if input.evaluation_summary:
                    command.extend(["--evaluation-summary", input.evaluation_summary[0]])
                if input.ml_summary:
                    command.extend(["--ml-summary", input.ml_summary[0]])
                if input.benchmark_summary:
                    command.extend(["--benchmark-summary", input.benchmark_summary[0]])
                if input.mosaic_truth_validation:
                    command.extend(["--truth-validation-summary", input.mosaic_truth_validation[0]])
                for path_value in input.branch_a_validation_summaries:
                    command.extend(["--branch-a-validation-summary", path_value])
                for path_value in input.events:
                    command.extend(["--event-tsv", path_value])
                for path_value in input.plots:
                    command.extend(["--plot-svg", path_value])
                for path_value in input.genders:
                    command.extend(["--gender-tsv", path_value])
                for path_value in input.qcs:
                    command.extend(["--qc-tsv", path_value])
                for path_value in input.a_branch:
                    command.extend(["--a-branch-bed", path_value])
                for path_value in input.branch_b_evidence_summaries:
                    command.extend(["--branch-b-evidence-summary", path_value])
                for path_value in input.branch_s_summaries:
                    command.extend(["--branch-s-summary", path_value])
                subprocess.run(command, check=True)
