from pathlib import Path
import shutil


REFERENCE_BAMS = [SORTED_BAM.format(sample=sample_id) for sample_id in REF_SAMPLE_IDS]
REFERENCE_SAMPLE_TEXT = "\n".join(REF_SAMPLE_IDS)


def load_sample_ids_from_file(path_value):
    return [
        line.strip()
        for line in Path(path_value).read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]


def select_allowed_samples(path_value, allowed_samples):
    selected = load_sample_ids_from_file(path_value)
    allowed = {str(sample_id).strip() for sample_id in allowed_samples if str(sample_id).strip()}
    if allowed:
        selected = [sample_id for sample_id in selected if sample_id in allowed]
    if not selected:
        raise ValueError(f"No samples available after applying allowed set to {path_value}")
    return selected


def resolve_bams_for_sample_ids(sample_ids):
    return [SORTED_BAM.format(sample=sample_id) for sample_id in sample_ids]


def publish_reference_output(source_path, output_path, logger):
    src = Path(source_path)
    dst = Path(output_path)
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, dst)
    logger.info("published default reference output: %s <- %s", dst, src)


if TUNING_ENABLED:
    rule select_reference_cohort:
        input:
            summary=BASELINE_QC_SUMMARY,
            metadata=RUN_METADATA
        output:
            sample_list=project_path("reference", "cohorts", "{cohort}", "selected_samples.txt")
        log:
            project_path("logs", "wisecondorx", "reference_select_{cohort}.log")
        benchmark:
            str(Path(BENCHMARK_ROOT) / "reference" / "select_cohort" / "{cohort}.tsv")
        params:
            decisions=lambda wildcards: ",".join(REF_SET_CFG[wildcards.cohort]["decisions"])
        threads: 1
        run:
            from pgta.core.logging import setup_logger, write_rule_audit_log
            from pgta.reference.cohort import load_selected_samples, parse_decisions

            write_rule_audit_log(log[0], input.metadata, [("REFERENCE COHORT", wildcards.cohort), ("QC DECISIONS", params.decisions)])
            logger = setup_logger("select_reference_cohort", log[0])
            selected = load_selected_samples(input.summary, parse_decisions(params.decisions))
            output_path = Path(output.sample_list)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            output_path.write_text("".join(f"{sample_id}\n" for sample_id in selected), encoding="utf-8")
            logger.info(
                "selected reference cohort=%s decisions=%s samples=%d",
                wildcards.cohort,
                params.decisions,
                len(selected),
            )

    if BUILD_REF_BY_SEX_ENABLED:
        rule reference_prefilter_cohort_sex:
            input:
                cohort_samples=lambda wildcards: REF_SET_SAMPLE_LISTS[wildcards.cohort],
                metadata=RUN_METADATA
            output:
                qc=project_path("reference", "cohorts", "{cohort}", "{sex}", "prefilter", "reference_sample_qc.tsv"),
                plot=project_path("reference", "cohorts", "{cohort}", "{sex}", "prefilter", "reference_sample_qc.svg"),
                inliers=project_path("reference", "cohorts", "{cohort}", "{sex}", "prefilter", "reference_inlier_samples.txt"),
                summary=project_path("reference", "cohorts", "{cohort}", "{sex}", "prefilter", "prefilter_summary.yaml")
            log:
                project_path("logs", "wisecondorx", "reference_prefilter_{cohort}_{sex}.log")
            benchmark:
                str(Path(BENCHMARK_ROOT) / "reference" / "reference_prefilter" / "{cohort}" / "{sex}.tsv")
            params:
                wise=config["biosoft"]["WisecondorX"],
                binsize=PREFILTER_BINSIZE,
                pca_min=TUNING_PCA_MIN,
                pca_max=TUNING_PCA_MAX,
                min_ref_samples=lambda wildcards: int(
                    REF_SET_CFG[wildcards.cohort].get("min_reference_samples_override") or TUNING_MIN_REF_SAMPLES
                ),
                min_reads=TUNING_MIN_READS,
                min_corr=TUNING_MIN_CORR,
                max_recon_z=TUNING_MAX_RECON_Z,
                max_noise_z=TUNING_MAX_NOISE_Z,
                max_iter=PREFILTER_MAX_ITER,
                allowed_samples=lambda wildcards: REF_SAMPLE_IDS_BY_SEX[wildcards.sex]
            threads: 4
            run:
                from pgta.core.logging import setup_logger, write_rule_audit_log
                from pgta.reference.prefilter import run_reference_prefilter_qc

                sample_ids = select_allowed_samples(input.cohort_samples, params.allowed_samples)
                sample_text = "\n".join(sample_ids)
                write_rule_audit_log(
                    log[0],
                    input.metadata,
                    [("REFERENCE COHORT", wildcards.cohort), ("SEX GROUP", wildcards.sex), ("REFERENCE SAMPLES", sample_text)],
                )
                logger = setup_logger("reference_prefilter_qc", log[0])
                run_reference_prefilter_qc(
                    wisecondorx=params.wise,
                    bams=resolve_bams_for_sample_ids(sample_ids),
                    sample_ids=sample_ids,
                    binsize=params.binsize,
                    pca_min_components=params.pca_min,
                    pca_max_components=params.pca_max,
                    min_reference_samples=params.min_ref_samples,
                    min_reads_per_sample=params.min_reads,
                    min_corr_to_median=params.min_corr,
                    max_reconstruction_error_z=params.max_recon_z,
                    max_noise_mad_z=params.max_noise_z,
                    max_iterations=params.max_iter,
                    threads=threads,
                    workdir=str(Path(output.qc).parent),
                    qc_output=output.qc,
                    plot_output=output.plot,
                    inlier_samples_output=output.inliers,
                    summary_output=output.summary,
                    logger=logger,
                )

        rule merge_reference_prefilter_inliers_by_cohort:
            input:
                xx_inliers=project_path("reference", "cohorts", "{cohort}", "XX", "prefilter", "reference_inlier_samples.txt"),
                xy_inliers=project_path("reference", "cohorts", "{cohort}", "XY", "prefilter", "reference_inlier_samples.txt"),
                metadata=RUN_METADATA
            output:
                inliers=project_path("reference", "cohorts", "{cohort}", "prefilter", "reference_inlier_samples.txt")
            log:
                project_path("logs", "wisecondorx", "reference_prefilter_merge_{cohort}.log")
            threads: 1
            run:
                from pgta.core.logging import setup_logger, write_rule_audit_log

                write_rule_audit_log(log[0], input.metadata, [("REFERENCE COHORT", wildcards.cohort), ("REFERENCE SAMPLES", REFERENCE_SAMPLE_TEXT)])
                logger = setup_logger("merge_reference_prefilter_inliers", log[0])
                merged = []
                for source in [input.xx_inliers, input.xy_inliers]:
                    merged.extend(load_sample_ids_from_file(source))
                merged = sorted(dict.fromkeys(merged))
                output_path = Path(output.inliers)
                output_path.parent.mkdir(parents=True, exist_ok=True)
                output_path.write_text("".join(f"{sample_id}\n" for sample_id in merged), encoding="utf-8")
                logger.info("merged prefilter inliers for cohort=%s: %d samples", wildcards.cohort, len(merged))

        rule tune_wisecondorx_reference_qc_by_cohort:
            input:
                cohort_samples=lambda wildcards: REF_SET_SAMPLE_LISTS[wildcards.cohort],
                prefilter_inliers=project_path("reference", "cohorts", "{cohort}", "prefilter", "reference_inlier_samples.txt"),
                metadata=RUN_METADATA
            output:
                summary=project_path("wisecondorx", "tuning", "{cohort}", "bin_pca_grid.tsv"),
                binsize_summary=project_path("wisecondorx", "tuning", "{cohort}", "binsize_ranking.tsv"),
                best=project_path("wisecondorx", "tuning", "{cohort}", "best_params.yaml"),
                qc=project_path("wisecondorx", "tuning", "{cohort}", "reference_sample_qc.tsv"),
                plot=project_path("wisecondorx", "tuning", "{cohort}", "best_bin_pca_elbow.svg"),
                qc_plot=project_path("wisecondorx", "tuning", "{cohort}", "reference_qc_metrics.svg"),
                inliers=project_path("wisecondorx", "tuning", "{cohort}", "reference_inlier_samples.txt")
            log:
                project_path("logs", "wisecondorx", "tuning_{cohort}.log")
            benchmark:
                str(Path(BENCHMARK_ROOT) / "reference" / "tune_wisecondorx_reference_qc" / "{cohort}.tsv")
            params:
                wise=config["biosoft"]["WisecondorX"],
                bin_sizes=",".join(str(item) for item in TUNING_BIN_SIZES),
                pca_min=TUNING_PCA_MIN,
                pca_max=TUNING_PCA_MAX,
                pca_min_var=TUNING_MIN_VAR,
                min_ref_samples=lambda wildcards: int(
                    REF_SET_CFG[wildcards.cohort].get("min_reference_samples_override") or TUNING_MIN_REF_SAMPLES
                ),
                max_outlier_fraction=TUNING_MAX_OUTLIER_FRAC,
                min_reads=TUNING_MIN_READS,
                min_corr=TUNING_MIN_CORR,
                max_recon_z=TUNING_MAX_RECON_Z,
                max_noise_z=TUNING_MAX_NOISE_Z
            threads: 4
            run:
                from pgta.core.logging import setup_logger, write_rule_audit_log
                from pgta.reference.tune import run_tune_wisecondorx

                sample_ids = load_sample_ids_from_file(input.cohort_samples)
                sample_text = "\n".join(sample_ids)
                write_rule_audit_log(
                    log[0],
                    input.metadata,
                    [("REFERENCE COHORT", wildcards.cohort), ("REFERENCE SAMPLES", sample_text)],
                )
                logger = setup_logger("tune_wisecondorx_bin_pca", log[0])
                run_tune_wisecondorx(
                    wisecondorx=params.wise,
                    bams=resolve_bams_for_sample_ids(sample_ids),
                    sample_ids=sample_ids,
                    allowed_samples_file=input.prefilter_inliers,
                    bin_sizes=params.bin_sizes,
                    pca_min_components=params.pca_min,
                    pca_max_components=params.pca_max,
                    pca_min_explained_variance=params.pca_min_var,
                    min_reference_samples=params.min_ref_samples,
                    max_outlier_fraction=params.max_outlier_fraction,
                    min_reads_per_sample=params.min_reads,
                    min_corr_to_median=params.min_corr,
                    max_reconstruction_error_z=params.max_recon_z,
                    max_noise_mad_z=params.max_noise_z,
                    threads=threads,
                    workdir=REF_SET_TUNING_WORKDIR[wildcards.cohort],
                    summary_output=output.summary,
                    binsize_summary_output=output.binsize_summary,
                    best_output=output.best,
                    qc_output=output.qc,
                    plot_output=output.plot,
                    qc_stats_plot_output=output.qc_plot,
                    inlier_samples_output=output.inliers,
                    reference_output=REF_SET_GENDER_OUTPUT[wildcards.cohort],
                    skip_build_reference=True,
                    logger=logger,
                )

        rule write_common_reference_binsize_from_tuning_by_cohort:
            input:
                best=project_path("wisecondorx", "tuning", "{cohort}", "best_params.yaml"),
                metadata=RUN_METADATA
            output:
                binsize=project_path("reference", "cohorts", "{cohort}", "gender", "common_best_binsize.txt")
            log:
                project_path("logs", "wisecondorx", "reference_common_binsize_{cohort}.log")
            threads: 1
            run:
                from pgta.core.logging import setup_logger, write_rule_audit_log
                from pgta.reference.build import load_best_binsize

                write_rule_audit_log(log[0], input.metadata, [("REFERENCE COHORT", wildcards.cohort)])
                logger = setup_logger("write_common_reference_binsize_from_tuning", log[0])
                best_binsize = load_best_binsize(input.best)
                output_path = Path(output.binsize)
                output_path.parent.mkdir(parents=True, exist_ok=True)
                output_path.write_text(f"{best_binsize}\n", encoding="utf-8")
                logger.info("common reference binsize written for cohort=%s: %s", wildcards.cohort, best_binsize)

        rule build_wisecondorx_reference_from_tuning_by_cohort_sex:
            wildcard_constraints:
                sex="XX|XY"
            input:
                best=project_path("wisecondorx", "tuning", "{cohort}", "best_params.yaml"),
                inliers=project_path("wisecondorx", "tuning", "{cohort}", "reference_inlier_samples.txt"),
                metadata=RUN_METADATA
            output:
                ref=project_path("reference", "cohorts", "{cohort}", "{sex}", "result", "ref_{sex}_best.npz")
            log:
                project_path("logs", "wisecondorx", "build_reference_{cohort}_{sex}.log")
            benchmark:
                str(Path(BENCHMARK_ROOT) / "reference" / "build_reference" / "{cohort}" / "{sex}.tsv")
            params:
                wise=config["biosoft"]["WisecondorX"],
                allowed_samples=lambda wildcards: ",".join(REF_SAMPLE_IDS_BY_SEX[wildcards.sex])
            threads: 4
            run:
                from pgta.core.logging import setup_logger, write_rule_audit_log
                from pgta.reference.build import build_reference_from_tuning

                sample_text = "\n".join(REF_SAMPLE_IDS_BY_SEX[wildcards.sex])
                write_rule_audit_log(
                    log[0],
                    input.metadata,
                    [("REFERENCE COHORT", wildcards.cohort), ("SEX GROUP", wildcards.sex), ("REFERENCE SAMPLES", sample_text)],
                )
                logger = setup_logger("build_reference_from_tuning", log[0])
                build_reference_from_tuning(
                    wisecondorx=params.wise,
                    best_yaml=input.best,
                    inlier_samples=input.inliers,
                    allowed_samples=params.allowed_samples,
                    tuning_workdir=REF_SET_TUNING_WORKDIR[wildcards.cohort],
                    reference_output=output.ref,
                    threads=threads,
                    logger=logger,
                )

        rule build_wisecondorx_gender_reference_from_tuning_by_cohort:
            input:
                best=project_path("wisecondorx", "tuning", "{cohort}", "best_params.yaml"),
                inliers=project_path("wisecondorx", "tuning", "{cohort}", "reference_inlier_samples.txt"),
                common_binsize=project_path("reference", "cohorts", "{cohort}", "gender", "common_best_binsize.txt"),
                metadata=RUN_METADATA
            output:
                ref=project_path("reference", "cohorts", "{cohort}", "gender", "result", "ref_gender_best.npz")
            log:
                project_path("logs", "wisecondorx", "build_reference_gender_{cohort}.log")
            benchmark:
                str(Path(BENCHMARK_ROOT) / "reference" / "build_gender_reference" / "{cohort}.tsv")
            params:
                wise=config["biosoft"]["WisecondorX"]
            threads: 4
            run:
                from pgta.core.logging import setup_logger, write_rule_audit_log
                from pgta.reference.build import build_reference_from_npz_paths, resolve_inlier_npz_paths

                write_rule_audit_log(log[0], input.metadata, [("REFERENCE COHORT", wildcards.cohort), ("REFERENCE SAMPLES", REFERENCE_SAMPLE_TEXT)])
                logger = setup_logger("build_gender_reference_from_tuning", log[0])
                expected_binsize = int(Path(input.common_binsize).read_text(encoding="utf-8").strip())
                best_binsize, inlier_ids, npz_paths = resolve_inlier_npz_paths(
                    best_yaml=input.best,
                    inlier_samples=input.inliers,
                    tuning_workdir=REF_SET_TUNING_WORKDIR[wildcards.cohort],
                )
                if best_binsize != expected_binsize:
                    raise ValueError(
                        f"Common binsize mismatch for cohort={wildcards.cohort}: expected={expected_binsize}, best={best_binsize}"
                    )
                logger.info("building gender reference for cohort=%s from inlier samples=%d", wildcards.cohort, len(inlier_ids))
                build_reference_from_npz_paths(
                    wisecondorx=params.wise,
                    npz_paths=npz_paths,
                    reference_output=output.ref,
                    binsize=expected_binsize,
                    threads=threads,
                    logger=logger,
                )

        rule publish_default_reference_output_xx:
            input:
                src=REF_SET_OUTPUTS_BY_SEX[DEFAULT_REF_SET]["XX"],
                metadata=RUN_METADATA
            output:
                ref=REF_OUTPUTS_BY_SEX["XX"]
            log:
                project_path("logs", "wisecondorx", "publish_default_reference_XX.log")
            threads: 1
            run:
                from pgta.core.logging import setup_logger, write_rule_audit_log

                write_rule_audit_log(log[0], input.metadata, [("DEFAULT REFERENCE COHORT", DEFAULT_REF_SET), ("SEX GROUP", "XX")])
                logger = setup_logger("publish_default_reference", log[0])
                publish_reference_output(input.src, output.ref, logger)

        rule publish_default_reference_output_xy:
            input:
                src=REF_SET_OUTPUTS_BY_SEX[DEFAULT_REF_SET]["XY"],
                metadata=RUN_METADATA
            output:
                ref=REF_OUTPUTS_BY_SEX["XY"]
            log:
                project_path("logs", "wisecondorx", "publish_default_reference_XY.log")
            threads: 1
            run:
                from pgta.core.logging import setup_logger, write_rule_audit_log

                write_rule_audit_log(log[0], input.metadata, [("DEFAULT REFERENCE COHORT", DEFAULT_REF_SET), ("SEX GROUP", "XY")])
                logger = setup_logger("publish_default_reference", log[0])
                publish_reference_output(input.src, output.ref, logger)

        rule publish_default_gender_reference_output:
            input:
                src=REF_SET_GENDER_OUTPUT[DEFAULT_REF_SET],
                metadata=RUN_METADATA
            output:
                ref=GENDER_REF_OUTPUT
            log:
                project_path("logs", "wisecondorx", "publish_default_reference_gender.log")
            threads: 1
            run:
                from pgta.core.logging import setup_logger, write_rule_audit_log

                write_rule_audit_log(log[0], input.metadata, [("DEFAULT REFERENCE COHORT", DEFAULT_REF_SET)])
                logger = setup_logger("publish_default_reference", log[0])
                publish_reference_output(input.src, output.ref, logger)

        rule publish_default_common_reference_binsize:
            input:
                src=REF_SET_COMMON_BINSIZE[DEFAULT_REF_SET],
                metadata=RUN_METADATA
            output:
                binsize=COMMON_REF_BINSIZE
            log:
                project_path("logs", "wisecondorx", "publish_default_common_binsize.log")
            threads: 1
            run:
                from pgta.core.logging import setup_logger, write_rule_audit_log

                write_rule_audit_log(log[0], input.metadata, [("DEFAULT REFERENCE COHORT", DEFAULT_REF_SET)])
                logger = setup_logger("publish_default_reference", log[0])
                publish_reference_output(input.src, output.binsize, logger)
    else:
        rule reference_prefilter:
            input:
                bams=REFERENCE_BAMS,
                metadata=RUN_METADATA
            output:
                qc=str(Path(REF_PREFILTER_DIR) / "reference_sample_qc.tsv"),
                plot=str(Path(REF_PREFILTER_DIR) / "reference_sample_qc.svg"),
                inliers=str(Path(REF_PREFILTER_DIR) / "reference_inlier_samples.txt"),
                summary=str(Path(REF_PREFILTER_DIR) / "prefilter_summary.yaml")
            log:
                project_path("logs", "wisecondorx", "reference_prefilter.log")
            benchmark:
                BENCH_REFERENCE_PREFILTER
            params:
                wise=config["biosoft"]["WisecondorX"],
                binsize=PREFILTER_BINSIZE,
                pca_min=TUNING_PCA_MIN,
                pca_max=TUNING_PCA_MAX,
                min_ref_samples=TUNING_MIN_REF_SAMPLES,
                min_reads=TUNING_MIN_READS,
                min_corr=TUNING_MIN_CORR,
                max_recon_z=TUNING_MAX_RECON_Z,
                max_noise_z=TUNING_MAX_NOISE_Z,
                max_iter=PREFILTER_MAX_ITER,
                workdir=REF_PREFILTER_DIR,
                sample_ids=REF_SAMPLE_IDS,
                sample_text=REFERENCE_SAMPLE_TEXT
            threads: 4
            run:
                from pgta.core.logging import setup_logger, write_rule_audit_log
                from pgta.reference.prefilter import run_reference_prefilter_qc

                write_rule_audit_log(log[0], input.metadata, [("REFERENCE SAMPLES", params.sample_text)])
                logger = setup_logger("reference_prefilter_qc", log[0])
                run_reference_prefilter_qc(
                    wisecondorx=params.wise,
                    bams=input.bams,
                    sample_ids=params.sample_ids,
                    binsize=params.binsize,
                    pca_min_components=params.pca_min,
                    pca_max_components=params.pca_max,
                    min_reference_samples=params.min_ref_samples,
                    min_reads_per_sample=params.min_reads,
                    min_corr_to_median=params.min_corr,
                    max_reconstruction_error_z=params.max_recon_z,
                    max_noise_mad_z=params.max_noise_z,
                    max_iterations=params.max_iter,
                    threads=threads,
                    workdir=params.workdir,
                    qc_output=output.qc,
                    plot_output=output.plot,
                    inlier_samples_output=output.inliers,
                    summary_output=output.summary,
                    logger=logger,
                )

        rule tune_wisecondorx_reference_qc:
            input:
                bams=REFERENCE_BAMS,
                prefilter_inliers=str(Path(REF_PREFILTER_DIR) / "reference_inlier_samples.txt"),
                metadata=RUN_METADATA
            output:
                summary=TUNING_SUMMARY,
                binsize_summary=TUNING_BINSIZE_SUMMARY,
                best=TUNING_BEST,
                qc=TUNING_QC,
                plot=TUNING_PLOT,
                qc_plot=TUNING_QC_STATS_PLOT,
                inliers=TUNING_INLIERS
            log:
                project_path("logs", "wisecondorx", "tuning.log")
            benchmark:
                BENCH_TUNE_WISECONDORX_REFERENCE_QC
            params:
                wise=config["biosoft"]["WisecondorX"],
                bin_sizes=",".join(str(item) for item in TUNING_BIN_SIZES),
                sample_ids=REF_SAMPLE_IDS,
                pca_min=TUNING_PCA_MIN,
                pca_max=TUNING_PCA_MAX,
                pca_min_var=TUNING_MIN_VAR,
                min_ref_samples=TUNING_MIN_REF_SAMPLES,
                max_outlier_fraction=TUNING_MAX_OUTLIER_FRAC,
                min_reads=TUNING_MIN_READS,
                min_corr=TUNING_MIN_CORR,
                max_recon_z=TUNING_MAX_RECON_Z,
                max_noise_z=TUNING_MAX_NOISE_Z,
                workdir=TUNING_WORKDIR,
                reference_output=REF_OUTPUT,
                sample_text=REFERENCE_SAMPLE_TEXT
            threads: 4
            run:
                from pgta.core.logging import setup_logger, write_rule_audit_log
                from pgta.reference.tune import run_tune_wisecondorx

                write_rule_audit_log(log[0], input.metadata, [("REFERENCE SAMPLES", params.sample_text)])
                logger = setup_logger("tune_wisecondorx_bin_pca", log[0])
                run_tune_wisecondorx(
                    wisecondorx=params.wise,
                    bams=input.bams,
                    sample_ids=params.sample_ids,
                    allowed_samples_file=input.prefilter_inliers,
                    bin_sizes=params.bin_sizes,
                    pca_min_components=params.pca_min,
                    pca_max_components=params.pca_max,
                    pca_min_explained_variance=params.pca_min_var,
                    min_reference_samples=params.min_ref_samples,
                    max_outlier_fraction=params.max_outlier_fraction,
                    min_reads_per_sample=params.min_reads,
                    min_corr_to_median=params.min_corr,
                    max_reconstruction_error_z=params.max_recon_z,
                    max_noise_mad_z=params.max_noise_z,
                    threads=threads,
                    workdir=params.workdir,
                    summary_output=output.summary,
                    binsize_summary_output=output.binsize_summary,
                    best_output=output.best,
                    qc_output=output.qc,
                    plot_output=output.plot,
                    qc_stats_plot_output=output.qc_plot,
                    inlier_samples_output=output.inliers,
                    reference_output=params.reference_output,
                    skip_build_reference=True,
                    logger=logger,
                )

        rule build_wisecondorx_reference_from_tuning:
            input:
                best=TUNING_BEST,
                inliers=TUNING_INLIERS,
                metadata=RUN_METADATA
            output:
                ref=REF_OUTPUT
            log:
                project_path("logs", "wisecondorx", "build_reference.log")
            benchmark:
                BENCH_BUILD_REFERENCE
            params:
                wise=config["biosoft"]["WisecondorX"],
                workdir=TUNING_WORKDIR,
                allowed_samples=",".join(REF_SAMPLE_IDS),
                sample_text=REFERENCE_SAMPLE_TEXT
            threads: 4
            run:
                from pgta.core.logging import setup_logger, write_rule_audit_log
                from pgta.reference.build import build_reference_from_tuning

                write_rule_audit_log(log[0], input.metadata, [("REFERENCE SAMPLES", params.sample_text)])
                logger = setup_logger("build_reference_from_tuning", log[0])
                build_reference_from_tuning(
                    wisecondorx=params.wise,
                    best_yaml=input.best,
                    inlier_samples=input.inliers,
                    allowed_samples=params.allowed_samples,
                    tuning_workdir=params.workdir,
                    reference_output=output.ref,
                    threads=threads,
                    logger=logger,
                )
else:
    if BUILD_REF_BY_SEX_ENABLED:
        rule write_reference_common_binsize_fixed:
            input:
                metadata=RUN_METADATA
            output:
                binsize=COMMON_REF_BINSIZE
            log:
                project_path("logs", "wisecondorx", "reference_common_binsize.log")
            params:
                binsize=WISE_CFG["binsize"]
            threads: 1
            run:
                from pgta.core.logging import setup_logger, write_rule_audit_log

                write_rule_audit_log(log[0], input.metadata)
                logger = setup_logger("write_reference_common_binsize_fixed", log[0])
                output_path = Path(output.binsize)
                output_path.parent.mkdir(parents=True, exist_ok=True)
                output_path.write_text(f"{int(params.binsize)}\n", encoding="utf-8")
                logger.info("fixed common reference binsize written: %s", params.binsize)

        rule build_wisecondorx_reference_fixed_xx:
            input:
                bams=[SORTED_BAM.format(sample=sample_id) for sample_id in REF_SAMPLE_IDS_BY_SEX["XX"]],
                metadata=RUN_METADATA
            output:
                ref=REF_OUTPUTS_BY_SEX["XX"]
            log:
                project_path("logs", "wisecondorx", "build_reference_XX.log")
            benchmark:
                BENCH_BUILD_REFERENCE_BY_SEX["XX"]
            params:
                wise=config["biosoft"]["WisecondorX"],
                binsize=WISE_CFG["binsize"],
                converted_dir=project_path("wisecondorx", "converted", "XX"),
                sample_ids=REF_SAMPLE_IDS_BY_SEX["XX"],
                sample_text="\n".join(REF_SAMPLE_IDS_BY_SEX["XX"])
            threads: 4
            run:
                from pgta.core.logging import setup_logger, write_rule_audit_log
                from pgta.reference.tune import build_reference, convert_all_bams

                write_rule_audit_log(log[0], input.metadata, [("SEX GROUP", "XX"), ("REFERENCE SAMPLES", params.sample_text)])
                logger = setup_logger("build_wisecondorx_reference_fixed", log[0])
                npz_paths = convert_all_bams(
                    wisecondorx=params.wise,
                    bams=input.bams,
                    sample_ids=params.sample_ids,
                    binsize=params.binsize,
                    output_dir=Path(params.converted_dir),
                    threads=threads,
                    logger=logger,
                )
                build_reference(
                    wisecondorx=params.wise,
                    binsize=params.binsize,
                    npz_paths=npz_paths,
                    reference_output=Path(output.ref),
                    threads=threads,
                    logger=logger,
                )

        rule build_wisecondorx_reference_fixed_xy:
            input:
                bams=[SORTED_BAM.format(sample=sample_id) for sample_id in REF_SAMPLE_IDS_BY_SEX["XY"]],
                metadata=RUN_METADATA
            output:
                ref=REF_OUTPUTS_BY_SEX["XY"]
            log:
                project_path("logs", "wisecondorx", "build_reference_XY.log")
            benchmark:
                BENCH_BUILD_REFERENCE_BY_SEX["XY"]
            params:
                wise=config["biosoft"]["WisecondorX"],
                binsize=WISE_CFG["binsize"],
                converted_dir=project_path("wisecondorx", "converted", "XY"),
                sample_ids=REF_SAMPLE_IDS_BY_SEX["XY"],
                sample_text="\n".join(REF_SAMPLE_IDS_BY_SEX["XY"])
            threads: 4
            run:
                from pgta.core.logging import setup_logger, write_rule_audit_log
                from pgta.reference.tune import build_reference, convert_all_bams

                write_rule_audit_log(log[0], input.metadata, [("SEX GROUP", "XY"), ("REFERENCE SAMPLES", params.sample_text)])
                logger = setup_logger("build_wisecondorx_reference_fixed", log[0])
                npz_paths = convert_all_bams(
                    wisecondorx=params.wise,
                    bams=input.bams,
                    sample_ids=params.sample_ids,
                    binsize=params.binsize,
                    output_dir=Path(params.converted_dir),
                    threads=threads,
                    logger=logger,
                )
                build_reference(
                    wisecondorx=params.wise,
                    binsize=params.binsize,
                    npz_paths=npz_paths,
                    reference_output=Path(output.ref),
                    threads=threads,
                    logger=logger,
                )

        rule build_wisecondorx_gender_reference_fixed:
            input:
                bams=REFERENCE_BAMS,
                metadata=RUN_METADATA
            output:
                ref=GENDER_REF_OUTPUT
            log:
                project_path("logs", "wisecondorx", "build_reference_gender.log")
            benchmark:
                BENCH_BUILD_GENDER_REFERENCE
            params:
                wise=config["biosoft"]["WisecondorX"],
                binsize=WISE_CFG["binsize"],
                converted_dir=project_path("wisecondorx", "converted", "gender"),
                sample_ids=REF_SAMPLE_IDS,
                sample_text=REFERENCE_SAMPLE_TEXT
            threads: 4
            run:
                from pgta.core.logging import setup_logger, write_rule_audit_log
                from pgta.reference.tune import build_reference, convert_all_bams

                write_rule_audit_log(log[0], input.metadata, [("REFERENCE SAMPLES", params.sample_text)])
                logger = setup_logger("build_wisecondorx_gender_reference_fixed", log[0])
                npz_paths = convert_all_bams(
                    wisecondorx=params.wise,
                    bams=input.bams,
                    sample_ids=params.sample_ids,
                    binsize=params.binsize,
                    output_dir=Path(params.converted_dir),
                    threads=threads,
                    logger=logger,
                )
                build_reference(
                    wisecondorx=params.wise,
                    binsize=params.binsize,
                    npz_paths=npz_paths,
                    reference_output=Path(output.ref),
                    threads=threads,
                    logger=logger,
                )
    else:
        rule build_wisecondorx_reference_fixed:
            input:
                bams=REFERENCE_BAMS,
                metadata=RUN_METADATA
            output:
                ref=REF_OUTPUT
            log:
                project_path("logs", "wisecondorx", "build_reference.log")
            benchmark:
                BENCH_BUILD_REFERENCE
            params:
                wise=config["biosoft"]["WisecondorX"],
                binsize=WISE_CFG["binsize"],
                converted_dir=project_path("wisecondorx", "converted"),
                sample_ids=REF_SAMPLE_IDS,
                sample_text=REFERENCE_SAMPLE_TEXT
            threads: 4
            run:
                from pgta.core.logging import setup_logger, write_rule_audit_log
                from pgta.reference.tune import build_reference, convert_all_bams

                write_rule_audit_log(log[0], input.metadata, [("REFERENCE SAMPLES", params.sample_text)])
                logger = setup_logger("build_wisecondorx_reference_fixed", log[0])
                npz_paths = convert_all_bams(
                    wisecondorx=params.wise,
                    bams=input.bams,
                    sample_ids=params.sample_ids,
                    binsize=params.binsize,
                    output_dir=Path(params.converted_dir),
                    threads=threads,
                    logger=logger,
                )
                build_reference(
                    wisecondorx=params.wise,
                    binsize=params.binsize,
                    npz_paths=npz_paths,
                    reference_output=Path(output.ref),
                    threads=threads,
                    logger=logger,
                )
