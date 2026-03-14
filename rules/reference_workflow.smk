if TUNING_ENABLED:
    pass
    if "XX" in REF_SEXES:
        pass
        rule reference_prefilter_xx:
            input:
                bams=[SORTED_BAM.format(sample=s) for s in REF_SAMPLE_IDS_BY_SEX["XX"]],
                metadata=RUN_METADATA
            output:
                qc=str(Path(REF_PREFILTER_DIR_BY_SEX["XX"]) / "reference_sample_qc.tsv"),
                plot=str(Path(REF_PREFILTER_DIR_BY_SEX["XX"]) / "reference_sample_qc.svg"),
                inliers=str(Path(REF_PREFILTER_DIR_BY_SEX["XX"]) / "reference_inlier_samples.txt"),
                summary=str(Path(REF_PREFILTER_DIR_BY_SEX["XX"]) / "prefilter_summary.yaml")
            log:
                project_path("logs", "wisecondorx", "reference_prefilter_XX.log")
            params:
                python_bin=config["biosoft"]["python"],
                wise=config["biosoft"]["WisecondorX"],
                sample_ids=REF_SAMPLE_IDS_BY_SEX["XX"],
                binsize=PREFILTER_BINSIZE,
                pca_min=TUNING_PCA_MIN,
                pca_max=TUNING_PCA_MAX,
                min_ref_samples=TUNING_MIN_REF_SAMPLES,
                min_reads=TUNING_MIN_READS,
                min_corr=TUNING_MIN_CORR,
                max_recon_z=TUNING_MAX_RECON_Z,
                max_noise_z=TUNING_MAX_NOISE_Z,
                max_iter=PREFILTER_MAX_ITER,
                workdir=REF_PREFILTER_DIR_BY_SEX["XX"]
            threads: 4
            shell:
                r"""
                mkdir -p "{params.workdir}" "$(dirname {log})"
                (
                    echo "=== PIPELINE AUDIT ==="
                    cat {input.metadata:q}
                    echo "=== SEX GROUP ==="
                    echo "XX"
                    echo "=== COMMAND ==="
                ) > {log:q}
                {params.python_bin:q} {SCRIPT_REFERENCE_PREFILTER_QC:q} \
                    --wisecondorx {params.wise:q} \
                    --bams {input.bams:q} \
                    --sample-ids {params.sample_ids:q} \
                    --binsize {params.binsize} \
                    --pca-min-components {params.pca_min} \
                    --pca-max-components {params.pca_max} \
                    --min-reference-samples {params.min_ref_samples} \
                    --min-reads-per-sample {params.min_reads} \
                    --min-corr-to-median {params.min_corr} \
                    --max-reconstruction-error-z {params.max_recon_z} \
                    --max-noise-mad-z {params.max_noise_z} \
                    --max-iterations {params.max_iter} \
                    --threads {threads} \
                    --workdir {params.workdir:q} \
                    --qc-output {output.qc:q} \
                    --plot-output {output.plot:q} \
                    --inlier-samples-output {output.inliers:q} \
                    --summary-output {output.summary:q} \
                    --log {log:q}
                """

        rule tune_wisecondorx_reference_qc_xx:
            input:
                bams=[SORTED_BAM.format(sample=s) for s in REF_SAMPLE_IDS_BY_SEX["XX"]],
                prefilter_inliers=str(Path(REF_PREFILTER_DIR_BY_SEX["XX"]) / "reference_inlier_samples.txt"),
                metadata=RUN_METADATA
            output:
                summary=str(Path(REF_TUNING_DIR_BY_SEX["XX"]) / "bin_pca_grid.tsv"),
                best=str(Path(REF_TUNING_DIR_BY_SEX["XX"]) / "best_params.yaml"),
                qc=str(Path(REF_TUNING_DIR_BY_SEX["XX"]) / "reference_sample_qc.tsv"),
                plot=str(Path(REF_TUNING_DIR_BY_SEX["XX"]) / "best_bin_pca_elbow.svg"),
                qc_plot=str(Path(REF_TUNING_DIR_BY_SEX["XX"]) / "reference_qc_metrics.svg"),
                inliers=str(Path(REF_TUNING_DIR_BY_SEX["XX"]) / "reference_inlier_samples.txt")
            log:
                project_path("logs", "wisecondorx", "tuning_XX.log")
            params:
                python_bin=config["biosoft"]["python"],
                wise=config["biosoft"]["WisecondorX"],
                bin_sizes=",".join(str(item) for item in TUNING_BIN_SIZES),
                sample_ids=REF_SAMPLE_IDS_BY_SEX["XX"],
                pca_min=TUNING_PCA_MIN,
                pca_max=TUNING_PCA_MAX,
                pca_min_var=TUNING_MIN_VAR,
                min_ref_samples=TUNING_MIN_REF_SAMPLES,
                max_outlier_fraction=TUNING_MAX_OUTLIER_FRAC,
                min_reads=TUNING_MIN_READS,
                min_corr=TUNING_MIN_CORR,
                max_recon_z=TUNING_MAX_RECON_Z,
                max_noise_z=TUNING_MAX_NOISE_Z,
                workdir=REF_TUNING_DIR_BY_SEX["XX"],
                reference_output=REF_OUTPUTS_BY_SEX["XX"]
            threads: 4
            shell:
                r"""
                mkdir -p "{params.workdir}" "$(dirname {log})"
                (
                    echo "=== PIPELINE AUDIT ==="
                    cat {input.metadata:q}
                    echo "=== SEX GROUP ==="
                    echo "XX"
                    echo "=== COMMAND ==="
                ) > {log:q}
                {params.python_bin:q} {SCRIPT_TUNE_WISECONDORX:q} \
                    --wisecondorx {params.wise:q} \
                    --bams {input.bams:q} \
                    --sample-ids {params.sample_ids:q} \
                    --allowed-samples-file {input.prefilter_inliers:q} \
                    --bin-sizes {params.bin_sizes} \
                    --pca-min-components {params.pca_min} \
                    --pca-max-components {params.pca_max} \
                    --pca-min-explained-variance {params.pca_min_var} \
                    --min-reference-samples {params.min_ref_samples} \
                    --max-outlier-fraction {params.max_outlier_fraction} \
                    --min-reads-per-sample {params.min_reads} \
                    --min-corr-to-median {params.min_corr} \
                    --max-reconstruction-error-z {params.max_recon_z} \
                    --max-noise-mad-z {params.max_noise_z} \
                    --threads {threads} \
                    --workdir {params.workdir:q} \
                    --summary-output {output.summary:q} \
                    --best-output {output.best:q} \
                    --qc-output {output.qc:q} \
                    --plot-output {output.plot:q} \
                    --qc-stats-plot-output {output.qc_plot:q} \
                    --inlier-samples-output {output.inliers:q} \
                    --reference-output {params.reference_output:q} \
                    --skip-build-reference \
                    --log {log:q}
                """

        rule build_wisecondorx_reference_from_tuning_xx:
            input:
                best=str(Path(REF_TUNING_DIR_BY_SEX["XX"]) / "best_params.yaml"),
                inliers=str(Path(REF_TUNING_DIR_BY_SEX["XX"]) / "reference_inlier_samples.txt"),
                metadata=RUN_METADATA
            output:
                ref=REF_OUTPUTS_BY_SEX["XX"]
            log:
                project_path("logs", "wisecondorx", "build_reference_XX.log")
            params:
                python_bin=config["biosoft"]["python"],
                wise=config["biosoft"]["WisecondorX"],
                workdir=REF_TUNING_DIR_BY_SEX["XX"],
                allowed_samples=",".join(REF_SAMPLE_IDS_BY_SEX["XX"])
            threads: 4
            shell:
                r"""
                mkdir -p "$(dirname {output.ref})" "$(dirname {log})"
                (
                    echo "=== PIPELINE AUDIT ==="
                    cat {input.metadata:q}
                    echo "=== SEX GROUP ==="
                    echo "XX"
                    echo "=== COMMAND ==="
                ) > {log:q}
                {params.python_bin:q} {SCRIPT_BUILD_REF_FROM_TUNING:q} \
                    --wisecondorx {params.wise:q} \
                    --best-yaml {input.best:q} \
                    --inlier-samples {input.inliers:q} \
                    --allowed-samples {params.allowed_samples:q} \
                    --tuning-workdir {params.workdir:q} \
                    --reference-output {output.ref:q} \
                    --threads {threads} \
                    --log {log:q}
                """

    if "XY" in REF_SEXES:
        pass
        rule reference_prefilter_xy:
            input:
                bams=[SORTED_BAM.format(sample=s) for s in REF_SAMPLE_IDS_BY_SEX["XY"]],
                metadata=RUN_METADATA
            output:
                qc=str(Path(REF_PREFILTER_DIR_BY_SEX["XY"]) / "reference_sample_qc.tsv"),
                plot=str(Path(REF_PREFILTER_DIR_BY_SEX["XY"]) / "reference_sample_qc.svg"),
                inliers=str(Path(REF_PREFILTER_DIR_BY_SEX["XY"]) / "reference_inlier_samples.txt"),
                summary=str(Path(REF_PREFILTER_DIR_BY_SEX["XY"]) / "prefilter_summary.yaml")
            log:
                project_path("logs", "wisecondorx", "reference_prefilter_XY.log")
            params:
                python_bin=config["biosoft"]["python"],
                wise=config["biosoft"]["WisecondorX"],
                sample_ids=REF_SAMPLE_IDS_BY_SEX["XY"],
                binsize=PREFILTER_BINSIZE,
                pca_min=TUNING_PCA_MIN,
                pca_max=TUNING_PCA_MAX,
                min_ref_samples=TUNING_MIN_REF_SAMPLES,
                min_reads=TUNING_MIN_READS,
                min_corr=TUNING_MIN_CORR,
                max_recon_z=TUNING_MAX_RECON_Z,
                max_noise_z=TUNING_MAX_NOISE_Z,
                max_iter=PREFILTER_MAX_ITER,
                workdir=REF_PREFILTER_DIR_BY_SEX["XY"]
            threads: 4
            shell:
                r"""
                mkdir -p "{params.workdir}" "$(dirname {log})"
                (
                    echo "=== PIPELINE AUDIT ==="
                    cat {input.metadata:q}
                    echo "=== SEX GROUP ==="
                    echo "XY"
                    echo "=== COMMAND ==="
                ) > {log:q}
                {params.python_bin:q} {SCRIPT_REFERENCE_PREFILTER_QC:q} \
                    --wisecondorx {params.wise:q} \
                    --bams {input.bams:q} \
                    --sample-ids {params.sample_ids:q} \
                    --binsize {params.binsize} \
                    --pca-min-components {params.pca_min} \
                    --pca-max-components {params.pca_max} \
                    --min-reference-samples {params.min_ref_samples} \
                    --min-reads-per-sample {params.min_reads} \
                    --min-corr-to-median {params.min_corr} \
                    --max-reconstruction-error-z {params.max_recon_z} \
                    --max-noise-mad-z {params.max_noise_z} \
                    --max-iterations {params.max_iter} \
                    --threads {threads} \
                    --workdir {params.workdir:q} \
                    --qc-output {output.qc:q} \
                    --plot-output {output.plot:q} \
                    --inlier-samples-output {output.inliers:q} \
                    --summary-output {output.summary:q} \
                    --log {log:q}
                """

        rule tune_wisecondorx_reference_qc_xy:
            input:
                bams=[SORTED_BAM.format(sample=s) for s in REF_SAMPLE_IDS_BY_SEX["XY"]],
                prefilter_inliers=str(Path(REF_PREFILTER_DIR_BY_SEX["XY"]) / "reference_inlier_samples.txt"),
                metadata=RUN_METADATA
            output:
                summary=str(Path(REF_TUNING_DIR_BY_SEX["XY"]) / "bin_pca_grid.tsv"),
                best=str(Path(REF_TUNING_DIR_BY_SEX["XY"]) / "best_params.yaml"),
                qc=str(Path(REF_TUNING_DIR_BY_SEX["XY"]) / "reference_sample_qc.tsv"),
                plot=str(Path(REF_TUNING_DIR_BY_SEX["XY"]) / "best_bin_pca_elbow.svg"),
                qc_plot=str(Path(REF_TUNING_DIR_BY_SEX["XY"]) / "reference_qc_metrics.svg"),
                inliers=str(Path(REF_TUNING_DIR_BY_SEX["XY"]) / "reference_inlier_samples.txt")
            log:
                project_path("logs", "wisecondorx", "tuning_XY.log")
            params:
                python_bin=config["biosoft"]["python"],
                wise=config["biosoft"]["WisecondorX"],
                bin_sizes=",".join(str(item) for item in TUNING_BIN_SIZES),
                sample_ids=REF_SAMPLE_IDS_BY_SEX["XY"],
                pca_min=TUNING_PCA_MIN,
                pca_max=TUNING_PCA_MAX,
                pca_min_var=TUNING_MIN_VAR,
                min_ref_samples=TUNING_MIN_REF_SAMPLES,
                max_outlier_fraction=TUNING_MAX_OUTLIER_FRAC,
                min_reads=TUNING_MIN_READS,
                min_corr=TUNING_MIN_CORR,
                max_recon_z=TUNING_MAX_RECON_Z,
                max_noise_z=TUNING_MAX_NOISE_Z,
                workdir=REF_TUNING_DIR_BY_SEX["XY"],
                reference_output=REF_OUTPUTS_BY_SEX["XY"]
            threads: 4
            shell:
                r"""
                mkdir -p "{params.workdir}" "$(dirname {log})"
                (
                    echo "=== PIPELINE AUDIT ==="
                    cat {input.metadata:q}
                    echo "=== SEX GROUP ==="
                    echo "XY"
                    echo "=== COMMAND ==="
                ) > {log:q}
                {params.python_bin:q} {SCRIPT_TUNE_WISECONDORX:q} \
                    --wisecondorx {params.wise:q} \
                    --bams {input.bams:q} \
                    --sample-ids {params.sample_ids:q} \
                    --allowed-samples-file {input.prefilter_inliers:q} \
                    --bin-sizes {params.bin_sizes} \
                    --pca-min-components {params.pca_min} \
                    --pca-max-components {params.pca_max} \
                    --pca-min-explained-variance {params.pca_min_var} \
                    --min-reference-samples {params.min_ref_samples} \
                    --max-outlier-fraction {params.max_outlier_fraction} \
                    --min-reads-per-sample {params.min_reads} \
                    --min-corr-to-median {params.min_corr} \
                    --max-reconstruction-error-z {params.max_recon_z} \
                    --max-noise-mad-z {params.max_noise_z} \
                    --threads {threads} \
                    --workdir {params.workdir:q} \
                    --summary-output {output.summary:q} \
                    --best-output {output.best:q} \
                    --qc-output {output.qc:q} \
                    --plot-output {output.plot:q} \
                    --qc-stats-plot-output {output.qc_plot:q} \
                    --inlier-samples-output {output.inliers:q} \
                    --reference-output {params.reference_output:q} \
                    --skip-build-reference \
                    --log {log:q}
                """

        rule build_wisecondorx_reference_from_tuning_xy:
            input:
                best=str(Path(REF_TUNING_DIR_BY_SEX["XY"]) / "best_params.yaml"),
                inliers=str(Path(REF_TUNING_DIR_BY_SEX["XY"]) / "reference_inlier_samples.txt"),
                metadata=RUN_METADATA
            output:
                ref=REF_OUTPUTS_BY_SEX["XY"]
            log:
                project_path("logs", "wisecondorx", "build_reference_XY.log")
            params:
                python_bin=config["biosoft"]["python"],
                wise=config["biosoft"]["WisecondorX"],
                workdir=REF_TUNING_DIR_BY_SEX["XY"],
                allowed_samples=",".join(REF_SAMPLE_IDS_BY_SEX["XY"])
            threads: 4
            shell:
                r"""
                mkdir -p "$(dirname {output.ref})" "$(dirname {log})"
                (
                    echo "=== PIPELINE AUDIT ==="
                    cat {input.metadata:q}
                    echo "=== SEX GROUP ==="
                    echo "XY"
                    echo "=== COMMAND ==="
                ) > {log:q}
                {params.python_bin:q} {SCRIPT_BUILD_REF_FROM_TUNING:q} \
                    --wisecondorx {params.wise:q} \
                    --best-yaml {input.best:q} \
                    --inlier-samples {input.inliers:q} \
                    --allowed-samples {params.allowed_samples:q} \
                    --tuning-workdir {params.workdir:q} \
                    --reference-output {output.ref:q} \
                    --threads {threads} \
                    --log {log:q}
                """

    if not REF_SEXES:
        pass
        rule tune_wisecondorx_reference_qc:
            input:
                bams=expand(SORTED_BAM, sample=SAMPLES),
                metadata=RUN_METADATA
            output:
                summary=TUNING_SUMMARY,
                best=TUNING_BEST,
                qc=TUNING_QC,
                plot=TUNING_PLOT,
                qc_plot=TUNING_QC_STATS_PLOT,
                inliers=TUNING_INLIERS
            log:
                project_path("logs", "wisecondorx", "tuning.log")
            params:
                python_bin=config["biosoft"]["python"],
                wise=config["biosoft"]["WisecondorX"],
                bin_sizes=",".join(str(item) for item in TUNING_BIN_SIZES),
                sample_ids=SAMPLES,
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
                reference_output=REF_OUTPUT
            threads: 4
            shell:
                r"""
                mkdir -p "$(dirname {output.summary})" "$(dirname {log})"
                (
                    echo "=== PIPELINE AUDIT ==="
                    cat {input.metadata:q}
                    echo "=== COMMAND ==="
                ) > {log:q}
                {params.python_bin:q} {SCRIPT_TUNE_WISECONDORX:q} \
                    --wisecondorx {params.wise:q} \
                    --bams {input.bams:q} \
                    --sample-ids {params.sample_ids:q} \
                    --bin-sizes {params.bin_sizes} \
                    --pca-min-components {params.pca_min} \
                    --pca-max-components {params.pca_max} \
                    --pca-min-explained-variance {params.pca_min_var} \
                    --min-reference-samples {params.min_ref_samples} \
                    --max-outlier-fraction {params.max_outlier_fraction} \
                    --min-reads-per-sample {params.min_reads} \
                    --min-corr-to-median {params.min_corr} \
                    --max-reconstruction-error-z {params.max_recon_z} \
                    --max-noise-mad-z {params.max_noise_z} \
                    --threads {threads} \
                    --workdir {params.workdir:q} \
                    --summary-output {output.summary:q} \
                    --best-output {output.best:q} \
                    --qc-output {output.qc:q} \
                    --plot-output {output.plot:q} \
                    --qc-stats-plot-output {output.qc_plot:q} \
                    --inlier-samples-output {output.inliers:q} \
                    --reference-output {params.reference_output:q} \
                    --skip-build-reference \
                    --log {log:q}
                """

        rule build_wisecondorx_reference_from_tuning:
            input:
                best=TUNING_BEST,
                inliers=TUNING_INLIERS,
                metadata=RUN_METADATA
            output:
                ref=REF_OUTPUT
            log:
                project_path("logs", "wisecondorx", "build_reference.log")
            params:
                python_bin=config["biosoft"]["python"],
                wise=config["biosoft"]["WisecondorX"],
                workdir=TUNING_WORKDIR
            threads: 4
            shell:
                r"""
                mkdir -p "$(dirname {output.ref})" "$(dirname {log})"
                (
                    echo "=== PIPELINE AUDIT ==="
                    cat {input.metadata:q}
                    echo "=== COMMAND ==="
                ) > {log:q}
                {params.python_bin:q} {SCRIPT_BUILD_REF_FROM_TUNING:q} \
                    --wisecondorx {params.wise:q} \
                    --best-yaml {input.best:q} \
                    --inlier-samples {input.inliers:q} \
                    --tuning-workdir {params.workdir:q} \
                    --reference-output {output.ref:q} \
                    --threads {threads} \
                    --log {log:q}
                """
else:
    pass
    if "XX" in REF_SEXES:
        pass
        rule build_wisecondorx_reference_fixed_xx:
            input:
                bams=[SORTED_BAM.format(sample=s) for s in REF_SAMPLE_IDS_BY_SEX["XX"]],
                metadata=RUN_METADATA
            output:
                ref=REF_OUTPUTS_BY_SEX["XX"]
            log:
                project_path("logs", "wisecondorx", "build_reference_XX.log")
            params:
                wise=config["biosoft"]["WisecondorX"],
                binsize=WISE_CFG["binsize"],
                converted_dir=project_path("wisecondorx", "converted", "XX")
            threads: 4
            shell:
                r"""
                mkdir -p "{params.converted_dir}" "$(dirname {output.ref})" "$(dirname {log})"
                (
                    echo "=== PIPELINE AUDIT ==="
                    cat {input.metadata:q}
                    echo "=== SEX GROUP ==="
                    echo "XX"
                    echo "=== COMMAND ==="
                ) > {log:q}
                for bam in {input.bams:q}; do
                    sample=$(basename "$bam" .sorted.bam)
                    npz="{params.converted_dir}/${{sample}}.npz"
                    {params.wise:q} convert "$bam" "$npz" --binsize {params.binsize} >> {log:q} 2>&1
                done
                {params.wise:q} newref {params.converted_dir:q}/*.npz {output.ref:q} --binsize {params.binsize} --cpus {threads} >> {log:q} 2>&1
                """

    if "XY" in REF_SEXES:
        pass
        rule build_wisecondorx_reference_fixed_xy:
            input:
                bams=[SORTED_BAM.format(sample=s) for s in REF_SAMPLE_IDS_BY_SEX["XY"]],
                metadata=RUN_METADATA
            output:
                ref=REF_OUTPUTS_BY_SEX["XY"]
            log:
                project_path("logs", "wisecondorx", "build_reference_XY.log")
            params:
                wise=config["biosoft"]["WisecondorX"],
                binsize=WISE_CFG["binsize"],
                converted_dir=project_path("wisecondorx", "converted", "XY")
            threads: 4
            shell:
                r"""
                mkdir -p "{params.converted_dir}" "$(dirname {output.ref})" "$(dirname {log})"
                (
                    echo "=== PIPELINE AUDIT ==="
                    cat {input.metadata:q}
                    echo "=== SEX GROUP ==="
                    echo "XY"
                    echo "=== COMMAND ==="
                ) > {log:q}
                for bam in {input.bams:q}; do
                    sample=$(basename "$bam" .sorted.bam)
                    npz="{params.converted_dir}/${{sample}}.npz"
                    {params.wise:q} convert "$bam" "$npz" --binsize {params.binsize} >> {log:q} 2>&1
                done
                {params.wise:q} newref {params.converted_dir:q}/*.npz {output.ref:q} --binsize {params.binsize} --cpus {threads} >> {log:q} 2>&1
                """
    if not REF_SEXES:
        pass
        rule build_wisecondorx_reference_fixed:
            input:
                bams=expand(SORTED_BAM, sample=SAMPLES),
                metadata=RUN_METADATA
            output:
                ref=REF_OUTPUT
            log:
                project_path("logs", "wisecondorx", "build_reference.log")
            params:
                wise=config["biosoft"]["WisecondorX"],
                binsize=WISE_CFG["binsize"],
                converted_dir=project_path("wisecondorx", "converted")
            threads: 4
            shell:
                r"""
                mkdir -p "{params.converted_dir}" "$(dirname {output.ref})" "$(dirname {log})"
                (
                    echo "=== PIPELINE AUDIT ==="
                    cat {input.metadata:q}
                    echo "=== COMMAND ==="
                ) > {log:q}
                for bam in {input.bams:q}; do
                    sample=$(basename "$bam" .sorted.bam)
                    npz="{params.converted_dir}/${{sample}}.npz"
                    {params.wise:q} convert "$bam" "$npz" --binsize {params.binsize} >> {log:q} 2>&1
                done
                {params.wise:q} newref {params.converted_dir:q}/*.npz {output.ref:q} --binsize {params.binsize} --cpus {threads} >> {log:q} 2>&1
                """
