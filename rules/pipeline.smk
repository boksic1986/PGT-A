if TUNING_ENABLED:
    FINAL_REF_TARGETS = [REF_OUTPUT, TUNING_SUMMARY, TUNING_BEST]
else:
    FINAL_REF_TARGETS = [REF_OUTPUT]


rule all:
    input:
        expand(SORTED_BAM, sample=SAMPLES),
        expand(SORTED_BAI, sample=SAMPLES),
        FINAL_REF_TARGETS


rule fastp_bwa:
    input:
        r1=lambda wildcards: config["samples"][wildcards.sample]["R1"],
        r2=lambda wildcards: config["samples"][wildcards.sample]["R2"]
    output:
        clean_r1=FASTP_R1,
        clean_r2=FASTP_R2,
        html=FASTP_HTML,
        json=FASTP_JSON,
        bam=SORTED_BAM,
        bai=SORTED_BAI
    log:
        fastp=project_path("logs", "fastp", "{sample}.log"),
        bwa=project_path("logs", "bwa", "{sample}.log")
    threads: 8
    params:
        fastp=config["biosoft"]["fastp"],
        bwa=config["biosoft"]["bwa"],
        samtools=config["biosoft"]["samtools"],
        reference=config["core"]["reference_genome"]
    shell:
        r"""
        mkdir -p "$(dirname {output.clean_r1})" "$(dirname {output.bam})" "$(dirname {log.fastp})" "$(dirname {log.bwa})"
        {params.fastp} \
            -i {input.r1} -I {input.r2} \
            -o {output.clean_r1} -O {output.clean_r2} \
            -h {output.html} -j {output.json} \
            -w {threads} \
            > {log.fastp} 2>&1
        {params.bwa} mem -t {threads} {params.reference} {output.clean_r1} {output.clean_r2} 2> {log.bwa} \
            | {params.samtools} sort -@ {threads} -o {output.bam}
        {params.samtools} index {output.bam}
        """


if TUNING_ENABLED:
    rule tune_wisecondorx_bin_pca:
        input:
            bams=expand(SORTED_BAM, sample=SAMPLES)
        output:
            ref=REF_OUTPUT,
            summary=TUNING_SUMMARY,
            best=TUNING_BEST
        log:
            project_path("logs", "wisecondorx", "tuning.log")
        params:
            python_bin=config["biosoft"]["python"],
            wise=config["biosoft"]["WisecondorX"],
            bin_sizes=",".join(str(item) for item in TUNING_BIN_SIZES),
            pca_components=",".join(str(item) for item in TUNING_PCA_COMPONENTS),
            sample_ids=SAMPLES,
            workdir=TUNING_WORKDIR
        threads: 4
        shell:
            r"""
            mkdir -p "$(dirname {output.ref})" "$(dirname {output.summary})" "$(dirname {log})"
            {params.python_bin:q} scripts/tune_wisecondorx_bin_pca.py \
                --wisecondorx {params.wise:q} \
                --bams {input.bams:q} \
                --sample-ids {params.sample_ids:q} \
                --bin-sizes {params.bin_sizes} \
                --pca-components {params.pca_components} \
                --threads {threads} \
                --workdir {params.workdir:q} \
                --summary-output {output.summary:q} \
                --best-output {output.best:q} \
                --reference-output {output.ref:q} \
                --log {log:q}
            """
else:
    rule build_wisecondorx_reference:
        input:
            bams=expand(SORTED_BAM, sample=SAMPLES)
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
            for bam in {input.bams:q}; do
                sample=$(basename "$bam")
                sample=${sample%.sorted.bam}
                npz="{params.converted_dir}/${sample}.npz"
                {params.wise:q} convert "$bam" "$npz" --binsize {params.binsize} >> {log:q} 2>&1
            done
            {params.wise:q} newref {params.converted_dir:q}/*.npz {output.ref:q} --binsize {params.binsize} >> {log:q} 2>&1
            """
