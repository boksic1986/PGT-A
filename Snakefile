configfile: "config.yaml"

from pathlib import Path


SAMPLES = sorted(config["samples"].keys())
if not SAMPLES:
    raise ValueError("No samples found in config['samples'].")

PROJECT = Path(config["core"]["project_path"])
WISE_CFG = config["core"]["wisecondorx"]


def project_path(*parts):
    return str(PROJECT.joinpath(*parts))


def resolve_path(path_value):
    path_obj = Path(path_value)
    if path_obj.is_absolute():
        return str(path_obj)
    return str(PROJECT / path_obj)


FASTP_R1 = project_path("fastp", "{sample}.R1.clean.fastq.gz")
FASTP_R2 = project_path("fastp", "{sample}.R2.clean.fastq.gz")
FASTP_HTML = project_path("fastp", "{sample}.fastp.html")
FASTP_JSON = project_path("fastp", "{sample}.fastp.json")
SORTED_BAM = project_path("mapping", "{sample}.sorted.bam")
SORTED_BAI = project_path("mapping", "{sample}.sorted.bam.bai")
REF_OUTPUT = resolve_path(WISE_CFG["reference_output"])


rule all:
    input:
        expand(SORTED_BAM, sample=SAMPLES),
        expand(SORTED_BAI, sample=SAMPLES),
        REF_OUTPUT


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
        for bam in {input.bams}; do
            sample=$(basename "$bam")
            sample=${sample%.sorted.bam}
            npz="{params.converted_dir}/${sample}.npz"
            {params.wise} convert "$bam" "$npz" --binsize {params.binsize} >> {log} 2>&1
        done
        {params.wise} newref {params.converted_dir}/*.npz {output.ref} --binsize {params.binsize} >> {log} 2>&1
        """
