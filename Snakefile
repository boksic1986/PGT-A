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

include: "rules/pipeline.smk"
