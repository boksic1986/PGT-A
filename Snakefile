configfile: "config.yaml"

from pathlib import Path


SAMPLES = sorted(config["samples"].keys())
if not SAMPLES:
    raise ValueError("No samples found in config['samples'].")

PROJECT = Path(config["core"]["project_path"])
WISE_CFG = config["core"]["wisecondorx"]
TUNE_CFG = WISE_CFG.get("tuning", {})


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
TUNING_ENABLED = bool(TUNE_CFG.get("enable", False))
TUNING_BIN_SIZES = [int(item) for item in TUNE_CFG.get("bin_sizes", [int(WISE_CFG["binsize"])])]
TUNING_PCA_COMPONENTS = [int(item) for item in TUNE_CFG.get("pca_components", [2, 3, 4, 5])]
TUNING_WORKDIR = project_path("wisecondorx", "tuning")
TUNING_SUMMARY = project_path("wisecondorx", "tuning", "bin_pca_grid.tsv")
TUNING_BEST = project_path("wisecondorx", "tuning", "best_params.yaml")

include: "rules/pipeline.smk"
