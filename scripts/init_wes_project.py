#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author: jiucheng
# Updated: 2026/03/08
import argparse
import shutil
from collections import defaultdict
from pathlib import Path

import yaml


def detect_samples(fq_dir: Path) -> dict:
    fq_dir = Path(fq_dir)
    fq_files = list(fq_dir.rglob("*.fastq.gz")) + list(fq_dir.rglob("*.fq.gz"))

    temp_dict = defaultdict(list)

    for fq in fq_files:
        fname = fq.name

        if "_R1" in fname or "_1" in fname:
            direction = "R1"
        elif "_R2" in fname or "_2" in fname:
            direction = "R2"
        else:
            continue

        sample_id = fname.split("_")[0]
        temp_dict[sample_id].append((direction, fq))

    sample_dict = {}
    for sample_id, items in temp_dict.items():
        r1_list = sorted([fq for d, fq in items if d == "R1"])
        r2_list = sorted([fq for d, fq in items if d == "R2"])
        count = min(len(r1_list), len(r2_list))

        for i in range(count):
            sample_key = f"{sample_id}_{i + 1}" if count > 1 else sample_id
            sample_dict[sample_key] = {
                "R1": r1_list[i],
                "R2": r2_list[i],
            }

    return sample_dict


def load_template_config() -> dict:
    return {
        "core": {
            "project_path": "results",
            "reference_genome": "/data/Database/index/hg19/hg19.fa",
            "chromosome_list": [
                "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
                "chr21", "chr22", "chrX", "chrY",
            ],
            "wisecondorx": {
                "binsize": 500000,
                "reference_output": "wisecondorx/reference/ref_500kb.npz",
                "use_chr_prefix": True,
            },
        },
        "biosoft": {
            "fastp": "/biosoftware/bin/fastp",
            "bwa": "/biosoftware/bin/bwa",
            "samtools": "/biosoftware/bin/samtools",
            "WisecondorX": "/biosoftware/miniconda/envs/wise_env/bin/WisecondorX",
            "python": "/biosoftware/miniconda/envs/wise_env/bin/python",
        },
        "samples": {},
    }


def main(project, fq_dir, output_config, template_snakefile=None):
    fq_dir = Path(fq_dir).resolve()
    project = Path(project).resolve()
    config_path = Path(output_config).resolve()

    project.mkdir(parents=True, exist_ok=True)

    # 1) Create only directories required by fastp and WisecondorX reference building
    for sub in [
        "data/fastq",
        "fastp",
        "mapping",
        "wisecondorx/converted",
        "wisecondorx/reference",
        "logs/fastp",
        "logs/bwa",
        "logs/wisecondorx",
    ]:
        (project / sub).mkdir(parents=True, exist_ok=True)

    # 2) Detect paired-end samples
    samples = detect_samples(fq_dir)

    # 3) Symlink source FASTQs into project/data/fastq and write absolute paths to config
    for sample, fqs in samples.items():
        for direction in ["R1", "R2"]:
            src = fqs[direction].resolve()
            dst = project / "data/fastq" / src.name

            if dst.exists() or dst.is_symlink():
                dst.unlink()
            dst.symlink_to(src)
            samples[sample][direction] = str(dst.resolve())

    # 4) Merge config
    config = load_template_config()
    config["core"]["project_path"] = str(project)
    config["samples"] = samples

    # 5) Write config
    config_path.parent.mkdir(parents=True, exist_ok=True)
    with open(config_path, "w", encoding="utf-8") as f:
        yaml.safe_dump(config, f, sort_keys=False, allow_unicode=True)
    print(f"Config written to {config_path}")

    # 6) Optionally copy Snakefile template
    if template_snakefile:
        snakefile_src = Path(template_snakefile).resolve()
        snakefile_dst = project / "Snakefile"
        if snakefile_src.exists():
            shutil.copy(snakefile_src, snakefile_dst)
            print(f"Snakefile copied to {snakefile_dst}")
        else:
            print(f"Warning: Snakefile not found at {snakefile_src}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Init PGT-A project config (fastp + WisecondorX reference only)"
    )
    parser.add_argument("--project", required=True, help="Path to the new project directory")
    parser.add_argument("--fq_dir", required=True, help="Path to input fastq directory")
    parser.add_argument("--output_config", required=True, help="Path to output config.yaml")
    parser.add_argument(
        "--template_snakefile",
        required=False,
        help="Optional path to template Snakefile",
    )
    args = parser.parse_args()

    main(args.project, args.fq_dir, args.output_config, args.template_snakefile)
