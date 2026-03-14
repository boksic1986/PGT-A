#!/usr/bin/env python3
import argparse
import re
import shlex
import subprocess
from pathlib import Path

from pipeline_logging import setup_logger


def run_command(command, logger):
    logger.info("$ %s", " ".join(shlex.quote(part) for part in command))
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    if process.stdout is None:
        raise RuntimeError("Failed to capture command stdout.")
    for line in process.stdout:
        line = line.rstrip("\n")
        if not line:
            continue
        logger.info("[cmd] %s", line)
    return_code = process.wait()
    if return_code != 0:
        raise subprocess.CalledProcessError(return_code, command)


def load_best_binsize(best_yaml):
    text = Path(best_yaml).read_text(encoding="utf-8")
    match = re.search(r"^best_binsize:\s*(\d+)\s*$", text, flags=re.MULTILINE)
    if not match:
        raise ValueError(f"best_binsize not found in {best_yaml}")
    return int(match.group(1))


def load_inliers(inlier_file):
    ids = []
    for line in Path(inlier_file).read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line:
            ids.append(line)
    if not ids:
        raise ValueError("No inlier samples found.")
    return ids


def main():
    parser = argparse.ArgumentParser(description="Build WisecondorX reference from tuned binsize and inlier samples.")
    parser.add_argument("--wisecondorx", required=True)
    parser.add_argument("--best-yaml", required=True)
    parser.add_argument("--inlier-samples", required=True)
    parser.add_argument("--allowed-samples", default="", help="Comma-separated sample IDs allowed for this reference.")
    parser.add_argument("--tuning-workdir", required=True)
    parser.add_argument("--reference-output", required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--log", required=True)
    args = parser.parse_args()
    logger = setup_logger("build_reference_from_tuning", args.log)

    binsize = load_best_binsize(args.best_yaml)
    logger.info("loaded best binsize=%s from %s", binsize, args.best_yaml)
    inlier_ids = load_inliers(args.inlier_samples)
    logger.info("loaded inlier sample ids=%d", len(inlier_ids))
    allowed = {item.strip() for item in args.allowed_samples.split(",") if item.strip()}
    if allowed:
        inlier_ids = [sample_id for sample_id in inlier_ids if sample_id in allowed]
        if not inlier_ids:
            raise ValueError("No inlier samples left after applying --allowed-samples filter.")
        logger.info("after --allowed-samples filter, inliers=%d", len(inlier_ids))
    converted_dir = Path(args.tuning_workdir) / f"bin_{binsize}" / "converted"

    npz_paths = [converted_dir / f"{sample_id}.npz" for sample_id in inlier_ids]
    missing = [str(path) for path in npz_paths if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing NPZ files for reference build: {', '.join(missing)}")

    Path(args.reference_output).parent.mkdir(parents=True, exist_ok=True)
    run_command(
        [
            args.wisecondorx,
            "newref",
            *[str(path) for path in npz_paths],
            str(args.reference_output),
            "--binsize",
            str(binsize),
            "--cpus",
            str(args.threads),
        ],
        logger,
    )
    logger.info(
        "reference build completed: output=%s samples=%d binsize=%s",
        args.reference_output,
        len(npz_paths),
        binsize,
    )


if __name__ == "__main__":
    main()
