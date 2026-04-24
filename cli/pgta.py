#!/usr/bin/env python3
"""Small CLI scaffold for planning legacy-compatible workflow commands."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


SNAKEMAKE_BIN = "/biosoftware/miniconda/envs/snakemake_env/bin/snakemake"
PYTHON_BIN = "/biosoftware/miniconda/envs/snakemake_env/bin/python"
WISECONDORX_BIN = "/biosoftware/miniconda/envs/wise_env/bin/WisecondorX"

STAGE_TARGETS = {
    "qc": ["mapping", "metadata", "baseline_qc"],
    "reference": ["mapping", "metadata", "reference_qc", "reference"],
    "predict": ["mapping", "metadata", "cnv_qc", "cnv"],
}

DEFAULT_CONFIGS = {
    "qc": "config_qc.yaml",
    "reference": "config_reference.yaml",
    "predict": "config_predict.yaml",
}


def repo_root() -> Path:
    return Path(__file__).resolve().parent.parent


def legacy_snakefile() -> Path:
    return repo_root() / "Snakefile"


def server_validation_script() -> Path:
    return repo_root() / "tests" / "server_validation" / "03_snakemake_dryrun.sh"


def build_snakemake_command(stage: str, config_file: str | None = None, cores: int = 32) -> str:
    if stage not in STAGE_TARGETS:
        raise ValueError(f"Unknown stage: {stage}")
    config_path = config_file or DEFAULT_CONFIGS[stage]
    targets = " ".join(STAGE_TARGETS[stage])
    return (
        f"{SNAKEMAKE_BIN} "
        f"-s {legacy_snakefile()} "
        f"--configfile {repo_root() / config_path} "
        f"--cores {cores} -p {targets}"
    )


def command_layout(_: argparse.Namespace) -> int:
    root = repo_root()
    required = [
        "config",
        "cli",
        "tests",
        "architecture.md",
        "validation_plan.md",
        "Snakefile",
        "build_samples.yaml",
        "predict_samples.yaml",
        "config_qc.yaml",
        "config_reference.yaml",
        "config_predict.yaml",
        "tests/server_validation/03_snakemake_dryrun.sh",
    ]
    missing = [item for item in required if not (root / item).exists()]
    print(f"repo_root: {root}")
    print(f"snakefile: {legacy_snakefile()}")
    print(f"python: {PYTHON_BIN}")
    print(f"wisecondorx: {WISECONDORX_BIN}")
    if missing:
        print("missing:")
        for item in missing:
            print(f"  - {item}")
        return 1
    print("layout_ok: true")
    return 0


def command_plan(args: argparse.Namespace) -> int:
    command = build_snakemake_command(args.stage, args.config, args.cores)
    print("flow: main")
    print(f"stage: {args.stage}")
    print(f"targets: {' '.join(STAGE_TARGETS[args.stage])}")
    print(f"command: {command}")
    return 0


def command_test_plan(args: argparse.Namespace) -> int:
    repo = Path(args.repo).resolve() if args.repo else repo_root()
    cfg_qc = Path(args.config_qc).resolve() if args.config_qc else repo / DEFAULT_CONFIGS["qc"]
    cfg_ref = Path(args.config_reference).resolve() if args.config_reference else repo / DEFAULT_CONFIGS["reference"]
    cfg_pred = Path(args.config_predict).resolve() if args.config_predict else repo / DEFAULT_CONFIGS["predict"]
    print("flow: test")
    print("scope: server_validation")
    print(
        "command: "
        f"bash {server_validation_script()} "
        f"{repo} {cfg_qc} {cfg_ref} {cfg_pred}"
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="PGT-A workflow entry CLI")
    subparsers = parser.add_subparsers(dest="command", required=True)

    layout_parser = subparsers.add_parser("layout", help="Check scaffold layout")
    layout_parser.set_defaults(func=command_layout)

    plan_parser = subparsers.add_parser("plan", help="Print legacy-compatible Snakemake command")
    plan_parser.add_argument("--stage", required=True, choices=sorted(STAGE_TARGETS))
    plan_parser.add_argument("--config", default=None, help="Override config file path")
    plan_parser.add_argument("--cores", type=int, default=32)
    plan_parser.set_defaults(func=command_plan)

    test_plan_parser = subparsers.add_parser("test-plan", help="Print server-side validation command")
    test_plan_parser.add_argument("--repo", default=None, help="Override repository root path")
    test_plan_parser.add_argument("--config-qc", dest="config_qc", default=None, help="Override QC config path")
    test_plan_parser.add_argument("--config-reference", dest="config_reference", default=None, help="Override reference config path")
    test_plan_parser.add_argument("--config-predict", dest="config_predict", default=None, help="Override predict config path")
    test_plan_parser.set_defaults(func=command_test_plan)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
