from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path

SCRIPT_MODULE_MAP = {
    "init_pgta_project.py": "pgta.core.init_project",
    "init_predict_config.py": "pgta.core.init_predict",
    "runtime.py": "pgta.core.runtime_tracking",
}

DISPATCHER_MODULE_MAP = {
    "predict.py": {
        "cnv_artifact_rules": "pgta.predict.branch_b.artifact_rules",
        "cnv_benchmark": "pgta.predict.benchmark",
        "cnv_calibration": "pgta.predict.branch_b.calibration",
        "cnv_calling": "pgta.predict.branch_b.calling",
        "cnv_correction": "pgta.predict.branch_b.correction",
        "cnv_evaluation": "pgta.predict.evaluation",
        "cnv_ml": "pgta.predict.ml",
        "cnv_mosaic_fraction": "pgta.predict.branch_b.mosaic_fraction",
        "cnv_qc": "pgta.predict.cnv_qc",
        "cnv_report": "pgta.predict.report",
        "wisecondorx_gender": "pgta.predict.sex_routing",
    },
    "qc.py": {
        "aggregate_baseline_qc": "pgta.qc.aggregate",
        "bam_uniformity_qc": "pgta.qc.sample_qc",
        "baseline_qc_report": "pgta.qc.report",
        "batch_bam_uniformity_qc": "pgta.qc.batch_entry",
        "run_qc_framework": "pgta.qc.framework",
    },
    "reference.py": {
        "build_bin_annotations": "pgta.reference.assets",
        "build_frozen_reference_package": "pgta.reference.package",
        "build_reference_from_tuning": "pgta.reference.build",
        "reference_prefilter_qc": "pgta.reference.prefilter",
        "select_reference_cohorts": "pgta.reference.cohort",
        "tune_wisecondorx_bin_pca": "pgta.reference.tune",
    },
    "validation.py": {
        "compare_outputs": "pgta.validation.compare_outputs",
        "validate_fraction_truth": "pgta.validation.validate_fraction_truth",
    },
}


def bootstrap_repo_root() -> Path:
    root = Path(__file__).resolve().parents[1]
    root_str = str(root)
    if root_str not in sys.path:
        sys.path.insert(0, root_str)
    return root


def export_module(namespace: dict[str, object], module_name: str):
    bootstrap_repo_root()
    module = importlib.import_module(module_name)
    if "snakemake" in namespace:
        setattr(module, "snakemake", namespace["snakemake"])
    public_names = getattr(module, "__all__", None)
    if public_names is None:
        public_names = [name for name in dir(module) if not name.startswith("_")]
    for name in public_names:
        namespace[name] = getattr(module, name)
    return getattr(module, "main", None)


def export_script(namespace: dict[str, object], script_path: str | Path):
    script_name = Path(script_path).name
    try:
        module_name = SCRIPT_MODULE_MAP[script_name]
    except KeyError as exc:
        raise KeyError(
            f"no compatibility module mapping configured for {script_name}"
        ) from exc
    return export_module(namespace, module_name)


def run_dispatcher(script_path: str | Path, argv: list[str] | None = None):
    bootstrap_repo_root()
    script_name = Path(script_path).name
    try:
        action_map = DISPATCHER_MODULE_MAP[script_name]
    except KeyError as exc:
        raise KeyError(
            f"no dispatcher action mapping configured for {script_name}"
        ) from exc

    parser = argparse.ArgumentParser(
        description=f"Grouped compatibility dispatcher for {script_name}."
    )
    parser.add_argument(
        "action",
        nargs="?",
        choices=sorted(action_map),
        help="Action to execute.",
    )
    parsed, remaining = parser.parse_known_args(list(sys.argv[1:] if argv is None else argv))
    if not parsed.action:
        parser.print_help()
        return 0

    module = importlib.import_module(action_map[parsed.action])
    main = getattr(module, "main", None)
    if main is None:
        raise AttributeError(f"module {module.__name__} does not expose main()")

    original_argv = sys.argv[:]
    try:
        sys.argv = [f"{script_name} {parsed.action}", *remaining]
        return main()
    finally:
        sys.argv = original_argv
