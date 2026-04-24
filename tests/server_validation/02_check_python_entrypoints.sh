#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: bash tests/server_validation/02_check_python_entrypoints.sh /abs/path/to/PGT-A" >&2
    exit 1
fi

REPO="$1"
PYTHON_BIN="/biosoftware/miniconda/envs/snakemake_env/bin/python"

if [[ ! -d "$REPO" ]]; then
    echo "[fail] repo directory not found: $REPO" >&2
    exit 1
fi

echo "[check] py_compile"
"$PYTHON_BIN" -m py_compile \
    "$REPO/scripts/init_pgta_project.py" \
    "$REPO/scripts/init_predict_config.py" \
    "$REPO/scripts/predict.py" \
    "$REPO/scripts/qc.py" \
    "$REPO/scripts/reference.py" \
    "$REPO/scripts/runtime.py" \
    "$REPO/scripts/validation.py"

echo
echo "[check] --help entrypoints"
"$PYTHON_BIN" "$REPO/scripts/init_pgta_project.py" --help >/dev/null
"$PYTHON_BIN" "$REPO/scripts/init_predict_config.py" --help >/dev/null
"$PYTHON_BIN" "$REPO/scripts/predict.py" --help >/dev/null
"$PYTHON_BIN" "$REPO/scripts/predict.py" cnv_report --help >/dev/null
"$PYTHON_BIN" "$REPO/scripts/qc.py" --help >/dev/null
"$PYTHON_BIN" "$REPO/scripts/qc.py" aggregate_baseline_qc --help >/dev/null
"$PYTHON_BIN" "$REPO/scripts/reference.py" --help >/dev/null
"$PYTHON_BIN" "$REPO/scripts/reference.py" build_bin_annotations --help >/dev/null
"$PYTHON_BIN" "$REPO/scripts/validation.py" --help >/dev/null
"$PYTHON_BIN" "$REPO/scripts/validation.py" validate_fraction_truth --help >/dev/null

echo
echo "[check] compatibility mapping"
"$PYTHON_BIN" "$REPO/tests/unit/test_script_compat.py"

echo
echo "[check] rule import boundary"
"$PYTHON_BIN" "$REPO/tests/unit/test_rule_import_boundary.py"

echo
echo "[done] python entrypoint validation passed"
