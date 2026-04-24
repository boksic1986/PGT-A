#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo "Usage: bash tests/server_validation/03_snakemake_dryrun.sh /abs/path/to/PGT-A /abs/path/to/config_qc.yaml /abs/path/to/config_reference.yaml /abs/path/to/config_predict.yaml" >&2
    exit 1
fi

REPO="$1"
CFG_QC="$2"
CFG_REF="$3"
CFG_PRED="$4"
SNAKEMAKE_BIN="/biosoftware/miniconda/envs/snakemake_env/bin/snakemake"

for path_value in "$REPO" "$CFG_QC" "$CFG_REF" "$CFG_PRED"; do
    if [[ ! -e "$path_value" ]]; then
        echo "[fail] required path not found: $path_value" >&2
        exit 1
    fi
done

echo "[check] snakemake --list with QC config"
"$SNAKEMAKE_BIN" \
    -s "$REPO/Snakefile" \
    --configfile "$CFG_QC" \
    --cores 1 \
    --list >/dev/null

echo
echo "[check] dry-run baseline_qc"
"$SNAKEMAKE_BIN" \
    -s "$REPO/Snakefile" \
    --configfile "$CFG_QC" \
    --cores 1 \
    -n -p baseline_qc

echo
echo "[check] dry-run reference_qc + reference"
"$SNAKEMAKE_BIN" \
    -s "$REPO/Snakefile" \
    --configfile "$CFG_REF" \
    --cores 1 \
    -n -p reference_qc reference

echo
echo "[check] dry-run cnv_qc + cnv"
"$SNAKEMAKE_BIN" \
    -s "$REPO/Snakefile" \
    --configfile "$CFG_PRED" \
    --cores 1 \
    -n -p cnv_qc cnv

echo
echo "[done] snakemake dry-run validation passed"
