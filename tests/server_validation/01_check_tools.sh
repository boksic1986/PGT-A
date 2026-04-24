#!/usr/bin/env bash
set -euo pipefail

PYTHON_BIN="/biosoftware/miniconda/envs/snakemake_env/bin/python"
SNAKEMAKE_BIN="/biosoftware/miniconda/envs/snakemake_env/bin/snakemake"
WISECONDORX_BIN="/biosoftware/miniconda/envs/wise_env/bin/WisecondorX"
FASTP_BIN="/biosoftware/bin/fastp"
BWA_BIN="/biosoftware/bin/bwa"
SAMTOOLS_BIN="/biosoftware/bin/samtools"

echo "[check] executables"
for exe in \
    "$PYTHON_BIN" \
    "$SNAKEMAKE_BIN" \
    "$WISECONDORX_BIN" \
    "$FASTP_BIN" \
    "$BWA_BIN" \
    "$SAMTOOLS_BIN"
do
    if [[ ! -x "$exe" ]]; then
        echo "[fail] executable not found or not executable: $exe" >&2
        exit 1
    fi
    echo "[ok] $exe"
done

echo
echo "[check] versions"
"$PYTHON_BIN" --version
"$SNAKEMAKE_BIN" --version
"$WISECONDORX_BIN" --version
"$FASTP_BIN" --version
"$SAMTOOLS_BIN" --version
"$BWA_BIN" 2>&1 | head -n 1

echo
echo "[done] tool validation passed"
