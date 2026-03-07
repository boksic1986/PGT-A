#!/usr/bin/env python3
import argparse
import shlex
import subprocess
from pathlib import Path

import numpy as np


def parse_int_list(value):
    values = []
    for item in str(value).split(","):
        item = item.strip()
        if not item:
            continue
        values.append(int(item))
    if not values:
        raise ValueError("Empty integer list.")
    return values


def run_command(command, log_path):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "a", encoding="utf-8") as log_handle:
        log_handle.write("$ " + " ".join(shlex.quote(part) for part in command) + "\n")
        log_handle.flush()
        subprocess.run(command, check=True, stdout=log_handle, stderr=log_handle)


def load_numeric_arrays(npz_path):
    arrays = {}
    with np.load(npz_path, allow_pickle=True) as data:
        for key in data.files:
            raw = np.asarray(data[key])
            if raw.size < 100 or not np.issubdtype(raw.dtype, np.number):
                continue
            flat = np.nan_to_num(raw.astype(np.float64).reshape(-1), nan=0.0, posinf=0.0, neginf=0.0)
            if flat.size < 100:
                continue
            arrays[key] = flat
    if not arrays:
        raise ValueError(f"No usable numeric arrays found in {npz_path}")
    return arrays


def build_matrix(npz_paths):
    loaded = [load_numeric_arrays(path) for path in npz_paths]
    common_keys = set(loaded[0].keys())
    for item in loaded[1:]:
        common_keys &= set(item.keys())

    best_key = None
    best_size = -1
    for key in sorted(common_keys):
        lengths = {len(item[key]) for item in loaded}
        if len(lengths) != 1:
            continue
        length = next(iter(lengths))
        if length > best_size:
            best_size = length
            best_key = key

    if best_key is not None:
        matrix = np.vstack([item[best_key] for item in loaded])
        return matrix, best_key

    fallback = []
    for item in loaded:
        longest_key = max(item.keys(), key=lambda key: len(item[key]))
        fallback.append(item[longest_key])
    min_len = min(len(vector) for vector in fallback)
    if min_len < 100:
        raise ValueError("Unable to build a stable feature matrix from converted NPZ files.")
    matrix = np.vstack([vector[:min_len] for vector in fallback])
    return matrix, "fallback_longest_numeric_vector"


def normalize_matrix(matrix):
    matrix = np.clip(matrix, a_min=0.0, a_max=None)
    totals = matrix.sum(axis=1, keepdims=True)
    totals[totals <= 0.0] = 1.0
    matrix = matrix / totals * 1e6
    return np.log1p(matrix)


def build_folds(sample_count, fold_count, seed=42):
    rng = np.random.default_rng(seed)
    indices = rng.permutation(sample_count)
    return [fold for fold in np.array_split(indices, fold_count) if len(fold) > 0]


def pca_reconstruction_mse(matrix, n_components, fold_count=5):
    sample_count, feature_count = matrix.shape
    folds = build_folds(sample_count, min(fold_count, sample_count))
    fold_errors = []

    for test_indices in folds:
        train_mask = np.ones(sample_count, dtype=bool)
        train_mask[test_indices] = False
        train_indices = np.where(train_mask)[0]
        if len(train_indices) < 2:
            continue

        train_data = matrix[train_indices]
        test_data = matrix[test_indices]

        mean = train_data.mean(axis=0)
        std = train_data.std(axis=0)
        std[std == 0.0] = 1.0

        train_scaled = (train_data - mean) / std
        test_scaled = (test_data - mean) / std

        max_components = min(len(train_indices) - 1, feature_count)
        if max_components < 1 or n_components > max_components:
            return None

        _, _, vt = np.linalg.svd(train_scaled, full_matrices=False)
        basis = vt[:n_components]

        projected = test_scaled @ basis.T
        reconstructed = projected @ basis
        mse = np.mean((test_scaled - reconstructed) ** 2)
        fold_errors.append(mse)

    if not fold_errors:
        return None
    return float(np.mean(fold_errors))


def convert_all_bams(wisecondorx, bams, sample_ids, binsize, output_dir, threads, log_path):
    output_dir.mkdir(parents=True, exist_ok=True)
    npz_paths = []

    for sample_id, bam_path in zip(sample_ids, bams):
        npz_path = output_dir / f"{sample_id}.npz"
        npz_paths.append(npz_path)
        if npz_path.exists():
            continue

        command = [
            wisecondorx,
            "convert",
            str(bam_path),
            str(npz_path),
            "--binsize",
            str(binsize),
            "--cpus",
            str(threads),
        ]
        run_command(command, log_path)

    return npz_paths


def build_reference(wisecondorx, binsize, npz_paths, reference_output, threads, log_path):
    reference_output.parent.mkdir(parents=True, exist_ok=True)
    command = [
        wisecondorx,
        "newref",
        *[str(path) for path in npz_paths],
        str(reference_output),
        "--binsize",
        str(binsize),
        "--cpus",
        str(threads),
    ]
    run_command(command, log_path)


def write_summary(summary_output, rows):
    summary_output.parent.mkdir(parents=True, exist_ok=True)
    with open(summary_output, "w", encoding="utf-8") as handle:
        handle.write("binsize\tpca_components\tcv_reconstruction_mse\tsignal_key\tsamples\tfeatures\n")
        for row in rows:
            handle.write(
                f"{row['binsize']}\t{row['pca_components']}\t{row['score']:.10f}\t"
                f"{row['signal_key']}\t{row['samples']}\t{row['features']}\n"
            )


def write_best_yaml(best_output, best_row):
    best_output.parent.mkdir(parents=True, exist_ok=True)
    with open(best_output, "w", encoding="utf-8") as handle:
        handle.write(f"best_binsize: {best_row['binsize']}\n")
        handle.write(f"best_pca_components: {best_row['pca_components']}\n")
        handle.write(f"best_cv_reconstruction_mse: {best_row['score']:.10f}\n")
        handle.write(f"signal_key: {best_row['signal_key']}\n")


def main():
    parser = argparse.ArgumentParser(description="Tune WisecondorX bin size and PCA components by CV error.")
    parser.add_argument("--wisecondorx", required=True, help="Path to WisecondorX binary")
    parser.add_argument("--bams", nargs="+", required=True, help="Sorted BAM files")
    parser.add_argument("--sample-ids", nargs="+", required=True, help="Sample IDs, same order as BAMs")
    parser.add_argument("--bin-sizes", required=True, help="Comma-separated bin sizes")
    parser.add_argument("--pca-components", required=True, help="Comma-separated PCA component candidates")
    parser.add_argument("--threads", type=int, default=4, help="CPUs for convert/newref")
    parser.add_argument("--workdir", required=True, help="Tuning workspace directory")
    parser.add_argument("--summary-output", required=True, help="Output TSV for all combinations")
    parser.add_argument("--best-output", required=True, help="Output YAML for best parameters")
    parser.add_argument("--reference-output", required=True, help="Final reference output path")
    parser.add_argument("--log", required=True, help="Log file path")
    args = parser.parse_args()

    bams = [Path(path) for path in args.bams]
    sample_ids = args.sample_ids
    if len(bams) != len(sample_ids):
        raise ValueError("--bams and --sample-ids must have the same number of items.")
    if len(bams) < 3:
        raise ValueError("At least 3 samples are required for bin/PCA gradient analysis.")

    bin_sizes = sorted(set(parse_int_list(args.bin_sizes)))
    pca_candidates = sorted(set(parse_int_list(args.pca_components)))
    if min(pca_candidates) < 1:
        raise ValueError("PCA components must be >= 1.")

    workdir = Path(args.workdir)
    summary_output = Path(args.summary_output)
    best_output = Path(args.best_output)
    reference_output = Path(args.reference_output)
    log_path = Path(args.log)

    result_rows = []
    converted_by_bin = {}

    for binsize in bin_sizes:
        converted_dir = workdir / f"bin_{binsize}" / "converted"
        npz_paths = convert_all_bams(
            wisecondorx=args.wisecondorx,
            bams=bams,
            sample_ids=sample_ids,
            binsize=binsize,
            output_dir=converted_dir,
            threads=args.threads,
            log_path=log_path,
        )
        converted_by_bin[binsize] = npz_paths

        matrix, signal_key = build_matrix(npz_paths)
        matrix = normalize_matrix(matrix)
        samples, features = matrix.shape

        for n_components in pca_candidates:
            score = pca_reconstruction_mse(matrix, n_components)
            if score is None:
                continue
            result_rows.append(
                {
                    "binsize": binsize,
                    "pca_components": n_components,
                    "score": score,
                    "signal_key": signal_key,
                    "samples": samples,
                    "features": features,
                }
            )

    if not result_rows:
        raise ValueError("No valid bin/PCA combinations were scored. Check sample size and candidate settings.")

    result_rows.sort(key=lambda row: row["score"])
    best_row = result_rows[0]

    write_summary(summary_output, result_rows)
    write_best_yaml(best_output, best_row)

    build_reference(
        wisecondorx=args.wisecondorx,
        binsize=best_row["binsize"],
        npz_paths=converted_by_bin[best_row["binsize"]],
        reference_output=reference_output,
        threads=args.threads,
        log_path=log_path,
    )


if __name__ == "__main__":
    main()
