#!/biosoftware/miniconda/envs/snakemake_env/bin/python
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


SEX_CHROM_ALIASES = {
    "chrX": ("chrX", "X", "23", 23),
    "X": ("chrX", "X", "23", 23),
    "23": ("chrX", "X", "23", 23),
    "chrY": ("chrY", "Y", "24", 24),
    "Y": ("chrY", "Y", "24", 24),
    "24": ("chrY", "Y", "24", 24),
    "chrM": ("chrM", "M", "MT", "25", 25),
    "M": ("chrM", "M", "MT", "25", 25),
    "MT": ("chrM", "M", "MT", "25", 25),
}


def chrom_key_candidates(chrom):
    chrom = str(chrom)
    candidates = [chrom]
    if chrom.startswith("chr"):
        bare = chrom[3:]
        candidates.append(bare)
        if bare.isdigit():
            candidates.append(int(bare))
    elif chrom.isdigit():
        candidates.extend([int(chrom), f"chr{chrom}"])
    candidates.extend(SEX_CHROM_ALIASES.get(chrom, ()))
    unique = []
    for item in candidates:
        if item not in unique:
            unique.append(item)
    return unique


def load_wisecondorx_npz(path_value):
    path = Path(path_value)
    with np.load(path, allow_pickle=True, encoding="latin1") as data:
        if "sample" not in data.files:
            raise ValueError(f"Missing 'sample' key in {path}")
        sample = data["sample"].item()
        if not isinstance(sample, dict):
            raise ValueError(f"'sample' is not dict-like in {path}")
        binsize = int(np.asarray(data["binsize"]).item()) if "binsize" in data.files else 0
        payload = {"binsize": binsize, "sample": sample}
        if "quality" in data.files:
            payload["quality"] = data["quality"]
    return payload


def validate_wisecondorx_count_like_npz(sample_dict):
    for chrom, raw_values in sample_dict.items():
        values = np.asarray(raw_values)
        if values.size == 0:
            continue
        if not np.issubdtype(values.dtype, np.number):
            raise ValueError(f"Non-numeric count vector for chromosome {chrom!r}")
        numeric = values.astype(np.float64, copy=False)
        if not np.all(np.isfinite(numeric)):
            raise ValueError(f"Non-finite count vector values for chromosome {chrom!r}")
        if np.any(numeric < 0.0):
            raise ValueError(f"Negative count vector values for chromosome {chrom!r}")
        if np.all(numeric == 0.0) and str(chrom) not in {"chrM", "M", "MT", "25"}:
            raise ValueError(f"All-zero count vector for chromosome {chrom!r}")


def build_gc_rc_fit_mask(bins_df, include_mask_labels={"pass"}):
    frame = bins_df.copy()
    if "is_autosome" in frame.columns:
        is_autosome = pd.to_numeric(frame["is_autosome"], errors="coerce").fillna(0).astype(int).eq(1)
    else:
        is_autosome = frame["chrom"].astype(str).str.match(r"^(chr)?([1-9]|1[0-9]|2[0-2])$")
    mask_labels = frame.get("mask_label", pd.Series("pass", index=frame.index)).fillna("pass").astype(str)
    required_numeric = ["normalized_signal", "gc_fraction", "mappability_score"]
    valid = is_autosome & mask_labels.isin(set(include_mask_labels))
    for column in required_numeric:
        if column not in frame.columns:
            return pd.Series(False, index=frame.index)
        values = pd.to_numeric(frame[column], errors="coerce")
        valid &= np.isfinite(values.to_numpy(dtype=np.float64))
    return valid


def _select_mask_rows(combined_mask, binsize):
    frame = combined_mask.copy()
    if "mask_label" not in frame.columns:
        raise ValueError("combined mask table must contain mask_label")
    if "bin_size" in frame.columns and binsize:
        same_size = pd.to_numeric(frame["bin_size"], errors="coerce").fillna(0).astype(int).eq(int(binsize))
        if same_size.any():
            frame = frame[same_size].copy()
    return frame[frame["mask_label"].astype(str).eq("hard")].copy()


def mask_npz_sample(sample_dict, binsize, combined_mask_df):
    if binsize <= 0:
        raise ValueError("WisecondorX NPZ binsize must be positive for mask-only preprocessing")

    masked = {key: np.asarray(value).copy() for key, value in sample_dict.items()}
    original_shapes = {key: np.asarray(value).shape for key, value in sample_dict.items()}
    hard_rows = _select_mask_rows(combined_mask_df, binsize)
    soft_count = int(combined_mask_df.get("mask_label", pd.Series(dtype=str)).astype(str).eq("soft").sum())
    hard_masked = 0

    for row in hard_rows.itertuples(index=False):
        row_dict = row._asdict()
        start = int(row_dict["start"])
        end = int(row_dict["end"])
        left_bin = max(int(start // binsize), 0)
        right_bin = max(int(np.ceil(float(end) / float(binsize))), left_bin + 1)
        for key in chrom_key_candidates(row_dict["chrom"]):
            if key not in masked:
                continue
            vector = np.asarray(masked[key]).copy()
            stop = min(right_bin, vector.size)
            if left_bin < stop:
                vector[left_bin:stop] = 0
                hard_masked += int(stop - left_bin)
                masked[key] = vector
            break

    for key, shape in original_shapes.items():
        if np.asarray(masked[key]).shape != shape:
            raise ValueError(f"Shape changed for chromosome {key!r} during mask-only preprocessing")

    validate_wisecondorx_count_like_npz(masked)
    metadata = {
        "preprocess_version": "build_ref_v2",
        "preprocess_strategy": "mask_only",
        "mask_policy": "hard_mask_to_zero_soft_keep",
        "hard_masked_bin_count": int(hard_masked),
        "soft_masked_bin_count": int(soft_count),
        "chrM_excluded": True,
        "par_policy": "annotation_only",
        "xy_homology_policy": "annotation_only",
        "count_like_contract_validated": True,
    }
    return masked, metadata


def write_wisecondorx_npz(output_npz, binsize, sample_dict, quality=None, metadata=None):
    output = Path(output_npz)
    output.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "binsize": np.asarray(binsize),
        "sample": sample_dict,
    }
    if quality is not None:
        payload["quality"] = quality
    if metadata is not None:
        payload["pgta_preprocess"] = metadata
    np.savez_compressed(output, **payload)


def run_mask_npz(args):
    payload = load_wisecondorx_npz(args.input_npz)
    combined_mask = pd.read_csv(args.combined_mask, sep="\t")
    masked, metadata = mask_npz_sample(payload["sample"], payload["binsize"], combined_mask)
    write_wisecondorx_npz(
        args.output_npz,
        payload["binsize"],
        masked,
        quality=payload.get("quality"),
        metadata=metadata,
    )
    summary_path = Path(args.output_summary)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")


def parse_args():
    parser = argparse.ArgumentParser(description="Build Ref V2 WisecondorX-compatible NPZ preprocessing.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    mask_npz = subparsers.add_parser("mask_npz", help="Apply hard-mask-to-zero mask-only preprocessing to a WisecondorX NPZ.")
    mask_npz.add_argument("--input-npz", required=True)
    mask_npz.add_argument("--annotations", default="")
    mask_npz.add_argument("--combined-mask", required=True)
    mask_npz.add_argument("--output-npz", required=True)
    mask_npz.add_argument("--output-summary", required=True)
    mask_npz.set_defaults(func=run_mask_npz)
    return parser.parse_args()


def main():
    args = parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
