#!/biosoftware/miniconda/envs/snakemake_env/bin/python
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from pgta.predict.branch_b.common import STATE_TO_SVTYPE, write_json, write_table


def parse_args():
    parser = argparse.ArgumentParser(description="Assemble high-recall Branch A candidates from WisecondorX BED output.")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--input-bed", required=True)
    parser.add_argument("--output-candidates", required=True)
    parser.add_argument("--output-summary", required=True)
    parser.add_argument("--merge-gap-bp", type=int, default=0)
    parser.add_argument("--strong-z", type=float, default=10.0)
    parser.add_argument("--log", default="")
    return parser.parse_args()


def normalize_chrom(value):
    text = str(value).strip()
    if not text:
        return text
    return text if text.startswith("chr") else f"chr{text}"


def normalize_state(value):
    text = str(value).strip().lower()
    if text in {"gain", "dup", "duplication"}:
        return "gain"
    if text in {"loss", "del", "deletion"}:
        return "loss"
    return text


def _chrom_sort_key(chrom):
    text = str(chrom)
    raw = text[3:] if text.startswith("chr") else text
    if raw.isdigit():
        return (0, int(raw))
    if raw == "X":
        return (1, 23)
    if raw == "Y":
        return (1, 24)
    return (2, raw)


def load_wisecondorx_bed(path_value, sample_id):
    path = Path(path_value)
    if not path.exists():
        raise FileNotFoundError(f"Branch A BED not found: {path}")
    if path.stat().st_size == 0:
        return pd.DataFrame(columns=["sample_id", "chrom", "start", "end", "state", "a_ratio", "a_zscore"])
    df = pd.read_csv(path, sep="\t")
    rename_map = {"chr": "chrom", "type": "state", "ratio": "a_ratio", "zscore": "a_zscore"}
    df = df.rename(columns={key: value for key, value in rename_map.items() if key in df.columns})
    required = {"chrom", "start", "end", "state", "a_ratio", "a_zscore"}
    missing = sorted(required.difference(df.columns))
    if missing:
        raise ValueError(f"Branch A BED missing columns: {','.join(missing)}")
    out = df[list(required)].copy()
    out["sample_id"] = sample_id
    out["chrom"] = out["chrom"].map(normalize_chrom)
    out["state"] = out["state"].map(normalize_state)
    out["start"] = pd.to_numeric(out["start"], errors="coerce").astype("Int64")
    out["end"] = pd.to_numeric(out["end"], errors="coerce").astype("Int64")
    out["a_ratio"] = pd.to_numeric(out["a_ratio"], errors="coerce")
    out["a_zscore"] = pd.to_numeric(out["a_zscore"], errors="coerce")
    out = out.dropna(subset=["start", "end", "a_zscore", "chrom", "state"]).copy()
    out["start"] = out["start"].astype(int)
    out["end"] = out["end"].astype(int)
    out = out[out["end"] > out["start"]].copy()
    out["a_abs_zscore"] = out["a_zscore"].abs()
    out["_chrom_sort"] = out["chrom"].map(_chrom_sort_key)
    return out.sort_values(["_chrom_sort", "state", "start", "end"]).drop(columns=["_chrom_sort"]).reset_index(drop=True)


def merge_a_branch_candidates(raw_df, merge_gap_bp=0):
    if raw_df.empty:
        return pd.DataFrame(
            columns=[
                "sample_id",
                "chrom",
                "start",
                "end",
                "state",
                "a_ratio",
                "a_zscore",
                "a_abs_zscore",
                "a_source_event_count",
            ]
        )
    merged = []
    gap = max(int(merge_gap_bp), 0)
    for (_sample_id, chrom, state), group in raw_df.groupby(["sample_id", "chrom", "state"], sort=False):
        current = None
        for row in group.sort_values(["start", "end"]).itertuples(index=False):
            span = max(int(row.end) - int(row.start), 1)
            payload = {
                "sample_id": str(row.sample_id),
                "chrom": str(row.chrom),
                "start": int(row.start),
                "end": int(row.end),
                "state": str(row.state),
                "_weighted_ratio_sum": float(row.a_ratio) * span if np.isfinite(float(row.a_ratio)) else 0.0,
                "_ratio_weight_sum": span if np.isfinite(float(row.a_ratio)) else 0,
                "a_zscore": float(row.a_zscore),
                "a_abs_zscore": abs(float(row.a_zscore)),
                "a_source_event_count": 1,
            }
            if current is None:
                current = payload
                continue
            if payload["start"] <= current["end"] + gap + 1:
                current["end"] = max(current["end"], payload["end"])
                current["_weighted_ratio_sum"] += payload["_weighted_ratio_sum"]
                current["_ratio_weight_sum"] += payload["_ratio_weight_sum"]
                current["a_source_event_count"] += 1
                if payload["a_abs_zscore"] > current["a_abs_zscore"]:
                    current["a_zscore"] = payload["a_zscore"]
                    current["a_abs_zscore"] = payload["a_abs_zscore"]
            else:
                merged.append(current)
                current = payload
        if current is not None:
            merged.append(current)
    out = pd.DataFrame(merged)
    if out.empty:
        return out
    denom = out["_ratio_weight_sum"].replace(0, np.nan)
    out["a_ratio"] = out["_weighted_ratio_sum"] / denom
    out["a_ratio"] = out["a_ratio"].fillna(0.0)
    out = out.drop(columns=["_weighted_ratio_sum", "_ratio_weight_sum"])
    out["_chrom_sort"] = out["chrom"].map(_chrom_sort_key)
    return out.sort_values(["a_abs_zscore", "_chrom_sort", "start"], ascending=[False, True, True]).drop(columns=["_chrom_sort"]).reset_index(drop=True)


def finalize_a_candidates(merged_df, strong_z=10.0):
    out = merged_df.copy()
    if out.empty:
        for column in ["candidate_id", "event_id", "svtype", "a_rank", "a_support_level", "caller", "caller_stage", "branch"]:
            out[column] = []
        return out
    out["a_rank"] = np.arange(1, len(out) + 1, dtype=int)
    out["candidate_id"] = out.apply(lambda row: f"{row.sample_id}.A{int(row.a_rank):04d}", axis=1)
    out["event_id"] = out["candidate_id"]
    out["svtype"] = out["state"].map(STATE_TO_SVTYPE).fillna("")
    out["a_support_level"] = np.where(out["a_abs_zscore"].astype(float) >= float(strong_z), "strong", "sensitive")
    out["caller"] = "wisecondorx"
    out["caller_stage"] = "branch_a_merged"
    out["branch"] = "A"
    ordered = [
        "candidate_id",
        "event_id",
        "sample_id",
        "branch",
        "caller",
        "caller_stage",
        "chrom",
        "start",
        "end",
        "state",
        "svtype",
        "a_ratio",
        "a_zscore",
        "a_abs_zscore",
        "a_rank",
        "a_support_level",
        "a_source_event_count",
    ]
    return out[ordered]


def assemble_a_branch_candidates(input_bed, sample_id, merge_gap_bp=0, strong_z=10.0):
    raw = load_wisecondorx_bed(input_bed, sample_id)
    merged = merge_a_branch_candidates(raw, merge_gap_bp=merge_gap_bp)
    return finalize_a_candidates(merged, strong_z=strong_z), raw


def main():
    args = parse_args()
    candidates, raw = assemble_a_branch_candidates(
        args.input_bed,
        args.sample_id,
        merge_gap_bp=args.merge_gap_bp,
        strong_z=args.strong_z,
    )
    write_table(args.output_candidates, candidates)
    write_json(
        args.output_summary,
        {
            "sample_id": args.sample_id,
            "raw_event_count": int(len(raw)),
            "candidate_event_count": int(len(candidates)),
            "strong_event_count": int(candidates["a_support_level"].eq("strong").sum()) if not candidates.empty else 0,
            "merge_gap_bp": int(args.merge_gap_bp),
            "strong_z": float(args.strong_z),
        },
    )


if __name__ == "__main__":
    main()
