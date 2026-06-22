#!/biosoftware/miniconda/envs/snakemake_env/bin/python
from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd

from pgta.predict.branch_b.common import load_sample_bins, read_table, robust_z, write_json, write_table


REF_STABILITY_VERSION = "branch_b_v2_ref_stability_v1"
DEFAULT_MODERATE_MAD_Z = 2.0
DEFAULT_HIGH_MAD_Z = 4.0


def parse_args():
    parser = argparse.ArgumentParser(description="Compute reference cohort per-bin MAD stability context.")
    parser.add_argument("--input-npz", action="append", default=[])
    parser.add_argument("--sample-id", action="append", default=[])
    parser.add_argument("--sample-sex", action="append", default=[])
    parser.add_argument("--input-bins", default="")
    parser.add_argument("--input-events", default="")
    parser.add_argument("--output-bins", default="")
    parser.add_argument("--output-events", default="")
    parser.add_argument("--output-summary", required=True)
    parser.add_argument("--moderate-mad-z", type=float, default=DEFAULT_MODERATE_MAD_Z)
    parser.add_argument("--high-mad-z", type=float, default=DEFAULT_HIGH_MAD_Z)
    return parser.parse_args()


def _safe_float(value, default=math.nan):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return float(default)
    return number if math.isfinite(number) else float(default)


def _label_ref_stability(ref_mad_z, moderate_mad_z, high_mad_z):
    if math.isfinite(ref_mad_z) and ref_mad_z >= high_mad_z:
        return "high_ref_mad"
    if math.isfinite(ref_mad_z) and ref_mad_z >= moderate_mad_z:
        return "moderate_ref_mad"
    return "stable"


def _event_context(median_z, high_fraction, moderate_mad_z, high_mad_z):
    if math.isfinite(high_fraction) and high_fraction >= 0.25:
        return "REF_STABILITY_HIGH_MAD_REVIEW"
    if math.isfinite(median_z) and median_z >= high_mad_z:
        return "REF_STABILITY_HIGH_MAD_REVIEW"
    if math.isfinite(high_fraction) and high_fraction > 0.0:
        return "REF_STABILITY_MODERATE_REVIEW"
    if math.isfinite(median_z) and median_z >= moderate_mad_z:
        return "REF_STABILITY_MODERATE_REVIEW"
    return "REF_STABILITY_STABLE"


def _summary_payload(sample_ids, binsize, bins, high_mad_z, sample_sexes=None):
    high_fraction = 0.0
    if not bins.empty:
        high_fraction = float((pd.to_numeric(bins["ref_mad_z"], errors="coerce") >= high_mad_z).mean())
    return {
        "version": REF_STABILITY_VERSION,
        "ref_sample_count": int(len(sample_ids)),
        "ref_sample_ids": [str(item) for item in sample_ids],
        "ref_sample_sexes": [str(item) for item in sample_sexes or []],
        "binsize": int(binsize),
        "bin_count": int(len(bins)),
        "high_ref_mad_bin_fraction": high_fraction,
        "final_report_impact": "development_review_only",
    }


def compute_ref_bin_stability(
    npz_paths,
    sample_ids=None,
    sample_sexes=None,
    moderate_mad_z=DEFAULT_MODERATE_MAD_Z,
    high_mad_z=DEFAULT_HIGH_MAD_Z,
):
    paths = [Path(path) for path in npz_paths]
    if not paths:
        raise ValueError("at least one reference NPZ is required")
    if sample_ids and len(sample_ids) != len(paths):
        raise ValueError("--sample-id count must match --input-npz count when supplied")
    if sample_sexes and len(sample_sexes) != len(paths):
        raise ValueError("--sample-sex count must match --input-npz count when supplied")
    resolved_sample_ids = [str(sample_ids[index]) if sample_ids else path.stem for index, path in enumerate(paths)]
    resolved_sample_sexes = [str(sample_sexes[index]).strip().upper() if sample_sexes else "" for index, _path in enumerate(paths)]

    frames = []
    binsize_seen = None
    for path, sample_id, sample_sex in zip(paths, resolved_sample_ids, resolved_sample_sexes):
        bins, binsize, _quality = load_sample_bins(path)
        if binsize_seen is None:
            binsize_seen = binsize
        elif binsize != binsize_seen:
            raise ValueError(f"inconsistent binsize in reference NPZ inputs: {binsize_seen} vs {binsize}")
        sample_bins = bins[["chrom", "bin_index", "start", "end", "normalized_signal"]].copy()
        sample_bins["ref_sample_id"] = sample_id
        sample_bins["ref_sample_sex"] = sample_sex if sample_sex in {"XX", "XY"} else ""
        frames.append(sample_bins)

    combined = pd.concat(frames, ignore_index=True)
    def _mad(series):
        values = pd.to_numeric(series, errors="coerce").dropna().to_numpy(dtype=float)
        if values.size == 0:
            return math.nan
        median = float(np.median(values))
        return float(np.median(np.abs(values - median)))

    def _group_stability(input_frame, ref_group):
        grouped = (
            input_frame.groupby(["chrom", "bin_index", "start", "end"], as_index=False)
            .agg(
                ref_sample_count=("normalized_signal", "count"),
                ref_median=("normalized_signal", "median"),
                ref_std=("normalized_signal", "std"),
                ref_mad=("normalized_signal", _mad),
            )
            .sort_values(["chrom", "bin_index"])
            .reset_index(drop=True)
        )
        grouped["ref_group"] = ref_group
        return grouped

    grouped_frames = [_group_stability(combined, "mixed")]
    if any(sex in {"XX", "XY"} for sex in resolved_sample_sexes):
        for sex in ("XX", "XY"):
            sex_frame = combined[combined["ref_sample_sex"].astype(str).eq(sex)].copy()
            if not sex_frame.empty:
                grouped_frames.append(_group_stability(sex_frame, sex))

    grouped = pd.concat(grouped_frames, ignore_index=True, sort=False)
    grouped = (
        grouped
        .sort_values(["chrom", "bin_index"])
        .reset_index(drop=True)
    )
    grouped["ref_std"] = pd.to_numeric(grouped["ref_std"], errors="coerce").fillna(0.0)
    grouped["ref_cv"] = grouped["ref_std"] / grouped["ref_median"].abs().clip(lower=1e-6)
    grouped["ref_mad_z"] = robust_z(grouped["ref_mad"].fillna(0.0).to_numpy(dtype=float))
    grouped["ref_stability_label"] = [
        _label_ref_stability(value, moderate_mad_z, high_mad_z)
        for value in grouped["ref_mad_z"].to_numpy(dtype=float)
    ]
    summary = _summary_payload(resolved_sample_ids, int(binsize_seen or 0), grouped, high_mad_z, resolved_sample_sexes)
    if "ref_group" in grouped.columns:
        summary["ref_group_counts"] = {
            str(key): int(value)
            for key, value in grouped["ref_group"].fillna("mixed").astype(str).value_counts().sort_index().to_dict().items()
        }
    return grouped, summary


def summarize_event_ref_stability(
    events,
    ref_bins,
    moderate_mad_z=DEFAULT_MODERATE_MAD_Z,
    high_mad_z=DEFAULT_HIGH_MAD_Z,
):
    frame = events.copy()
    if frame.empty:
        for column in [
            "event_ref_mad_median",
            "event_ref_mad_p90",
            "high_ref_mad_bin_fraction",
            "ref_stability_context",
        ]:
            frame[column] = []
        return frame

    ref_frame = ref_bins.copy()
    ref_frame["chrom"] = ref_frame["chrom"].astype(str)
    ref_frame["start"] = pd.to_numeric(ref_frame["start"], errors="coerce")
    ref_frame["end"] = pd.to_numeric(ref_frame["end"], errors="coerce")
    ref_frame["ref_mad"] = pd.to_numeric(ref_frame["ref_mad"], errors="coerce")
    ref_frame["ref_mad_z"] = pd.to_numeric(ref_frame["ref_mad_z"], errors="coerce")

    payloads = []
    for _, event in frame.iterrows():
        chrom = str(event.get("chrom", ""))
        if chrom and not chrom.lower().startswith("chr"):
            chrom = f"chr{chrom}"
        start = int(_safe_float(event.get("start", 0), default=0))
        end = int(_safe_float(event.get("end", 0), default=0))
        bins = ref_frame.loc[
            ref_frame["chrom"].eq(chrom)
            & ref_frame["start"].notna()
            & ref_frame["end"].notna()
            & (ref_frame["end"] > start)
            & (ref_frame["start"] < end)
        ]
        if bins.empty:
            payloads.append(
                {
                    "event_ref_mad_median": math.nan,
                    "event_ref_mad_p90": math.nan,
                    "high_ref_mad_bin_fraction": math.nan,
                    "ref_stability_context": "REF_STABILITY_NO_BIN_CONTEXT",
                }
            )
            continue
        mad = pd.to_numeric(bins["ref_mad"], errors="coerce").dropna()
        mad_z = pd.to_numeric(bins["ref_mad_z"], errors="coerce").dropna()
        median_mad = float(mad.median()) if not mad.empty else math.nan
        p90_mad = float(mad.quantile(0.90)) if not mad.empty else math.nan
        high_fraction = float((mad_z >= high_mad_z).mean()) if not mad_z.empty else math.nan
        median_z = float(mad_z.median()) if not mad_z.empty else math.nan
        payloads.append(
            {
                "event_ref_mad_median": round(median_mad, 6) if math.isfinite(median_mad) else math.nan,
                "event_ref_mad_p90": round(p90_mad, 6) if math.isfinite(p90_mad) else math.nan,
                "high_ref_mad_bin_fraction": round(high_fraction, 6) if math.isfinite(high_fraction) else math.nan,
                "ref_stability_context": _event_context(median_z, high_fraction, moderate_mad_z, high_mad_z),
            }
        )

    return pd.concat([frame.reset_index(drop=True), pd.DataFrame(payloads)], axis=1)


def main():
    args = parse_args()
    if args.input_bins:
        bins = read_table(args.input_bins, empty_ok=False)
        summary = _summary_payload(args.sample_id, 0, bins, args.high_mad_z, args.sample_sex)
    else:
        bins, summary = compute_ref_bin_stability(
            args.input_npz,
            sample_ids=args.sample_id or None,
            sample_sexes=args.sample_sex or None,
            moderate_mad_z=args.moderate_mad_z,
            high_mad_z=args.high_mad_z,
        )
    if args.output_bins:
        write_table(args.output_bins, bins)
    if args.input_events and args.output_events:
        events = read_table(args.input_events, empty_ok=True)
        event_summary = summarize_event_ref_stability(
            events,
            bins,
            moderate_mad_z=args.moderate_mad_z,
            high_mad_z=args.high_mad_z,
        )
        write_table(args.output_events, event_summary)
        summary["event_count"] = int(len(event_summary))
        if "ref_stability_context" in event_summary.columns:
            summary["event_ref_stability_context_counts"] = {
                str(key): int(value)
                for key, value in event_summary["ref_stability_context"]
                .fillna("UNKNOWN")
                .astype(str)
                .value_counts()
                .sort_index()
                .to_dict()
                .items()
            }
    write_json(args.output_summary, summary)


if __name__ == "__main__":
    main()
