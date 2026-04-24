#!/biosoftware/miniconda/envs/snakemake_env/bin/python
import argparse
import json
import math
from pathlib import Path

def add_annotations_parser(subparsers):
    parser = subparsers.add_parser("annotations", help="Build atomic/analysis/QC bins and bin annotations.")
    parser.add_argument("--reference-fasta", required=True)
    parser.add_argument("--chromosomes", nargs="+", required=True)
    parser.add_argument("--atomic-bin-size", type=int, required=True)
    parser.add_argument("--analysis-bin-size", type=int, required=True)
    parser.add_argument("--qc-bin-size", type=int, required=True)
    parser.add_argument("--atomic-bins-output", required=True)
    parser.add_argument("--analysis-bins-output", required=True)
    parser.add_argument("--qc-bins-output", required=True)
    parser.add_argument("--atomic-annotations-output", required=True)
    parser.add_argument("--analysis-annotations-output", required=True)
    parser.add_argument("--qc-annotations-output", required=True)
    parser.add_argument("--summary-json-output", required=True)
    parser.add_argument("--par-region", action="append", default=[])
    parser.add_argument("--xtr-region", action="append", default=[])
    parser.add_argument("--segmental-duplication-bed", default="")
    parser.add_argument("--low-mappability-bed", default="")
    parser.add_argument("--gap-centromere-telomere-bed", default="")
    parser.add_argument("--repeat-rich-bed", default="")
    parser.add_argument("--blacklist-bed", default="")
    parser.add_argument("--sex-homology-bed", default="")
    parser.add_argument("--ambiguous-alignment-bed", default="")
    parser.add_argument("--log", default="")
    parser.set_defaults(command="annotations")


def add_mask_parser(subparsers):
    parser = subparsers.add_parser("mask", help="Build hard/soft/dynamic masks from bin annotations and profiles.")
    parser.add_argument("--annotation-tsvs", nargs="+", required=True)
    parser.add_argument("--profile-tsvs", nargs="*", default=[])
    parser.add_argument("--hard-mask-output", required=True)
    parser.add_argument("--soft-mask-output", required=True)
    parser.add_argument("--dynamic-mask-output", required=True)
    parser.add_argument("--combined-mask-output", required=True)
    parser.add_argument("--hard-mask-json-output", required=True)
    parser.add_argument("--soft-mask-json-output", required=True)
    parser.add_argument("--dynamic-mask-json-output", required=True)
    parser.add_argument("--combined-mask-json-output", required=True)
    parser.add_argument("--summary-json-output", required=True)
    parser.add_argument("--hard-n-fraction", type=float, default=0.2)
    parser.add_argument("--soft-gc-low", type=float, default=0.2)
    parser.add_argument("--soft-gc-high", type=float, default=0.8)
    parser.add_argument("--dynamic-z-frac-threshold", type=float, default=0.05)
    parser.add_argument("--dynamic-median-abs-z-threshold", type=float, default=1.5)
    parser.add_argument("--log", default="")
    parser.set_defaults(command="mask")


def parse_args():
    parser = argparse.ArgumentParser(description="Build reference bin annotations or masks.")
    subparsers = parser.add_subparsers(dest="command", required=True)
    add_annotations_parser(subparsers)
    add_mask_parser(subparsers)
    return parser.parse_args()


def interval_overlap(start, end, left, right):
    return max(0, min(int(end), int(right)) - max(int(start), int(left)))


def normalize_chromosomes(fasta, requested):
    available = set(fasta.references)
    chroms = [chrom for chrom in requested if chrom in available]
    if not chroms:
        raise ValueError("None of the requested chromosomes exist in the reference FASTA")
    return chroms


def build_bins(chroms, fasta, bin_size, level):
    import pandas as pd

    rows = []
    for chrom in chroms:
        chrom_length = fasta.get_reference_length(chrom)
        n_bins = int(math.ceil(chrom_length / float(bin_size)))
        for idx in range(n_bins):
            start = idx * bin_size
            end = min((idx + 1) * bin_size, chrom_length)
            rows.append(
                {
                    "bin_level": level,
                    "bin_id": f"{level}:{chrom}:{start}-{end}",
                    "chrom": chrom,
                    "start": int(start),
                    "end": int(end),
                    "bin_size": int(bin_size),
                    "effective_size": int(end - start),
                }
            )
    return pd.DataFrame(rows)


def relabel_bins_for_level(df, level):
    relabeled = df.copy()
    relabeled["bin_level"] = level
    relabeled["bin_id"] = (
        level
        + ":"
        + relabeled["chrom"].astype(str)
        + ":"
        + relabeled["start"].astype(str)
        + "-"
        + relabeled["end"].astype(str)
    )
    return relabeled


def parse_region_specs(region_specs):
    parsed = {}
    for spec in region_specs:
        chrom_part, sep, range_part = str(spec).partition(":")
        start_part, sep2, end_part = range_part.partition("-")
        if not sep or not sep2:
            raise ValueError(f"Invalid region specification: {spec!r}")
        chrom = chrom_part.strip()
        parsed.setdefault(chrom, []).append((int(start_part), int(end_part)))
    for chrom in parsed:
        parsed[chrom] = sorted(parsed[chrom])
    return parsed


def load_bed_intervals(path_value):
    if not path_value:
        return {}
    path = Path(path_value)
    if not path.exists():
        raise FileNotFoundError(f"Annotation BED does not exist: {path}")
    intervals = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0].strip()
            intervals.setdefault(chrom, []).append((int(parts[1]), int(parts[2])))
    for chrom in intervals:
        intervals[chrom] = sorted(intervals[chrom])
    return intervals


def overlap_fractions_for_bins(bins_df, interval_map):
    fractions = [0.0] * len(bins_df)
    if bins_df.empty or not interval_map:
        return fractions
    bins_reset = bins_df.reset_index(drop=True)
    for chrom, chrom_bins in bins_reset.groupby("chrom", sort=False):
        intervals = interval_map.get(chrom, [])
        if not intervals:
            continue
        interval_idx = 0
        n_intervals = len(intervals)
        for row in chrom_bins.itertuples():
            start = int(row.start)
            end = int(row.end)
            event_length = max(end - start, 1)
            while interval_idx < n_intervals and intervals[interval_idx][1] <= start:
                interval_idx += 1
            overlap_bp = 0
            scan_idx = interval_idx
            while scan_idx < n_intervals and intervals[scan_idx][0] < end:
                left, right = intervals[scan_idx]
                overlap_bp += interval_overlap(start, end, int(left), int(right))
                scan_idx += 1
            fractions[int(row.Index)] = min(float(overlap_bp) / float(event_length), 1.0)
    return fractions


def annotate_bins(fasta, bins_df, region_maps):
    import pandas as pd

    overlap_maps = {
        "par": overlap_fractions_for_bins(bins_df, region_maps["par"]),
        "xtr": overlap_fractions_for_bins(bins_df, region_maps["xtr"]),
        "sex_homology": overlap_fractions_for_bins(bins_df, region_maps["sex_homology"]),
        "segmental_duplication": overlap_fractions_for_bins(bins_df, region_maps["segmental_duplication"]),
        "low_mappability": overlap_fractions_for_bins(bins_df, region_maps["low_mappability"]),
        "gap_centromere_telomere": overlap_fractions_for_bins(bins_df, region_maps["gap_centromere_telomere"]),
        "repeat_rich": overlap_fractions_for_bins(bins_df, region_maps["repeat_rich"]),
        "blacklist": overlap_fractions_for_bins(bins_df, region_maps["blacklist"]),
        "ambiguous_alignment": overlap_fractions_for_bins(bins_df, region_maps["ambiguous_alignment"]),
    }
    annotations = []
    for idx, row in enumerate(bins_df.itertuples(index=False)):
        seq = fasta.fetch(row.chrom, int(row.start), int(row.end)).upper()
        length = len(seq)
        a_count = seq.count("A")
        c_count = seq.count("C")
        g_count = seq.count("G")
        t_count = seq.count("T")
        n_count = seq.count("N")
        gc_count = g_count + c_count
        atgc_count = a_count + c_count + g_count + t_count
        par_overlap_fraction = overlap_maps["par"][idx]
        xtr_overlap_fraction = overlap_maps["xtr"][idx]
        sex_homology_overlap_fraction = overlap_maps["sex_homology"][idx]
        segmental_duplication_overlap_fraction = overlap_maps["segmental_duplication"][idx]
        low_mappability_overlap_fraction = overlap_maps["low_mappability"][idx]
        gap_centromere_telomere_overlap_fraction = overlap_maps["gap_centromere_telomere"][idx]
        repeat_rich_overlap_fraction = overlap_maps["repeat_rich"][idx]
        blacklist_overlap_fraction = overlap_maps["blacklist"][idx]
        ambiguous_alignment_overlap_fraction = overlap_maps["ambiguous_alignment"][idx]
        annotations.append(
            {
                "bin_level": row.bin_level,
                "bin_id": row.bin_id,
                "chrom": row.chrom,
                "start": int(row.start),
                "end": int(row.end),
                "bin_size": int(row.bin_size),
                "effective_size": int(row.effective_size),
                "gc_fraction": (float(gc_count) / float(atgc_count)) if atgc_count else float("nan"),
                "n_fraction": (float(n_count) / float(length)) if length else float("nan"),
                "atgc_fraction": (float(atgc_count) / float(length)) if length else float("nan"),
                "mappability_score": (float(atgc_count) / float(length)) if length else float("nan"),
                "is_autosome": row.chrom.startswith("chr") and row.chrom[3:].isdigit() or row.chrom.isdigit(),
                "par_overlap_fraction": par_overlap_fraction,
                "xtr_overlap_fraction": xtr_overlap_fraction,
                "sex_homology_overlap_fraction": max(par_overlap_fraction, xtr_overlap_fraction, sex_homology_overlap_fraction),
                "segmental_duplication_overlap_fraction": segmental_duplication_overlap_fraction,
                "low_mappability_overlap_fraction": low_mappability_overlap_fraction,
                "gap_centromere_telomere_overlap_fraction": gap_centromere_telomere_overlap_fraction,
                "repeat_rich_overlap_fraction": repeat_rich_overlap_fraction,
                "blacklist_overlap_fraction": blacklist_overlap_fraction,
                "ambiguous_alignment_overlap_fraction": ambiguous_alignment_overlap_fraction,
                "is_PAR": int(par_overlap_fraction > 0.0),
                "is_XTR": int(xtr_overlap_fraction > 0.0),
                "is_sex_homology": int(max(par_overlap_fraction, xtr_overlap_fraction, sex_homology_overlap_fraction) > 0.0),
                "is_segmental_duplication": int(segmental_duplication_overlap_fraction > 0.0),
                "is_low_mappability": int(low_mappability_overlap_fraction > 0.0),
                "is_gap_centromere_telomere": int(gap_centromere_telomere_overlap_fraction > 0.0),
                "is_repeat_rich": int(repeat_rich_overlap_fraction > 0.0),
                "is_blacklist_region": int(blacklist_overlap_fraction > 0.0),
                "is_ambiguous_alignment_region": int(ambiguous_alignment_overlap_fraction > 0.0),
            }
        )
    return pd.DataFrame(annotations)


def write_tsv(path_value, df):
    path = Path(path_value)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def run_annotations(args):
    import pysam

    with pysam.FastaFile(args.reference_fasta) as fasta:
        chroms = normalize_chromosomes(fasta, args.chromosomes)
        region_maps = {
            "par": parse_region_specs(args.par_region),
            "xtr": parse_region_specs(args.xtr_region),
            "segmental_duplication": load_bed_intervals(args.segmental_duplication_bed),
            "low_mappability": load_bed_intervals(args.low_mappability_bed),
            "gap_centromere_telomere": load_bed_intervals(args.gap_centromere_telomere_bed),
            "repeat_rich": load_bed_intervals(args.repeat_rich_bed),
            "blacklist": load_bed_intervals(args.blacklist_bed),
            "sex_homology": load_bed_intervals(args.sex_homology_bed),
            "ambiguous_alignment": load_bed_intervals(args.ambiguous_alignment_bed),
        }
        levels = [
            ("atomic", args.atomic_bin_size, args.atomic_bins_output, args.atomic_annotations_output),
            ("analysis", args.analysis_bin_size, args.analysis_bins_output, args.analysis_annotations_output),
            ("qc", args.qc_bin_size, args.qc_bins_output, args.qc_annotations_output),
        ]
        summary = {
            "reference_fasta": args.reference_fasta,
            "chromosomes": chroms,
            "levels": {},
        }
        cache_by_bin_size = {}
        for level, bin_size, bins_output, annotations_output in levels:
            if bin_size not in cache_by_bin_size:
                cached_bins_df = build_bins(chroms, fasta, bin_size, level)
                cached_annotations_df = annotate_bins(fasta, cached_bins_df, region_maps)
                cache_by_bin_size[bin_size] = (cached_bins_df, cached_annotations_df)
            base_bins_df, base_annotations_df = cache_by_bin_size[bin_size]
            bins_df = relabel_bins_for_level(base_bins_df, level)
            annotations_df = relabel_bins_for_level(base_annotations_df, level)
            write_tsv(bins_output, bins_df)
            write_tsv(annotations_output, annotations_df)
            summary["levels"][level] = {
                "bin_size": int(bin_size),
                "n_bins": int(len(bins_df)),
                "gc_fraction_median": float(annotations_df["gc_fraction"].median(skipna=True)),
                "n_fraction_median": float(annotations_df["n_fraction"].median(skipna=True)),
                "par_bin_count": int(annotations_df["is_PAR"].sum()),
                "xtr_bin_count": int(annotations_df["is_XTR"].sum()),
                "sex_homology_bin_count": int(annotations_df["is_sex_homology"].sum()),
                "segmental_duplication_bin_count": int(annotations_df["is_segmental_duplication"].sum()),
                "low_mappability_bin_count": int(annotations_df["is_low_mappability"].sum()),
                "blacklist_bin_count": int(annotations_df["is_blacklist_region"].sum()),
            }

    summary_path = Path(args.summary_json_output)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")


def read_annotations(paths):
    import pandas as pd

    frames = [pd.read_csv(path, sep="\t") for path in paths]
    return pd.concat(frames, ignore_index=True)


def build_dynamic_metrics(profile_tsvs):
    import pandas as pd

    if not profile_tsvs:
        return pd.DataFrame(columns=["chrom", "start", "end", "dynamic_z_frac", "dynamic_median_abs_z", "sample_count"])
    frames = []
    for path_value in profile_tsvs:
        profile_df = pd.read_csv(path_value, sep="\t")
        required = {"chrom", "start", "end", "z_score"}
        if not required.issubset(profile_df.columns):
            continue
        profile_df = profile_df[["chrom", "start", "end", "z_score"]].copy()
        profile_df["abs_z"] = profile_df["z_score"].abs()
        profile_df["gt3"] = profile_df["abs_z"] > 3.0
        frames.append(profile_df)
    if not frames:
        return pd.DataFrame(columns=["chrom", "start", "end", "dynamic_z_frac", "dynamic_median_abs_z", "sample_count"])
    merged = pd.concat(frames, ignore_index=True)
    grouped = (
        merged.groupby(["chrom", "start", "end"], as_index=False)
        .agg(
            dynamic_z_frac=("gt3", "mean"),
            dynamic_median_abs_z=("abs_z", "median"),
            sample_count=("abs_z", "size"),
        )
    )
    return grouped


def write_tsv_json(tsv_path, json_path, df):
    tsv = Path(tsv_path)
    tsv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(tsv, sep="\t", index=False)
    Path(json_path).write_text(df.to_json(orient="records", indent=2), encoding="utf-8")


def run_mask(args):
    annotations = read_annotations(args.annotation_tsvs)
    dynamic_metrics = build_dynamic_metrics(args.profile_tsvs)
    merged = annotations.merge(dynamic_metrics, on=["chrom", "start", "end"], how="left")

    hard_mask = merged[merged["n_fraction"].fillna(0.0) >= args.hard_n_fraction].copy()
    hard_mask["mask_type"] = "hard"
    hard_mask["mask_reason"] = hard_mask["n_fraction"].map(lambda value: f"n_fraction>={args.hard_n_fraction}: {value:.4f}")

    soft_condition = (
        merged["gc_fraction"].fillna(0.5).lt(args.soft_gc_low)
        | merged["gc_fraction"].fillna(0.5).gt(args.soft_gc_high)
    )
    soft_mask = merged[soft_condition].copy()
    soft_mask["mask_type"] = "soft"
    soft_mask["mask_reason"] = soft_mask["gc_fraction"].map(
        lambda value: f"gc_fraction_outside_[{args.soft_gc_low},{args.soft_gc_high}]: {value:.4f}"
    )

    dynamic_condition = (
        merged["dynamic_z_frac"].fillna(0.0).ge(args.dynamic_z_frac_threshold)
        | merged["dynamic_median_abs_z"].fillna(0.0).ge(args.dynamic_median_abs_z_threshold)
    )
    dynamic_mask = merged[dynamic_condition].copy()
    dynamic_mask["mask_type"] = "dynamic"
    dynamic_mask["mask_reason"] = dynamic_mask.apply(
        lambda row: (
            f"dynamic_z_frac={row['dynamic_z_frac']:.4f};dynamic_median_abs_z={row['dynamic_median_abs_z']:.4f}"
        ),
        axis=1,
    )

    combined = merged.copy()
    combined["hard_masked"] = combined["bin_id"].isin(set(hard_mask["bin_id"]))
    combined["soft_masked"] = combined["bin_id"].isin(set(soft_mask["bin_id"]))
    combined["dynamic_masked"] = combined["bin_id"].isin(set(dynamic_mask["bin_id"]))
    combined["mask_label"] = "pass"
    combined.loc[combined["dynamic_masked"], "mask_label"] = "dynamic"
    combined.loc[combined["soft_masked"], "mask_label"] = "soft"
    combined.loc[combined["hard_masked"], "mask_label"] = "hard"
    combined["mask_reason"] = ""
    reason_map = {}
    for frame in (dynamic_mask, soft_mask, hard_mask):
        for row in frame[["bin_id", "mask_reason"]].itertuples(index=False):
            reason_map[row.bin_id] = row.mask_reason
    combined["mask_reason"] = combined["bin_id"].map(reason_map).fillna("")

    write_tsv_json(args.hard_mask_output, args.hard_mask_json_output, hard_mask)
    write_tsv_json(args.soft_mask_output, args.soft_mask_json_output, soft_mask)
    write_tsv_json(args.dynamic_mask_output, args.dynamic_mask_json_output, dynamic_mask)
    write_tsv_json(args.combined_mask_output, args.combined_mask_json_output, combined)

    summary = {
        "n_bins_total": int(len(combined)),
        "n_hard_mask": int(len(hard_mask)),
        "n_soft_mask": int(len(soft_mask)),
        "n_dynamic_mask": int(len(dynamic_mask)),
        "mask_label_counts": combined["mask_label"].value_counts(dropna=False).to_dict(),
    }
    Path(args.summary_json_output).write_text(json.dumps(summary, indent=2), encoding="utf-8")


def main():
    args = parse_args()
    if args.command == "annotations":
        run_annotations(args)
        return
    if args.command == "mask":
        run_mask(args)
        return
    raise ValueError(f"Unsupported command: {args.command}")


if __name__ == "__main__":
    main()
