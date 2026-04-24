#!/biosoftware/miniconda/envs/snakemake_env/bin/python
import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pgta.core.logging import setup_logger


DECISION_RANK = {"PASS": 0, "WARN": 1, "FAIL": 2}
DECISION_COLOR = {"PASS": "#16a34a", "WARN": "#d97706", "FAIL": "#dc2626"}


def parse_args():
    parser = argparse.ArgumentParser(description="Aggregate per-sample baseline BAM QC outputs.")
    parser.add_argument("--qc-tsvs", nargs="+", required=True)
    parser.add_argument("--profile-tsvs", nargs="+", default=[])
    parser.add_argument("--summary-output", required=True)
    parser.add_argument("--pass-samples-output", required=True)
    parser.add_argument("--retained-samples-output", required=True)
    parser.add_argument("--outlier-samples-output", required=True)
    parser.add_argument("--report-output", required=True)
    parser.add_argument("--figures-dir", required=True)
    parser.add_argument("--log", default="", help="Optional log file path")
    return parser.parse_args()


def natural_key(value):
    parts = []
    token = ""
    is_digit = None
    for char in value:
        char_is_digit = char.isdigit()
        if is_digit is None or char_is_digit == is_digit:
            token += char
        else:
            parts.append(int(token) if is_digit else token)
            token = char
        is_digit = char_is_digit
    if token:
        parts.append(int(token) if is_digit else token)
    return tuple((0, item) if isinstance(item, int) else (1, item) for item in parts)


def sample_id_from_bam(target_bam):
    return Path(target_bam).name.replace(".sorted.bam", "")


def infer_batch_group(sample_id):
    if sample_id.startswith("E") and sample_id[1:].isdigit():
        return "batch5_E"
    if sample_id.isdigit():
        return "batch3_num"
    if len(sample_id) == 1 and sample_id in "ABCDEFGH":
        return "batch2_AH"
    return "other"


def _to_float(value):
    try:
        return float(value)
    except Exception:  # noqa: BLE001
        return float("nan")


def _median_text(values, digits=4):
    finite = [float(v) for v in values if np.isfinite(v)]
    if not finite:
        return "NA"
    return f"{float(np.median(finite)):.{digits}f}"


def load_profile_map(profile_tsvs):
    mapping = {}
    for profile_tsv in profile_tsvs:
        path = Path(profile_tsv)
        sample_id = path.parent.name
        mapping[sample_id] = str(path)
    return mapping


def load_qc_rows(qc_tsvs, profile_map):
    rows = []
    for qc_tsv in qc_tsvs:
        qc_path = Path(qc_tsv)
        with open(qc_path, "r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            file_rows = list(reader)
        if len(file_rows) != 1:
            raise ValueError(f"Expected exactly 1 row in {qc_path}, found {len(file_rows)}")
        row = file_rows[0]
        sample_id = sample_id_from_bam(row["target_bam"])
        row["sample_id"] = sample_id
        row["batch_group"] = infer_batch_group(sample_id)
        row["source_tsv"] = str(qc_path)
        row["profile_tsv"] = profile_map.get(sample_id, "")
        rows.append(row)
    return sorted(
        rows,
        key=lambda item: (
            DECISION_RANK.get(item["qc_decision"].strip().upper(), 99),
            natural_key(item["sample_id"]),
        ),
    )


def write_summary(summary_output, rows):
    summary_path = Path(summary_output)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "sample_id",
        "batch_group",
        "bin_size",
        "qc_decision",
        "qc_reason",
        "mapped_fragments",
        "usable_bins",
        "zero_bin_fraction",
        "bin_cv",
        "adjacent_diff_mad",
        "gini_coefficient",
        "pearson_r",
        "spearman_r",
        "median_abs_z",
        "outlier_frac_abs_z_gt_3",
        "outlier_frac_abs_z_gt_5",
        "gc_fraction_mean",
        "gc_signal_pearson_r",
        "gc_signal_spearman_r",
        "gc_signal_slope",
        "raw_gc_signal_pearson_r",
        "raw_gc_signal_spearman_r",
        "raw_gc_signal_slope",
        "gc_correction_method",
        "gc_correction_applied",
        "gc_correction_valid_bins",
        "target_bam",
        "source_tsv",
        "profile_tsv",
    ]
    with open(summary_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def write_pass_samples(pass_samples_output, rows):
    output_path = Path(pass_samples_output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    passed = [row["sample_id"] for row in rows if row["qc_decision"].strip().upper() == "PASS"]
    with open(output_path, "w", encoding="utf-8") as handle:
        for sample_id in sorted(passed, key=natural_key):
            handle.write(f"{sample_id}\n")


def write_decision_table(output_path, rows):
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "sample_id",
        "batch_group",
        "qc_decision",
        "qc_reason",
        "mapped_fragments",
        "bin_cv",
        "spearman_r",
        "median_abs_z",
        "outlier_frac_abs_z_gt_3",
        "target_bam",
        "source_tsv",
        "profile_tsv",
    ]
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def build_dataframe(rows):
    df = pd.DataFrame(rows).copy()
    numeric_columns = [
        "bin_size",
        "mapped_fragments",
        "usable_bins",
        "zero_bin_fraction",
        "bin_cv",
        "adjacent_diff_mad",
        "gini_coefficient",
        "pearson_r",
        "spearman_r",
        "median_abs_z",
        "outlier_frac_abs_z_gt_3",
        "outlier_frac_abs_z_gt_5",
        "gc_fraction_mean",
        "gc_signal_pearson_r",
        "gc_signal_spearman_r",
        "gc_signal_slope",
        "raw_gc_signal_pearson_r",
        "raw_gc_signal_spearman_r",
        "raw_gc_signal_slope",
        "gc_correction_applied",
        "gc_correction_valid_bins",
    ]
    for column in numeric_columns:
        df[column] = pd.to_numeric(df[column], errors="coerce")
    return df


def plot_decision_counts(df, out_png):
    counts = df["qc_decision"].str.upper().value_counts().reindex(["PASS", "WARN", "FAIL"], fill_value=0)
    fig, ax = plt.subplots(figsize=(6.0, 4.5))
    ax.bar(counts.index, counts.values, color=[DECISION_COLOR[item] for item in counts.index])
    ax.set_ylabel("Sample count")
    ax.set_title("Baseline QC decision counts")
    for idx, value in enumerate(counts.values):
        ax.text(idx, value + 0.1, str(int(value)), ha="center", va="bottom", fontsize=10)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def plot_mapped_fragments(df, out_png):
    plot_df = df.copy()
    plot_df["sample_sort_key"] = plot_df["sample_id"].map(lambda item: tuple(natural_key(item)))
    plot_df = plot_df.sort_values(["batch_group", "sample_sort_key"])
    colors = [DECISION_COLOR.get(item, "#6b7280") for item in plot_df["qc_decision"].str.upper()]
    fig, ax = plt.subplots(figsize=(max(10.0, 0.35 * len(plot_df)), 5.5))
    ax.bar(plot_df["sample_id"], plot_df["mapped_fragments"], color=colors)
    ax.set_ylabel("Mapped fragments")
    ax.set_xlabel("Sample")
    ax.set_title("Mapped fragments by baseline sample")
    ax.tick_params(axis="x", rotation=90)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def plot_metric_scatter(df, out_png):
    fig, ax = plt.subplots(figsize=(7.5, 6.0))
    for decision in ["PASS", "WARN", "FAIL"]:
        subset = df[df["qc_decision"].str.upper() == decision]
        if subset.empty:
            continue
        ax.scatter(
            subset["bin_cv"],
            subset["spearman_r"],
            s=50,
            alpha=0.9,
            color=DECISION_COLOR[decision],
            label=decision,
        )
    for _, row in df[df["qc_decision"].str.upper() != "PASS"].iterrows():
        ax.text(row["bin_cv"], row["spearman_r"], row["sample_id"], fontsize=8, ha="left", va="bottom")
    ax.set_xlabel("Bin CV")
    ax.set_ylabel("Spearman r to reference median")
    ax.set_title("QC metric scatter")
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def plot_batch_metric_boxplots(df, out_png):
    metrics = [
        ("mapped_fragments", "Mapped fragments"),
        ("bin_cv", "Bin CV"),
        ("spearman_r", "Spearman r"),
    ]
    groups = [group for group in ["batch2_AH", "batch3_num", "batch5_E", "other"] if group in set(df["batch_group"])]
    fig, axes = plt.subplots(1, len(metrics), figsize=(5.0 * len(metrics), 5.5))
    if len(metrics) == 1:
        axes = [axes]
    for axis, (metric, title) in zip(axes, metrics):
        data = [df.loc[df["batch_group"] == group, metric].dropna().tolist() for group in groups]
        axis.boxplot(data, labels=groups, patch_artist=True)
        axis.set_title(title)
        axis.tick_params(axis="x", rotation=25)
    fig.suptitle("Batch-level QC metric comparison", y=1.02)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)


def load_profile_signals(rows):
    signals = []
    sample_ids = []
    for row in rows:
        profile_tsv = row.get("profile_tsv", "")
        if not profile_tsv:
            continue
        profile_path = Path(profile_tsv)
        if not profile_path.exists():
            continue
        profile_df = pd.read_csv(profile_path, sep="\t")
        if "signal_log2cpm" not in profile_df.columns:
            continue
        signal = pd.to_numeric(profile_df["signal_log2cpm"], errors="coerce").to_numpy(dtype=float)
        signals.append(signal)
        sample_ids.append(row["sample_id"])
    return sample_ids, signals


def plot_sample_correlation_heatmap(rows, out_png):
    sample_ids, signals = load_profile_signals(rows)
    fig, ax = plt.subplots(figsize=(8.0, 7.0))
    if not signals:
        ax.text(0.5, 0.5, "No profile TSV available", ha="center", va="center")
        ax.set_axis_off()
    else:
        matrix = np.vstack(signals)
        corr = np.corrcoef(matrix)
        image = ax.imshow(corr, vmin=0.7, vmax=1.0, cmap="viridis")
        ax.set_xticks(range(len(sample_ids)))
        ax.set_yticks(range(len(sample_ids)))
        ax.set_xticklabels(sample_ids, rotation=90, fontsize=8)
        ax.set_yticklabels(sample_ids, fontsize=8)
        ax.set_title("Sample correlation heatmap")
        fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def write_markdown_report(report_output, rows, retained_rows, outlier_rows, figures_dir):
    report_path = Path(report_output)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    df = build_dataframe(rows)

    total = len(rows)
    pass_count = int((df["qc_decision"].str.upper() == "PASS").sum())
    warn_count = int((df["qc_decision"].str.upper() == "WARN").sum())
    fail_count = int((df["qc_decision"].str.upper() == "FAIL").sum())
    bin_sizes = sorted({str(int(value)) for value in df["bin_size"].dropna().unique()})
    correction_methods = sorted({str(item) for item in df.get("gc_correction_method", pd.Series(dtype=str)).dropna().unique() if str(item)})
    batch_lines = []
    for batch_name in ["batch2_AH", "batch3_num", "batch5_E", "other"]:
        batch_df = df[df["batch_group"] == batch_name]
        if batch_df.empty:
            continue
        batch_pass = int((batch_df["qc_decision"].str.upper() == "PASS").sum())
        batch_warn = int((batch_df["qc_decision"].str.upper() == "WARN").sum())
        batch_fail = int((batch_df["qc_decision"].str.upper() == "FAIL").sum())
        batch_lines.append(f"- {batch_name}: {len(batch_df)} samples, PASS/WARN/FAIL = {batch_pass}/{batch_warn}/{batch_fail}")

    lines = []
    lines.append("# Baseline QC Report")
    lines.append("")
    lines.append("## 1. Summary")
    lines.append("")
    lines.append(f"- sample_count: {total}")
    lines.append(f"- pass_warn_fail: {pass_count}/{warn_count}/{fail_count}")
    lines.append(f"- qc_bin_size_bp: {', '.join(bin_sizes) if bin_sizes else 'NA'}")
    lines.append(f"- retained_baseline_count: {len(retained_rows)}")
    lines.append(f"- outlier_or_warning_count: {len(outlier_rows)}")
    lines.append(f"- gc_correction_method: {', '.join(correction_methods) if correction_methods else 'NA'}")
    lines.append("")
    lines.append("## 2. Batch Distribution")
    lines.append("")
    lines.extend(batch_lines or ["- no_batch_information"])
    lines.append("")
    lines.append("## 3. Sample-level Results")
    lines.append("")
    lines.append("| sample_id | batch_group | decision | mapped_fragments | bin_cv | spearman_r | median_abs_z | outlier_frac_abs_z_gt_3 | reason |")
    lines.append("| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- |")
    for row in rows:
        lines.append(
            "| {sample_id} | {batch_group} | {decision} | {mapped_fragments} | {bin_cv:.4f} | {spearman_r:.4f} | {median_abs_z:.4f} | {outlier_frac:.4f} | {reason} |".format(
                sample_id=row["sample_id"],
                batch_group=row["batch_group"],
                decision=row["qc_decision"],
                mapped_fragments=int(_to_float(row["mapped_fragments"])),
                bin_cv=_to_float(row["bin_cv"]),
                spearman_r=_to_float(row["spearman_r"]),
                median_abs_z=_to_float(row["median_abs_z"]),
                outlier_frac=_to_float(row["outlier_frac_abs_z_gt_3"]),
                reason=row["qc_reason"],
            )
        )
    lines.append("")
    lines.append("## 4. Aggregate Metrics")
    lines.append("")
    lines.append(f"- mapped_fragments_median: {_median_text(df['mapped_fragments'], 0)}")
    lines.append(f"- bin_cv_median: {_median_text(df['bin_cv'])}")
    lines.append(f"- spearman_r_median: {_median_text(df['spearman_r'])}")
    lines.append(f"- median_abs_z_median: {_median_text(df['median_abs_z'])}")
    lines.append(f"- outlier_frac_abs_z_gt_3_median: {_median_text(df['outlier_frac_abs_z_gt_3'])}")
    lines.append(f"- raw_gc_signal_slope_median: {_median_text(df['raw_gc_signal_slope'])}")
    lines.append(f"- corrected_gc_signal_slope_median: {_median_text(df['gc_signal_slope'])}")
    lines.append("")
    lines.append("## 5. Retained / Outlier Lists")
    lines.append("")
    lines.append(f"- retained: {', '.join(row['sample_id'] for row in retained_rows) if retained_rows else 'None'}")
    lines.append(f"- outlier_or_warning: {', '.join(row['sample_id'] for row in outlier_rows) if outlier_rows else 'None'}")
    lines.append("")
    lines.append("## 6. Figures")
    lines.append("")
    lines.append(f"- figures dir: `{figures_dir}`")
    lines.append("- `qc_decision_counts.png`")
    lines.append("- `mapped_fragments_by_sample.png`")
    lines.append("- `metric_scatter.png`")
    lines.append("- `batch_metric_boxplots.png`")
    lines.append("- `sample_correlation_heatmap.png`")
    lines.append("")
    lines.append("Note: retained baseline currently follows PASS samples only; WARN and FAIL are written to outlier_samples.tsv for manual review.")
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    args = parse_args()
    logger = setup_logger("aggregate_baseline_qc", args.log or None)

    figures_dir = Path(args.figures_dir)
    figures_dir.mkdir(parents=True, exist_ok=True)

    profile_map = load_profile_map(args.profile_tsvs)
    rows = load_qc_rows(args.qc_tsvs, profile_map)
    retained_rows = [row for row in rows if row["qc_decision"].strip().upper() == "PASS"]
    outlier_rows = [row for row in rows if row["qc_decision"].strip().upper() != "PASS"]
    df = build_dataframe(rows)

    write_summary(args.summary_output, rows)
    write_pass_samples(args.pass_samples_output, retained_rows)
    write_decision_table(args.retained_samples_output, retained_rows)
    write_decision_table(args.outlier_samples_output, outlier_rows)

    plot_decision_counts(df, figures_dir / "qc_decision_counts.png")
    plot_mapped_fragments(df, figures_dir / "mapped_fragments_by_sample.png")
    plot_metric_scatter(df, figures_dir / "metric_scatter.png")
    plot_batch_metric_boxplots(df, figures_dir / "batch_metric_boxplots.png")
    plot_sample_correlation_heatmap(rows, figures_dir / "sample_correlation_heatmap.png")

    write_markdown_report(args.report_output, rows, retained_rows, outlier_rows, figures_dir)

    logger.info("baseline QC aggregated: samples=%d summary=%s", len(rows), args.summary_output)
    logger.info("baseline QC retained table written: %s", args.retained_samples_output)
    logger.info("baseline QC outlier table written: %s", args.outlier_samples_output)
    logger.info("baseline QC markdown report written: %s", args.report_output)


if __name__ == "__main__":
    main()
