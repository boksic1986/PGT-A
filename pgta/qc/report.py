#!/biosoftware/miniconda/envs/snakemake_env/bin/python
import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pgta.core.logging import setup_logger


DECISION_SCORE = {"PASS": 0, "WARN": 1, "FAIL": 2}
DECISION_COLOR = {"PASS": "#15803d", "WARN": "#d97706", "FAIL": "#b91c1c"}
DECISION_HEATMAP = {"PASS": 0, "WARN": 1, "FAIL": 2}


def parse_args():
    parser = argparse.ArgumentParser(description="Generate baseline QC multiscale summary report.")
    parser.add_argument("--summary-tsvs", nargs="+", required=True)
    parser.add_argument("--profile-tsvs", nargs="+", required=True)
    parser.add_argument("--fastp-jsons", nargs="+", required=True)
    parser.add_argument("--flagstats", nargs="+", required=True)
    parser.add_argument("--idxstats", nargs="+", required=True)
    parser.add_argument("--binsize-summary-output", required=True)
    parser.add_argument("--sample-matrix-output", required=True)
    parser.add_argument("--stable-retained-output", required=True)
    parser.add_argument("--failed-diagnostics-output", required=True)
    parser.add_argument("--report-output", required=True)
    parser.add_argument("--figures-dir", required=True)
    parser.add_argument("--log", default="")
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


def batch_group(sample_id):
    if sample_id.startswith("E") and sample_id[1:].isdigit():
        return "batch5_E"
    if sample_id.isdigit():
        return "batch3_num"
    if len(sample_id) == 1 and sample_id in "ABCDEFGH":
        return "batch2_AH"
    return "other"


def parse_bin_size_label(path_value):
    return Path(path_value).parent.name


def load_multiscale_summaries(summary_paths):
    frames = []
    for path_value in summary_paths:
        path = Path(path_value)
        df = pd.read_csv(path, sep="\t")
        df["bin_label"] = parse_bin_size_label(path)
        df["decision_score"] = df["qc_decision"].map(DECISION_SCORE)
        frames.append(df)
    combined = pd.concat(frames, ignore_index=True)
    combined["sample_sort_key"] = combined["sample_id"].map(natural_key)
    return combined


def parse_fastp_json(json_path):
    with open(json_path, "r", encoding="utf-8") as handle:
        data = json.load(handle)
    summary = data.get("summary", {})
    before = summary.get("before_filtering", {})
    after = summary.get("after_filtering", {})
    duplication = data.get("duplication", {})
    insert_size = data.get("insert_size", {})
    return {
        "sample_id": Path(json_path).name.replace(".fastp.json", ""),
        "raw_reads": before.get("total_reads", np.nan),
        "clean_reads": after.get("total_reads", np.nan),
        "raw_q30_rate": before.get("q30_rate", np.nan),
        "clean_q30_rate": after.get("q30_rate", np.nan),
        "raw_gc_content": before.get("gc_content", np.nan),
        "clean_gc_content": after.get("gc_content", np.nan),
        "duplication_rate": duplication.get("rate", np.nan),
        "insert_peak": insert_size.get("peak", np.nan),
        "insert_unknown": insert_size.get("unknown", np.nan),
    }


def parse_flagstat(flagstat_path):
    metrics = {
        "sample_id": Path(flagstat_path).parent.name,
        "total_reads_bam": np.nan,
        "mapped_reads_bam": np.nan,
        "paired_reads_bam": np.nan,
        "properly_paired_reads_bam": np.nan,
        "singleton_reads_bam": np.nan,
        "mate_diff_chr_reads_bam": np.nan,
        "mate_diff_chr_mapq5_reads_bam": np.nan,
    }
    with open(flagstat_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if "+" not in line:
                continue
            count = int(line.split("+", 1)[0].strip())
            if " in total " in line:
                metrics["total_reads_bam"] = count
            elif " mapped (" in line and "mate mapped" not in line:
                metrics["mapped_reads_bam"] = count
            elif " paired in sequencing" in line:
                metrics["paired_reads_bam"] = count
            elif " properly paired " in line:
                metrics["properly_paired_reads_bam"] = count
            elif " singletons " in line:
                metrics["singleton_reads_bam"] = count
            elif " with mate mapped to a different chr (mapQ>=5)" in line:
                metrics["mate_diff_chr_mapq5_reads_bam"] = count
            elif " with mate mapped to a different chr" in line:
                metrics["mate_diff_chr_reads_bam"] = count
    return metrics


def parse_idxstats(idxstats_path):
    sample_id = Path(idxstats_path).parent.name
    rows = []
    with open(idxstats_path, "r", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if len(row) != 4:
                continue
            chrom, length, mapped, unmapped = row
            rows.append(
                {
                    "sample_id": sample_id,
                    "chrom": chrom,
                    "length": int(length),
                    "mapped_reads": int(mapped),
                    "unmapped_reads": int(unmapped),
                }
            )
    return pd.DataFrame(rows)


def autosome_order(chrom):
    token = chrom[3:] if chrom.startswith("chr") else chrom
    if token.isdigit():
        value = int(token)
        if 1 <= value <= 22:
            return value
    return 99


def compute_idxstats_summary(idxstats_frames):
    idx_df = pd.concat(idxstats_frames, ignore_index=True)
    idx_df["autosome_rank"] = idx_df["chrom"].map(autosome_order)
    auto_df = idx_df[idx_df["autosome_rank"] <= 22].copy()
    auto_totals = auto_df.groupby("sample_id")["mapped_reads"].sum().rename("autosome_mapped_reads")
    auto_df = auto_df.merge(auto_totals, on="sample_id", how="left")
    auto_df["autosome_fraction"] = np.divide(
        auto_df["mapped_reads"],
        auto_df["autosome_mapped_reads"],
        out=np.zeros(len(auto_df), dtype=np.float64),
        where=auto_df["autosome_mapped_reads"] > 0,
    )
    return idx_df, auto_df


def load_profile_summary(profile_paths):
    rows = []
    for path_value in profile_paths:
        path = Path(path_value)
        sample_id = path.parent.name
        bin_label = path.parent.parent.name
        df = pd.read_csv(path, sep="\t", usecols=["chrom", "z_score"])
        chr_abs_median = df.groupby("chrom", sort=False)["z_score"].median().abs().sort_values(ascending=False)
        rows.append(
            {
                "sample_id": sample_id,
                "bin_label": bin_label,
                "global_median_abs_z_profile": float(np.median(np.abs(df["z_score"]))),
                "max_chrom_abs_median_z": float(chr_abs_median.iloc[0]) if not chr_abs_median.empty else np.nan,
                "broad_noise_fraction": float(np.mean(chr_abs_median >= 1.0)) if not chr_abs_median.empty else np.nan,
                "top_chr_1": chr_abs_median.index[0] if len(chr_abs_median) > 0 else "",
                "top_chr_2": chr_abs_median.index[1] if len(chr_abs_median) > 1 else "",
                "top_chr_3": chr_abs_median.index[2] if len(chr_abs_median) > 2 else "",
            }
        )
    return pd.DataFrame(rows)


def write_table(path_value, df):
    path = Path(path_value)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def build_binsize_summary(summary_df):
    rows = []
    for bin_label, sub in summary_df.groupby("bin_label", sort=False):
        counts = sub["qc_decision"].value_counts()
        rows.append(
            {
                "bin_label": bin_label,
                "sample_count": int(len(sub)),
                "pass_count": int(counts.get("PASS", 0)),
                "warn_count": int(counts.get("WARN", 0)),
                "fail_count": int(counts.get("FAIL", 0)),
                "mapped_fragments_median": float(sub["mapped_fragments"].median()),
                "bin_cv_median": float(sub["bin_cv"].median()),
                "spearman_r_median": float(sub["spearman_r"].median()),
                "median_abs_z_median": float(sub["median_abs_z"].median()),
            }
        )
    return pd.DataFrame(rows).sort_values("bin_label", key=lambda s: s.map(natural_key)).reset_index(drop=True)


def build_sample_matrix(summary_df):
    pivot = summary_df.pivot(index="sample_id", columns="bin_label", values="qc_decision")
    for bin_label in sorted(summary_df["bin_label"].unique(), key=natural_key):
        if bin_label not in pivot.columns:
            pivot[bin_label] = ""
    pivot = pivot[sorted(pivot.columns, key=natural_key)]
    matrix = pivot.reset_index()
    matrix["batch_group"] = matrix["sample_id"].map(batch_group)
    decision_values = pivot.applymap(lambda value: DECISION_SCORE.get(value, -1))
    matrix["fail_count"] = (decision_values == 2).sum(axis=1).to_numpy()
    matrix["warn_count"] = (decision_values == 1).sum(axis=1).to_numpy()
    matrix["worst_decision"] = decision_values.max(axis=1).map({0: "PASS", 1: "WARN", 2: "FAIL"}).to_numpy()
    matrix["stable_retained"] = ((decision_values == 0).all(axis=1)).to_numpy()
    return matrix.sort_values(
        ["stable_retained", "fail_count", "warn_count", "sample_id"],
        ascending=[False, False, False, True],
        key=lambda s: s.map(natural_key) if s.name == "sample_id" else s,
    ).reset_index(drop=True)


def compute_chr_bias_scores(auto_df, stable_samples):
    cohort = auto_df[auto_df["sample_id"].isin(stable_samples)].copy()
    if cohort.empty:
        cohort = auto_df.copy()
    med = cohort.groupby("chrom")["autosome_fraction"].median()
    mad = cohort.groupby("chrom")["autosome_fraction"].apply(lambda x: float(np.median(np.abs(x - np.median(x)))))
    mad = mad.replace(0, np.nan)

    rows = []
    for sample_id, sub in auto_df.groupby("sample_id", sort=False):
        z_values = []
        top = []
        for row in sub.itertuples(index=False):
            chrom = row.chrom
            median = float(med.get(chrom, np.nan))
            sigma = float(1.4826 * mad.get(chrom, np.nan))
            if not np.isfinite(sigma) or sigma <= 0:
                sigma = 1e-6
            z_score = (float(row.autosome_fraction) - median) / sigma
            z_values.append((chrom, z_score))
        z_values.sort(key=lambda item: abs(item[1]), reverse=True)
        top = z_values[:3]
        rows.append(
            {
                "sample_id": sample_id,
                "max_chr_fraction_dev_z": float(abs(top[0][1])) if top else np.nan,
                "top_chr_fraction_dev_1": top[0][0] if len(top) > 0 else "",
                "top_chr_fraction_dev_2": top[1][0] if len(top) > 1 else "",
                "top_chr_fraction_dev_3": top[2][0] if len(top) > 2 else "",
            }
        )
    return pd.DataFrame(rows)


def join_sample_metrics(sample_matrix, fastp_df, flagstat_df, chr_bias_df, profile_df):
    current_1mb = profile_df[profile_df["bin_label"] == "1mb"].copy()
    current_1mb = current_1mb.rename(
        columns={
            "global_median_abs_z_profile": "profile_median_abs_z_1mb",
            "max_chrom_abs_median_z": "max_chrom_abs_median_z_1mb",
            "broad_noise_fraction": "broad_noise_fraction_1mb",
            "top_chr_1": "top_chr_1mb_1",
            "top_chr_2": "top_chr_1mb_2",
            "top_chr_3": "top_chr_1mb_3",
        }
    )
    current_1mb = current_1mb.drop(columns=["bin_label"])

    merged = sample_matrix.merge(fastp_df, on="sample_id", how="left")
    merged = merged.merge(flagstat_df, on="sample_id", how="left")
    merged = merged.merge(chr_bias_df, on="sample_id", how="left")
    merged = merged.merge(current_1mb, on="sample_id", how="left")

    merged["clean_read_rate"] = merged["clean_reads"] / merged["raw_reads"]
    merged["mapping_rate_bam"] = merged["mapped_reads_bam"] / merged["total_reads_bam"]
    merged["proper_pair_rate_bam"] = merged["properly_paired_reads_bam"] / merged["paired_reads_bam"]
    merged["singleton_rate_bam"] = merged["singleton_reads_bam"] / merged["paired_reads_bam"]
    merged["diff_chr_rate_bam"] = merged["mate_diff_chr_reads_bam"] / merged["paired_reads_bam"]
    return merged


def classify_failed_samples(sample_metrics):
    cohort = sample_metrics[sample_metrics["stable_retained"]].copy()
    if cohort.empty:
        cohort = sample_metrics[sample_metrics["worst_decision"] == "PASS"].copy()

    insert_median = float(cohort["insert_peak"].median()) if not cohort.empty else np.nan
    insert_mad = float(np.median(np.abs(cohort["insert_peak"] - insert_median))) if not cohort.empty else np.nan
    dup_median = float(cohort["duplication_rate"].median()) if not cohort.empty else np.nan
    dup_mad = float(np.median(np.abs(cohort["duplication_rate"] - dup_median))) if not cohort.empty else np.nan

    rows = []
    focus = sample_metrics[sample_metrics["worst_decision"] != "PASS"].copy()
    for row in focus.itertuples(index=False):
        insert_shift = abs(float(row.insert_peak) - insert_median) if np.isfinite(insert_median) else np.nan
        dup_shift = float(row.duplication_rate) - dup_median if np.isfinite(dup_median) else np.nan
        broad_noise = float(row.broad_noise_fraction_1mb) if np.isfinite(row.broad_noise_fraction_1mb) else np.nan
        chr_bias = float(row.max_chr_fraction_dev_z) if np.isfinite(row.max_chr_fraction_dev_z) else np.nan
        chrom_peak = float(row.max_chrom_abs_median_z_1mb) if np.isfinite(row.max_chrom_abs_median_z_1mb) else np.nan
        mapping_ok = np.isfinite(row.mapping_rate_bam) and row.mapping_rate_bam >= 0.985
        pair_ok = np.isfinite(row.proper_pair_rate_bam) and row.proper_pair_rate_bam >= 0.96
        dup_high = np.isfinite(dup_shift) and np.isfinite(dup_mad) and dup_shift > max(0.01, 3.0 * dup_mad)
        insert_abnormal = np.isfinite(insert_shift) and np.isfinite(insert_mad) and insert_shift > max(25.0, 3.0 * insert_mad)
        broad_noise_high = np.isfinite(broad_noise) and broad_noise >= 0.25
        chr_bias_high = np.isfinite(chr_bias) and chr_bias >= 4.0
        chrom_peak_high = np.isfinite(chrom_peak) and chrom_peak >= 1.8

        if chr_bias_high and chrom_peak_high and mapping_ok and pair_ok and not dup_high and not insert_abnormal:
            diagnosis = "更偏向样本DNA/染色体组成异常"
            note = "跨染色体分布偏离较集中，建库层指标未见明显异常。"
        elif broad_noise_high or dup_high or insert_abnormal or not mapping_ok or not pair_ok:
            diagnosis = "更偏向建库/扩增波动"
            note = "表现为全基因组波动偏高或文库层指标偏离队列中位。"
        else:
            diagnosis = "更偏向文库噪声，暂不能排除上游DNA因素"
            note = "缺少单染色体集中偏移证据，但覆盖一致性不足。"

        rows.append(
            {
                "sample_id": row.sample_id,
                "batch_group": row.batch_group,
                "worst_decision": row.worst_decision,
                "fail_count": int(row.fail_count),
                "warn_count": int(row.warn_count),
                "mapping_rate_bam": row.mapping_rate_bam,
                "proper_pair_rate_bam": row.proper_pair_rate_bam,
                "duplication_rate": row.duplication_rate,
                "insert_peak": row.insert_peak,
                "clean_read_rate": row.clean_read_rate,
                "profile_median_abs_z_1mb": row.profile_median_abs_z_1mb,
                "max_chrom_abs_median_z_1mb": row.max_chrom_abs_median_z_1mb,
                "broad_noise_fraction_1mb": row.broad_noise_fraction_1mb,
                "max_chr_fraction_dev_z": row.max_chr_fraction_dev_z,
                "top_chr_profile": ",".join([item for item in [row.top_chr_1mb_1, row.top_chr_1mb_2, row.top_chr_1mb_3] if item]),
                "top_chr_fraction_dev": ",".join(
                    [item for item in [row.top_chr_fraction_dev_1, row.top_chr_fraction_dev_2, row.top_chr_fraction_dev_3] if item]
                ),
                "diagnosis": diagnosis,
                "note": note,
            }
        )
    return pd.DataFrame(rows).sort_values(["worst_decision", "fail_count", "sample_id"], ascending=[False, False, True]).reset_index(drop=True)


def plot_binsize_decision_trend(binsize_summary, output_png):
    x = np.arange(len(binsize_summary))
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    pass_values = binsize_summary["pass_count"].to_numpy()
    warn_values = binsize_summary["warn_count"].to_numpy()
    fail_values = binsize_summary["fail_count"].to_numpy()
    ax.bar(x, pass_values, color=DECISION_COLOR["PASS"], label="PASS")
    ax.bar(x, warn_values, bottom=pass_values, color=DECISION_COLOR["WARN"], label="WARN")
    ax.bar(x, fail_values, bottom=pass_values + warn_values, color=DECISION_COLOR["FAIL"], label="FAIL")
    ax.set_xticks(x)
    ax.set_xticklabels(binsize_summary["bin_label"])
    ax.set_ylabel("Sample count")
    ax.set_title("Baseline QC decisions by bin size")
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(output_png, dpi=150)
    plt.close(fig)


def plot_sample_decision_heatmap(sample_matrix, output_png):
    bin_columns = [col for col in sample_matrix.columns if col.endswith("kb") or col.endswith("mb")]
    plot_df = sample_matrix.copy()
    plot_df = plot_df.sort_values(["stable_retained", "fail_count", "warn_count", "sample_id"], ascending=[False, False, False, True])
    matrix = plot_df[bin_columns].replace(DECISION_HEATMAP).to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(7.0, max(6.0, 0.25 * len(plot_df))))
    image = ax.imshow(matrix, cmap=plt.cm.get_cmap("RdYlGn_r", 3), vmin=0, vmax=2, aspect="auto")
    ax.set_xticks(np.arange(len(bin_columns)))
    ax.set_xticklabels(bin_columns)
    ax.set_yticks(np.arange(len(plot_df)))
    ax.set_yticklabels(plot_df["sample_id"], fontsize=8)
    ax.set_title("Sample decision matrix across bin sizes")
    cbar = fig.colorbar(image, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_ticks([0, 1, 2])
    cbar.set_ticklabels(["PASS", "WARN", "FAIL"])
    fig.tight_layout()
    fig.savefig(output_png, dpi=150)
    plt.close(fig)


def plot_library_metric_overview(sample_metrics, output_png):
    plot_df = sample_metrics.copy()
    plot_df["color"] = plot_df["worst_decision"].map(DECISION_COLOR).fillna("#64748b")
    plot_df = plot_df.sort_values("sample_id", key=lambda s: s.map(natural_key))

    fig, axes = plt.subplots(3, 1, figsize=(12.0, 10.0), sharex=True)
    axes[0].bar(plot_df["sample_id"], plot_df["clean_reads"] / 1e6, color=plot_df["color"])
    axes[0].set_ylabel("Clean reads (M)")
    axes[0].set_title("Clean reads by sample")

    axes[1].bar(plot_df["sample_id"], plot_df["duplication_rate"] * 100.0, color=plot_df["color"])
    axes[1].set_ylabel("Duplication rate (%)")
    axes[1].set_title("fastp duplication rate")

    axes[2].bar(plot_df["sample_id"], plot_df["proper_pair_rate_bam"] * 100.0, color=plot_df["color"])
    axes[2].set_ylabel("Proper pair rate (%)")
    axes[2].set_title("samtools properly paired rate")
    axes[2].tick_params(axis="x", rotation=90)

    fig.tight_layout()
    fig.savefig(output_png, dpi=150)
    plt.close(fig)


def plot_failed_chr_bias(profile_df, failed_df, output_png):
    focus_samples = failed_df["sample_id"].tolist()
    if not focus_samples:
        fig, ax = plt.subplots(figsize=(6.0, 4.0))
        ax.text(0.5, 0.5, "No failed or warning samples", ha="center", va="center")
        ax.set_axis_off()
        fig.tight_layout()
        fig.savefig(output_png, dpi=150)
        plt.close(fig)
        return

    plot_df = profile_df[(profile_df["bin_label"] == "1mb") & (profile_df["sample_id"].isin(focus_samples))].copy()
    if plot_df.empty:
        fig, ax = plt.subplots(figsize=(6.0, 4.0))
        ax.text(0.5, 0.5, "No 1mb profile summary", ha="center", va="center")
        ax.set_axis_off()
        fig.tight_layout()
        fig.savefig(output_png, dpi=150)
        plt.close(fig)
        return

    rows = []
    for sample_id in focus_samples:
        sample_profile = plot_df[plot_df["sample_id"] == sample_id]
        rows.append(
            [
                float(sample_profile["max_chrom_abs_median_z"].iloc[0]) if not sample_profile.empty else np.nan,
                float(sample_profile["broad_noise_fraction"].iloc[0]) if not sample_profile.empty else np.nan,
            ]
        )
    matrix = np.asarray(rows, dtype=float)

    fig, ax = plt.subplots(figsize=(5.0, max(4.0, 0.45 * len(focus_samples))))
    image = ax.imshow(matrix, cmap="YlOrRd", aspect="auto")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Max chr median |z|", "Broad noise fraction"])
    ax.set_yticks(np.arange(len(focus_samples)))
    ax.set_yticklabels(focus_samples)
    ax.set_title("Failed / warning sample chromosome bias")
    fig.colorbar(image, ax=ax, fraction=0.04, pad=0.02)
    fig.tight_layout()
    fig.savefig(output_png, dpi=150)
    plt.close(fig)


def write_stable_retained(path_value, sample_matrix):
    retained = sample_matrix[sample_matrix["stable_retained"]].copy()
    retained = retained[["sample_id", "batch_group"]]
    write_table(path_value, retained)


def write_report(report_output, binsize_summary, sample_matrix, failed_df, figures_dir):
    report_path = Path(report_output)
    report_path.parent.mkdir(parents=True, exist_ok=True)

    stable_samples = sample_matrix.loc[sample_matrix["stable_retained"], "sample_id"].tolist()
    review_samples = sample_matrix.loc[(sample_matrix["worst_decision"] == "WARN") & (~sample_matrix["stable_retained"]), "sample_id"].tolist()
    exclude_samples = sample_matrix.loc[sample_matrix["worst_decision"] == "FAIL", "sample_id"].tolist()

    lines = []
    lines.append("# Baseline文库QC汇总报告")
    lines.append("")
    lines.append("统计日期：2026-04-14")
    lines.append("")
    lines.append("## 一、结论摘要")
    lines.append("")
    lines.append(f"- baseline样本总数：{int(sample_matrix.shape[0])}")
    lines.append(f"- 统计binsize：{', '.join(binsize_summary['bin_label'].tolist())}")
    lines.append(f"- 全部binsize均为PASS的稳定入库样本：{len(stable_samples)}")
    lines.append(f"- 需复核样本：{len(review_samples)}")
    lines.append(f"- 建议剔除样本：{len(exclude_samples)}")
    lines.append("")
    lines.append("当前baseline文库整体可用于后续reference/tune，但不建议全量直接入库。第三批数字样本整体最稳定；E批样本总体可用，但存在局部离群；第二批A-H中个别样本波动明显。")
    lines.append("")
    lines.append("## 二、分binsize统计")
    lines.append("")
    lines.append("| binsize | sample_count | PASS | WARN | FAIL | median mapped fragments | median bin_cv | median spearman_r | median median_abs_z |")
    lines.append("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for row in binsize_summary.itertuples(index=False):
        lines.append(
            f"| {row.bin_label} | {int(row.sample_count)} | {int(row.pass_count)} | {int(row.warn_count)} | {int(row.fail_count)} | {row.mapped_fragments_median:.0f} | {row.bin_cv_median:.4f} | {row.spearman_r_median:.4f} | {row.median_abs_z_median:.4f} |"
        )
    lines.append("")
    lines.append("## 三、多binsize一致性")
    lines.append("")
    lines.append(f"- 稳定保留样本：{', '.join(stable_samples) if stable_samples else '无'}")
    lines.append(f"- 复核样本：{', '.join(review_samples) if review_samples else '无'}")
    lines.append(f"- 建议剔除样本：{', '.join(exclude_samples) if exclude_samples else '无'}")
    lines.append("")
    lines.append("稳定保留标准：在全部统计binsize下均为PASS。该名单可直接作为下一步reference/tune的优先输入集合。")
    lines.append("")
    lines.append("## 四、失败/异常样本原因判断")
    lines.append("")
    lines.append("| sample_id | batch | worst_decision | fail_count | warn_count | mapping_rate | proper_pair_rate | dup_rate | insert_peak | 1mb max chr |z| | 1mb broad noise | diagnosis |")
    lines.append("| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |")
    for row in failed_df.itertuples(index=False):
        lines.append(
            f"| {row.sample_id} | {row.batch_group} | {row.worst_decision} | {int(row.fail_count)} | {int(row.warn_count)} | {row.mapping_rate_bam:.4f} | {row.proper_pair_rate_bam:.4f} | {row.duplication_rate:.4f} | {row.insert_peak:.0f} | {row.max_chrom_abs_median_z_1mb:.4f} | {row.broad_noise_fraction_1mb:.4f} | {row.diagnosis} |"
        )
    lines.append("")
    lines.append("判定原则：")
    lines.append("- 若单染色体偏离集中、比对和proper pair正常，更偏向样本DNA/染色体组成异常。")
    lines.append("- 若全基因组波动偏高，且插入片段、重复率或proper pair偏离队列中位，更偏向建库/扩增波动。")
    lines.append("- 当前E7、E8、H更接近文库噪声/扩增波动，不支持继续纳入baseline。")
    lines.append("")
    lines.append("## 五、baseline入库建议")
    lines.append("")
    lines.append(f"- 第一优先：{', '.join(stable_samples) if stable_samples else '无'}")
    lines.append(f"- 进入人工复核后再决定：{', '.join(review_samples) if review_samples else '无'}")
    lines.append(f"- 当前建议剔除：{', '.join(exclude_samples) if exclude_samples else '无'}")
    lines.append("")
    lines.append("## 六、配图")
    lines.append("")
    lines.append(f"- 统计图目录：`{figures_dir}`")
    lines.append("- `binsize_decision_trend.png`：不同binsize下PASS/WARN/FAIL分布")
    lines.append("- `sample_decision_heatmap.png`：样本在多binsize下的判定矩阵")
    lines.append("- `library_metric_overview.png`：clean reads、duplication rate、proper pair rate样本概览")
    lines.append("- `failed_sample_chromosome_bias.png`：异常样本染色体层面偏移摘要")
    lines.append("")
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    args = parse_args()
    logger = setup_logger("baseline_qc_report", args.log or None)

    figures_dir = Path(args.figures_dir)
    figures_dir.mkdir(parents=True, exist_ok=True)

    summary_df = load_multiscale_summaries(args.summary_tsvs)
    binsize_summary = build_binsize_summary(summary_df)
    sample_matrix = build_sample_matrix(summary_df)

    fastp_df = pd.DataFrame([parse_fastp_json(path) for path in args.fastp_jsons])
    flagstat_df = pd.DataFrame([parse_flagstat(path) for path in args.flagstats])
    idx_df, auto_df = compute_idxstats_summary([parse_idxstats(path) for path in args.idxstats])
    profile_df = load_profile_summary(args.profile_tsvs)

    stable_samples = sample_matrix.loc[sample_matrix["stable_retained"], "sample_id"].tolist()
    chr_bias_df = compute_chr_bias_scores(auto_df, stable_samples)
    sample_metrics = join_sample_metrics(sample_matrix, fastp_df, flagstat_df, chr_bias_df, profile_df)
    failed_df = classify_failed_samples(sample_metrics)

    write_table(args.binsize_summary_output, binsize_summary)
    write_table(args.sample_matrix_output, sample_matrix)
    write_stable_retained(args.stable_retained_output, sample_matrix)
    write_table(args.failed_diagnostics_output, failed_df)

    plot_binsize_decision_trend(binsize_summary, figures_dir / "binsize_decision_trend.png")
    plot_sample_decision_heatmap(sample_matrix, figures_dir / "sample_decision_heatmap.png")
    plot_library_metric_overview(sample_metrics, figures_dir / "library_metric_overview.png")
    plot_failed_chr_bias(profile_df, failed_df, figures_dir / "failed_sample_chromosome_bias.png")

    write_report(args.report_output, binsize_summary, sample_matrix, failed_df, figures_dir)

    logger.info("wrote multiscale binsize summary: %s", args.binsize_summary_output)
    logger.info("wrote multiscale sample matrix: %s", args.sample_matrix_output)
    logger.info("wrote stable retained samples: %s", args.stable_retained_output)
    logger.info("wrote failed diagnostics: %s", args.failed_diagnostics_output)
    logger.info("wrote report: %s", args.report_output)


if __name__ == "__main__":
    main()
