#!/biosoftware/miniconda/envs/snakemake_env/bin/python
import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


STATUS_ORDER = ["pass", "review", "fail", "no-call"]
STATUS_COLOR = {"pass": "#16a34a", "review": "#d97706", "fail": "#dc2626", "no-call": "#6b7280"}


def parse_args():
    parser = argparse.ArgumentParser(description="Structured run/sample/bin/event QC framework.")
    parser.add_argument("--baseline-summary", required=True)
    parser.add_argument("--profile-tsvs", nargs="*", default=[])
    parser.add_argument("--bin-annotations", required=True)
    parser.add_argument("--combined-mask", required=True)
    parser.add_argument("--run-metadata", default="")
    parser.add_argument("--run-qc-tsv", required=True)
    parser.add_argument("--run-qc-json", required=True)
    parser.add_argument("--sample-qc-tsv", required=True)
    parser.add_argument("--sample-qc-json", required=True)
    parser.add_argument("--bin-qc-tsv", required=True)
    parser.add_argument("--bin-qc-json", required=True)
    parser.add_argument("--event-qc-tsv", required=True)
    parser.add_argument("--event-qc-json", required=True)
    parser.add_argument("--sample-status-figure", required=True)
    parser.add_argument("--event-status-figure", required=True)
    parser.add_argument("--mask-status-figure", required=True)
    parser.add_argument("--metric-scatter-figure", required=True)
    return parser.parse_args()


def sample_status(decision):
    mapping = {"PASS": "pass", "WARN": "review", "FAIL": "fail", "NO_CALL": "no-call"}
    return mapping.get(str(decision).strip().upper(), "review")


def load_profiles(profile_tsvs):
    profile_map = {}
    for path_value in profile_tsvs:
        path = Path(path_value)
        sample_id = path.parent.name
        if not path.exists():
            continue
        df = pd.read_csv(path, sep="\t")
        required = {"chrom", "start", "end", "z_score"}
        if not required.issubset(df.columns):
            continue
        profile_map[sample_id] = df
    return profile_map


def derive_event_rows(summary_df, profile_map, combined_mask_df):
    event_rows = []
    masked_bins = combined_mask_df[combined_mask_df["mask_label"] != "pass"][["chrom", "start", "end"]].copy()
    masked_bins["is_masked"] = True
    for row in summary_df.itertuples(index=False):
        sample_id = row.sample_id
        if sample_id not in profile_map:
            event_rows.append({"sample_id": sample_id, "event_type": "profile_missing", "status": "no-call", "detail": "profile_tsv_missing"})
            continue
        profile_df = profile_map[sample_id].merge(masked_bins, on=["chrom", "start", "end"], how="left")
        profile_df["is_masked"] = profile_df["is_masked"].fillna(False)

        if float(row.spearman_r) < 0.83:
            event_rows.append({"sample_id": sample_id, "event_type": "correlation_drop", "status": "fail", "detail": f"spearman_r={row.spearman_r:.4f}"})
        elif float(row.spearman_r) < 0.90:
            event_rows.append({"sample_id": sample_id, "event_type": "correlation_drop", "status": "review", "detail": f"spearman_r={row.spearman_r:.4f}"})

        if float(row.median_abs_z) >= 1.0:
            severity = "fail" if float(row.median_abs_z) >= 1.5 else "review"
            event_rows.append({"sample_id": sample_id, "event_type": "broad_noise", "status": severity, "detail": f"median_abs_z={row.median_abs_z:.4f}"})

        masked_fraction = float(profile_df["is_masked"].mean()) if len(profile_df) else 0.0
        if masked_fraction >= 0.10:
            severity = "fail" if masked_fraction >= 0.20 else "review"
            event_rows.append({"sample_id": sample_id, "event_type": "masked_bin_burden", "status": severity, "detail": f"masked_fraction={masked_fraction:.4f}"})

        chrom_shift = (
            profile_df.groupby("chrom", as_index=False)["z_score"]
            .median()
            .assign(abs_median_z=lambda df: df["z_score"].abs())
            .sort_values("abs_median_z", ascending=False)
        )
        if not chrom_shift.empty and float(chrom_shift.iloc[0]["abs_median_z"]) >= 2.5:
            severity = "fail" if float(chrom_shift.iloc[0]["abs_median_z"]) >= 4.0 else "review"
            event_rows.append(
                {
                    "sample_id": sample_id,
                    "event_type": "chromosome_shift_hint",
                    "status": severity,
                    "detail": f"{chrom_shift.iloc[0]['chrom']} median_abs_z={chrom_shift.iloc[0]['abs_median_z']:.4f}",
                }
            )
    return pd.DataFrame(event_rows)


def plot_counts(series, out_png, title):
    counts = series.value_counts().reindex(STATUS_ORDER, fill_value=0)
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.bar(counts.index, counts.values, color=[STATUS_COLOR[item] for item in counts.index])
    ax.set_title(title)
    ax.set_ylabel("Count")
    fig.tight_layout()
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def plot_metric_scatter(sample_df, out_png):
    fig, ax = plt.subplots(figsize=(7.0, 6.0))
    for status in STATUS_ORDER:
        subset = sample_df[sample_df["status"] == status]
        if subset.empty:
            continue
        ax.scatter(subset["bin_cv"], subset["spearman_r"], label=status, color=STATUS_COLOR[status], alpha=0.85)
    ax.set_xlabel("bin_cv")
    ax.set_ylabel("spearman_r")
    ax.set_title("Sample QC metric scatter")
    ax.legend(loc="best")
    fig.tight_layout()
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def write_outputs(tsv_path, json_path, df):
    Path(tsv_path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(tsv_path, sep="\t", index=False)
    Path(json_path).write_text(df.to_json(orient="records", indent=2), encoding="utf-8")


def main():
    args = parse_args()
    summary_df = pd.read_csv(args.baseline_summary, sep="\t")
    summary_df["status"] = summary_df["qc_decision"].map(sample_status)
    profile_map = load_profiles(args.profile_tsvs)

    combined_mask_df = pd.read_csv(args.combined_mask, sep="\t")
    annotations_df = pd.read_csv(args.bin_annotations, sep="\t")
    bin_qc_df = annotations_df.merge(
        combined_mask_df[["bin_id", "mask_label", "mask_reason"]],
        on="bin_id",
        how="left",
    )
    bin_qc_df["mask_label"] = bin_qc_df["mask_label"].fillna("pass")
    bin_qc_df["status"] = "pass"
    bin_qc_df.loc[bin_qc_df["mask_label"].isin(["soft", "dynamic"]), "status"] = "review"
    bin_qc_df.loc[bin_qc_df["mask_label"] == "hard", "status"] = "no-call"

    event_df = derive_event_rows(summary_df, profile_map, combined_mask_df)
    if event_df.empty:
        event_df = pd.DataFrame(columns=["sample_id", "event_type", "status", "detail"])

    fail_samples = set(event_df.loc[event_df["status"] == "fail", "sample_id"])
    review_samples = set(event_df.loc[event_df["status"] == "review", "sample_id"])
    nocall_samples = set(event_df.loc[event_df["status"] == "no-call", "sample_id"])
    summary_df.loc[summary_df["sample_id"].isin(nocall_samples), "status"] = "no-call"
    summary_df.loc[summary_df["sample_id"].isin(review_samples) & (summary_df["status"] == "pass"), "status"] = "review"
    summary_df.loc[summary_df["sample_id"].isin(fail_samples), "status"] = "fail"

    run_row = pd.DataFrame(
        [
            {
                "run_id": Path(args.baseline_summary).parent.name or "baseline_qc",
                "sample_count": int(len(summary_df)),
                "pass_count": int((summary_df["status"] == "pass").sum()),
                "review_count": int((summary_df["status"] == "review").sum()),
                "fail_count": int((summary_df["status"] == "fail").sum()),
                "no_call_count": int((summary_df["status"] == "no-call").sum()),
                "bin_pass_count": int((bin_qc_df["status"] == "pass").sum()),
                "bin_review_count": int((bin_qc_df["status"] == "review").sum()),
                "bin_no_call_count": int((bin_qc_df["status"] == "no-call").sum()),
                "event_count": int(len(event_df)),
                "status": "fail" if (summary_df["status"] == "fail").any() else ("review" if (summary_df["status"] == "review").any() else "pass"),
            }
        ]
    )

    write_outputs(args.run_qc_tsv, args.run_qc_json, run_row)
    write_outputs(args.sample_qc_tsv, args.sample_qc_json, summary_df)
    write_outputs(args.bin_qc_tsv, args.bin_qc_json, bin_qc_df)
    write_outputs(args.event_qc_tsv, args.event_qc_json, event_df)

    plot_counts(summary_df["status"], args.sample_status_figure, "Sample QC Status")
    plot_counts(event_df["status"], args.event_status_figure, "Event QC Status")
    plot_counts(bin_qc_df["status"], args.mask_status_figure, "Bin/Mask QC Status")
    plot_metric_scatter(summary_df, args.metric_scatter_figure)


if __name__ == "__main__":
    main()
