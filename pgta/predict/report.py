#!/biosoftware/miniconda/envs/snakemake_env/bin/python
import argparse
import json
from pathlib import Path

import pandas as pd

from pgta.core.logging import setup_logger


SEX_CHROMS = {"chrX", "chrY"}


def parse_args():
    parser = argparse.ArgumentParser(description="Build project-level CNV technical and biological candidate reports.")
    parser.add_argument("--event-tsv", action="append", default=[])
    parser.add_argument("--gender-tsv", action="append", default=[])
    parser.add_argument("--qc-tsv", action="append", default=[])
    parser.add_argument("--a-branch-bed", action="append", default=[])
    parser.add_argument("--evaluation-summary", default="")
    parser.add_argument("--ml-summary", default="")
    parser.add_argument("--benchmark-summary", default="")
    parser.add_argument("--truth-validation-summary", default="")
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--output-json", required=True)
    parser.add_argument("--output-md", required=True)
    parser.add_argument("--output-html", required=True)
    parser.add_argument("--log", default="")
    return parser.parse_args()


def ensure_parent(path_value):
    path = Path(path_value)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def load_events(paths):
    frames = []
    for path_value in paths:
        path = Path(path_value)
        if not path.exists():
            continue
        df = pd.read_csv(path, sep="\t")
        if not df.empty:
            frames.append(df)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def load_one_row_tables(paths):
    rows = []
    for path_value in paths:
        path = Path(path_value)
        if not path.exists():
            continue
        df = pd.read_csv(path, sep="\t")
        if df.empty:
            continue
        row = df.iloc[0].to_dict()
        row.setdefault("sample_id", path.stem.split(".")[0])
        rows.append(row)
    return pd.DataFrame(rows)


def load_a_branch(paths):
    rows = []
    for path_value in paths:
        path = Path(path_value)
        if not path.exists():
            continue
        df = pd.read_csv(path, sep="\t")
        sample_id = path.name.split("_")[0]
        if df.empty:
            rows.append({"sample_id": sample_id, "a_branch_event_count": 0, "a_branch_top_event": ""})
            continue
        top = df.assign(abs_z=df["zscore"].abs()).sort_values(["abs_z", "end"], ascending=[False, False]).iloc[0]
        rows.append(
            {
                "sample_id": sample_id,
                "a_branch_event_count": int(len(df)),
                "a_branch_top_event": f"chr{top['chr']}:{int(top['start'])}-{int(top['end'])} {top['type']}",
            }
        )
    return pd.DataFrame(rows)


def read_optional_json(path_value):
    if not path_value:
        return {}
    path = Path(path_value)
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def text_or_empty(value):
    if pd.isna(value):
        return ""
    return str(value)


def is_suppressed_sex_review_event(row):
    artifact_status = str(row.get("artifact_status", "") or "").strip().lower()
    chrom = str(row.get("chrom", "") or "").strip()
    keep_event = int(row.get("keep_event", 0) or 0)
    return keep_event == 1 and artifact_status == "review" and chrom in SEX_CHROMS


def summarize_branch_b_events(events_df):
    ranked_df = events_df.copy()
    ranked_df["branch_b_top_display_suppressed"] = ranked_df.apply(is_suppressed_sex_review_event, axis=1)

    sample_df = (
        ranked_df.groupby("sample_id", dropna=False)
        .agg(
            branch_b_total_events=("event_id", "size"),
            branch_b_kept_events=("keep_event", "sum"),
            branch_b_pass_events=("artifact_status", lambda values: int((values == "pass").sum())),
            branch_b_review_events=("artifact_status", lambda values: int((values == "review").sum())),
            branch_b_top_priority_score=("priority_score", "max"),
            branch_b_suppressed_sex_review_events=("branch_b_top_display_suppressed", "sum"),
        )
        .reset_index()
    )

    display_events_df = ranked_df[
        (ranked_df["keep_event"] == 1) & (~ranked_df["branch_b_top_display_suppressed"])
    ].copy()
    if display_events_df.empty:
        return sample_df, pd.DataFrame(columns=["sample_id", "branch_b_top_event"])

    top_branch_b = (
        display_events_df.sort_values(
            ["sample_id", "priority_score", "n_bins", "end"],
            ascending=[True, False, False, False],
        )
        .groupby("sample_id", dropna=False)
        .head(1)[
            [
                "sample_id",
                "chrom",
                "start",
                "end",
                "state",
                "artifact_status",
                "technical_confidence",
                *[
                    column
                    for column in [
                        "biopsy_abnormal_cell_fraction_point",
                        "biopsy_abnormal_cell_fraction_ci_low",
                        "biopsy_abnormal_cell_fraction_ci_high",
                        "biopsy_abnormal_cell_fraction_status",
                        "biopsy_abnormal_cell_fraction_reliable",
                    ]
                    if column in display_events_df.columns
                ],
                *[
                    column
                    for column in ["artifact_flags", "downgrade_reason", "filter_reason", "retain_reason"]
                    if column in display_events_df.columns
                ],
            ]
        ]
        .copy()
    )
    top_branch_b["branch_b_top_event"] = top_branch_b.apply(
        lambda row: f"{row['chrom']}:{int(row['start'])}-{int(row['end'])} {row['state']} [{row['artifact_status']}/{row['technical_confidence']}]",
        axis=1,
    )
    if "biopsy_abnormal_cell_fraction_point" in top_branch_b.columns:
        top_branch_b["branch_b_top_fraction"] = top_branch_b.apply(
            lambda row: (
                f"{100.0 * float(row['biopsy_abnormal_cell_fraction_point']):.1f}%"
                if pd.notna(row["biopsy_abnormal_cell_fraction_point"])
                else ""
            ),
            axis=1,
        )
    if {"biopsy_abnormal_cell_fraction_ci_low", "biopsy_abnormal_cell_fraction_ci_high"}.issubset(top_branch_b.columns):
        top_branch_b["branch_b_top_fraction_ci"] = top_branch_b.apply(
            lambda row: (
                f"{100.0 * float(row['biopsy_abnormal_cell_fraction_ci_low']):.1f}% to {100.0 * float(row['biopsy_abnormal_cell_fraction_ci_high']):.1f}%"
                if pd.notna(row["biopsy_abnormal_cell_fraction_ci_low"]) and pd.notna(row["biopsy_abnormal_cell_fraction_ci_high"])
                else ""
            ),
            axis=1,
        )
    top_branch_b = top_branch_b.rename(
        columns={
            "artifact_flags": "branch_b_top_flags",
            "downgrade_reason": "branch_b_top_downgrade_reason",
            "filter_reason": "branch_b_top_filter_reason",
            "retain_reason": "branch_b_top_retain_reason",
            "biopsy_abnormal_cell_fraction_status": "branch_b_top_fraction_status",
            "biopsy_abnormal_cell_fraction_reliable": "branch_b_top_fraction_reliable",
        }
    )
    top_branch_b = top_branch_b[
        [
            "sample_id",
            "branch_b_top_event",
            *[
                column
                for column in [
                    "branch_b_top_fraction",
                    "branch_b_top_fraction_ci",
                    "branch_b_top_fraction_status",
                    "branch_b_top_fraction_reliable",
                ]
                if column in top_branch_b.columns
            ],
            *[
                column
                for column in [
                    "branch_b_top_flags",
                    "branch_b_top_downgrade_reason",
                    "branch_b_top_filter_reason",
                    "branch_b_top_retain_reason",
                ]
                if column in top_branch_b.columns
            ],
        ]
    ]
    return sample_df, top_branch_b


def format_technical_conclusion(row):
    suppressed_count = int(row.get("branch_b_suppressed_sex_review_events", 0) or 0)
    top_event = text_or_empty(row.get("branch_b_top_event", "")) or ("suppressed_sex_chromosome_review" if suppressed_count > 0 else "none")
    top_fraction = text_or_empty(row.get("branch_b_top_fraction", "")) or "not_estimated"
    top_flags = text_or_empty(row.get("branch_b_top_flags", "")) or ("sex_chromosome_event" if suppressed_count > 0 else "none")
    top_downgrade = text_or_empty(row.get("branch_b_top_downgrade_reason", "")) or ("sex_chromosome_review_suppressed" if suppressed_count > 0 else "none")
    return (
        f"Branch B kept {int(row['branch_b_kept_events'])} events; "
        f"top event: {top_event}; "
        f"fraction={top_fraction}; "
        f"flags={top_flags}; "
        f"downgrade={top_downgrade}; "
        f"sex_review_suppressed={suppressed_count}; "
        f"QC={text_or_empty(row.get('qc_status', 'NA')) or 'NA'}; sex={text_or_empty(row.get('sex_call', 'NA')) or 'NA'}"
    )


def format_biological_candidate_conclusion(row):
    branch_b_top_event = text_or_empty(row.get("branch_b_top_event", ""))
    if branch_b_top_event:
        return branch_b_top_event
    a_branch_top_event = text_or_empty(row.get("a_branch_top_event", ""))
    if a_branch_top_event:
        return a_branch_top_event
    if int(row.get("branch_b_suppressed_sex_review_events", 0) or 0) > 0:
        return "Branch B top event suppressed (sex-chromosome review only)"
    return "No A-branch event"


def main():
    args = parse_args()
    logger = setup_logger("cnv_report", args.log or None)
    events_df = load_events(args.event_tsv)
    gender_df = load_one_row_tables(args.gender_tsv)
    qc_df = load_one_row_tables(args.qc_tsv)
    a_branch_df = load_a_branch(args.a_branch_bed)
    evaluation_summary = read_optional_json(args.evaluation_summary)
    ml_summary = read_optional_json(args.ml_summary)
    benchmark_summary = read_optional_json(args.benchmark_summary)
    truth_validation_summary = read_optional_json(args.truth_validation_summary)
    fraction_benchmark = benchmark_summary.get("fraction_estimation", {}) if isinstance(benchmark_summary, dict) else {}
    low_fraction_detection = benchmark_summary.get("low_fraction_detection", []) if isinstance(benchmark_summary, dict) else []

    if events_df.empty:
        empty = pd.DataFrame(columns=["sample_id"])
        ensure_parent(args.output_tsv)
        empty.to_csv(args.output_tsv, sep="\t", index=False)
        ensure_parent(args.output_json).write_text(json.dumps({"status": "empty"}, indent=2), encoding="utf-8")
        ensure_parent(args.output_md).write_text("# CNV Report\n\nNo events available.\n", encoding="utf-8")
        ensure_parent(args.output_html).write_text("<html><body><h1>CNV Report</h1><p>No events available.</p></body></html>", encoding="utf-8")
        logger.info("report skipped: no events")
        return

    sample_df, top_branch_b = summarize_branch_b_events(events_df)

    sample_df = sample_df.merge(top_branch_b, on="sample_id", how="left")
    if not gender_df.empty:
        sample_df = sample_df.merge(
            gender_df[[column for column in ["sample_id", "sex_call", "sex_call_source", "predict_gender"] if column in gender_df.columns]],
            on="sample_id",
            how="left",
        )
    if not qc_df.empty:
        sample_df = sample_df.merge(
            qc_df[[column for column in ["sample_id", "status", "mad_log1p", "nonzero_fraction"] if column in qc_df.columns]].rename(
                columns={"status": "qc_status"}
            ),
            on="sample_id",
            how="left",
        )
    if not a_branch_df.empty:
        sample_df = sample_df.merge(a_branch_df, on="sample_id", how="left")

    sample_df["technical_conclusion"] = sample_df.apply(format_technical_conclusion, axis=1)
    sample_df["biological_candidate_conclusion"] = sample_df.apply(format_biological_candidate_conclusion, axis=1)

    ensure_parent(args.output_tsv)
    sample_df.to_csv(args.output_tsv, sep="\t", index=False)
    payload = {
        "status": "completed",
        "evaluation_summary": evaluation_summary,
        "ml_summary": ml_summary,
        "benchmark_summary": benchmark_summary,
        "samples": sample_df.to_dict(orient="records"),
    }
    if truth_validation_summary:
        payload["mosaic_truth_validation_summary"] = truth_validation_summary
    ensure_parent(args.output_json).write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")

    md_lines = [
        "# CNV Report",
        "",
        "## Project Summary",
        f"- Samples: `{sample_df['sample_id'].nunique()}`",
        f"- Branch B kept events: `{int(sample_df['branch_b_kept_events'].sum())}`",
        f"- Evaluation status: `{evaluation_summary.get('status', 'not_run')}`",
        f"- ML status: `{ml_summary.get('status', 'not_run')}`",
        f"- Benchmark status: `{benchmark_summary.get('status', 'not_run')}`",
        "",
        "## Sample Table",
        "",
        "| Sample | QC | Sex | Branch B Top Event | Technical Conclusion | Biological Candidate Conclusion |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    if fraction_benchmark:
        md_lines[8:8] = [
            f"- Fraction evaluable events: `{fraction_benchmark.get('evaluable_event_count', 0)}`",
            f"- Fraction MAE: `{fraction_benchmark.get('mae', 'NA')}`",
            f"- Fraction RMSE: `{fraction_benchmark.get('rmse', 'NA')}`",
            f"- Fraction CI coverage: `{fraction_benchmark.get('ci_coverage', 'NA')}`",
        ]
    if truth_validation_summary:
        md_lines[8:8] = [
            f"- Mosaic truth validation status: `{truth_validation_summary.get('status', 'unknown')}`",
            f"- Mosaic truth rows: `{truth_validation_summary.get('row_count', 'NA')}`",
            f"- Mosaic truth rows with filled fraction: `{truth_validation_summary.get('complete_fraction_row_count', 'NA')}`",
            f"- Mosaic truth rows missing fraction: `{truth_validation_summary.get('missing_fraction_row_count', 'NA')}`",
        ]
    if low_fraction_detection:
        md_lines.insert(
            8,
            "- Low-fraction detection: `"
            + "; ".join(
                [
                    f"<= {100.0 * float(item.get('fraction_threshold', 0.0)):.0f}%: B={item.get('branch_b_detection_rate', 'NA')}, A={item.get('a_branch_detection_rate', 'NA')}"
                    for item in low_fraction_detection
                ]
            )
            + "`",
        )
    for row in sample_df.itertuples(index=False):
        md_lines.append(
            f"| `{row.sample_id}` | `{getattr(row, 'qc_status', 'NA')}` | `{getattr(row, 'sex_call', 'NA')}` | "
            f"{getattr(row, 'branch_b_top_event', '') or 'none'} | {row.technical_conclusion} | {row.biological_candidate_conclusion} |"
        )
    ensure_parent(args.output_md).write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    html_rows = []
    for row in sample_df.itertuples(index=False):
        html_rows.append(
            "<tr>"
            f"<td>{row.sample_id}</td>"
            f"<td>{getattr(row, 'qc_status', 'NA')}</td>"
            f"<td>{getattr(row, 'sex_call', 'NA')}</td>"
            f"<td>{getattr(row, 'branch_b_top_event', '') or 'none'}</td>"
            f"<td>{row.technical_conclusion}</td>"
            f"<td>{row.biological_candidate_conclusion}</td>"
            "</tr>"
        )
    html = (
        "<html><head><meta charset='utf-8'><title>CNV Report</title>"
        "<style>body{font-family:Arial,sans-serif;margin:24px;}table{border-collapse:collapse;width:100%;}"
        "th,td{border:1px solid #ccc;padding:8px;vertical-align:top;}th{background:#f3f3f3;text-align:left;}</style>"
        "</head><body><h1>CNV Report</h1>"
        f"<p>Samples: {sample_df['sample_id'].nunique()} | Branch B kept events: {int(sample_df['branch_b_kept_events'].sum())}</p>"
        "<table><thead><tr><th>Sample</th><th>QC</th><th>Sex</th><th>Branch B Top Event</th>"
        "<th>Technical Conclusion</th><th>Biological Candidate Conclusion</th></tr></thead><tbody>"
        + "".join(html_rows)
        + "</tbody></table></body></html>"
    )
    ensure_parent(args.output_html).write_text(html, encoding="utf-8")
    logger.info("report completed: samples=%d", sample_df["sample_id"].nunique())


if __name__ == "__main__":
    main()
