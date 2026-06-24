#!/biosoftware/miniconda/envs/snakemake_env/bin/python
import argparse
import html as html_lib
import json
import os
from pathlib import Path

import pandas as pd

from pgta.core.logging import setup_logger


SEX_CHROMS = {"chrX", "chrY"}
A_BRANCH_REVIEW_SHORTLIST_SIZE = 3
A_BRANCH_STRONG_SIGNAL_Z = 10.0
A_BRANCH_INTERNAL_COLUMNS = {
    "a_branch_review_candidate_count",
    "a_branch_review_shortlist",
    "a_branch_strong_signal_count",
}
CNVSEQ_TIER_RANK = {
    "whole_chromosome": 0,
    "reportable": 1,
    "review_1_2mb": 2,
    "subreportable_lt1mb": 3,
}
BRANCH_B_SAMPLE_DEFAULTS = {
    "branch_b_total_events": 0,
    "branch_b_kept_events": 0,
    "branch_b_pass_events": 0,
    "branch_b_review_events": 0,
    "branch_b_reportable_events": 0,
    "branch_b_review_tier_events": 0,
    "branch_b_subreportable_events": 0,
    "branch_b_top_priority_score": 0.0,
    "branch_b_suppressed_sex_review_events": 0,
}


def parse_args():
    parser = argparse.ArgumentParser(description="Build project-level CNV technical and biological candidate reports.")
    parser.add_argument("--event-tsv", action="append", default=[])
    parser.add_argument("--gender-tsv", action="append", default=[])
    parser.add_argument("--qc-tsv", action="append", default=[])
    parser.add_argument("--a-branch-bed", action="append", default=[])
    parser.add_argument("--branch-a-validation-summary", action="append", default=[])
    parser.add_argument("--branch-b-evidence-summary", action="append", default=[])
    parser.add_argument("--branch-s-summary", action="append", default=[])
    parser.add_argument("--branch-b-v2-benchmark-summary", default="")
    parser.add_argument("--branch-b-v2-sample-summary", default="")
    parser.add_argument("--event-annotation-tsv", default="")
    parser.add_argument("--predict-bam-qc-summary", default="")
    parser.add_argument("--sample-report-qc", default="")
    parser.add_argument("--reference-id", default="")
    parser.add_argument("--wisecondorx-predict-command", default="")
    parser.add_argument("--evaluation-summary", default="")
    parser.add_argument("--ml-summary", default="")
    parser.add_argument("--benchmark-summary", default="")
    parser.add_argument("--truth-validation-summary", default="")
    parser.add_argument("--plot-svg", action="append", default=[])
    parser.add_argument("--copy-number-plot-svg", action="append", default=[])
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


def apply_event_annotations(events_df, annotation_tsv=""):
    if events_df is None or events_df.empty or not annotation_tsv:
        return events_df
    path = Path(annotation_tsv)
    if not path.exists() or path.stat().st_size == 0:
        return events_df
    annotations = pd.read_csv(path, sep="\t")
    if annotations.empty:
        return events_df
    annotation_columns = [
        "cytoband",
        "genes",
        "gene_number",
        "gene_location",
        "omim_genes",
        "omim_phenotypes",
        "hpo_terms",
        "region_context",
        "gene_source_status",
        "omim_source_status",
        "hpo_source_status",
        "annotation_backend",
        "annotation_status",
        "annotation_bundle_id",
        "genome_build",
    ]
    frame = events_df.copy()
    if "event_id" in frame.columns and "event_id" in annotations.columns:
        keys = ["sample_id", "event_id"] if "sample_id" in frame.columns and "sample_id" in annotations.columns else ["event_id"]
    else:
        keys = [
            column
            for column in ["sample_id", "chrom", "start", "end", "state"]
            if column in frame.columns and column in annotations.columns
        ]
    if not keys:
        return frame
    merge_columns = keys + [column for column in annotation_columns if column in annotations.columns]
    right = annotations[merge_columns].drop_duplicates(subset=keys, keep="first")
    merged = frame.merge(right, on=keys, how="left")
    if len(merged) != len(frame):
        raise ValueError("event annotation merge changed event row count")
    for column in annotation_columns:
        if column not in merged.columns:
            merged[column] = ""
    return merged


def summarize_event_annotation_status(annotation_tsv=""):
    if not annotation_tsv:
        return {
            "event_annotation_status": "not_configured",
            "event_annotation_row_count": 0,
            "event_annotation_backend": "not_configured",
            "event_annotation_bundle_id": "",
            "event_annotation_gene_source_status": "not_configured",
            "event_annotation_omim_source_status": "not_configured",
            "event_annotation_hpo_source_status": "not_configured",
        }
    path = Path(annotation_tsv)
    if not path.exists() or path.stat().st_size == 0:
        return {
            "event_annotation_status": "missing",
            "event_annotation_row_count": 0,
            "event_annotation_backend": "pgta_sqlite",
            "event_annotation_bundle_id": "",
            "event_annotation_gene_source_status": "missing",
            "event_annotation_omim_source_status": "missing",
            "event_annotation_hpo_source_status": "missing",
        }
    annotations = pd.read_csv(path, sep="\t")
    if annotations.empty:
        return {
            "event_annotation_status": "empty",
            "event_annotation_row_count": 0,
            "event_annotation_backend": "pgta_sqlite",
            "event_annotation_bundle_id": "",
            "event_annotation_gene_source_status": "empty",
            "event_annotation_omim_source_status": "empty",
            "event_annotation_hpo_source_status": "empty",
        }
    statuses = annotations.get("annotation_status", pd.Series("", index=annotations.index)).fillna("").astype(str)
    bundle_ids = annotations.get("annotation_bundle_id", pd.Series("", index=annotations.index)).dropna().astype(str)
    backends = annotations.get("annotation_backend", pd.Series("", index=annotations.index)).dropna().astype(str)
    gene_source_statuses = annotations.get("gene_source_status", pd.Series("", index=annotations.index)).dropna().astype(str)
    omim_source_statuses = annotations.get("omim_source_status", pd.Series("", index=annotations.index)).dropna().astype(str)
    hpo_source_statuses = annotations.get("hpo_source_status", pd.Series("", index=annotations.index)).dropna().astype(str)
    return {
        "event_annotation_status": ",".join(sorted(set(statuses))) if not statuses.empty else "unknown",
        "event_annotation_row_count": int(len(annotations)),
        "event_annotation_backend": ",".join(sorted(set(backends))) if not backends.empty else "pgta_sqlite",
        "event_annotation_bundle_id": ",".join(sorted(set(bundle_ids))) if not bundle_ids.empty else "",
        "event_annotation_gene_source_status": ",".join(sorted(set(gene_source_statuses))) if not gene_source_statuses.empty else "unknown",
        "event_annotation_omim_source_status": ",".join(sorted(set(omim_source_statuses))) if not omim_source_statuses.empty else "unknown",
        "event_annotation_hpo_source_status": ",".join(sorted(set(hpo_source_statuses))) if not hpo_source_statuses.empty else "unknown",
    }


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


def load_sample_report_qc(path_value=""):
    columns = [
        "sample_id",
        "sample_report_qc_status",
        "sample_report_qc_reasons",
        "rebuild_or_resample_recommendation",
        "library_qc_status",
        "global_wave_context",
        "multi_chromosome_cnv_context",
        "possible_contamination_or_mixture_context",
    ]
    if not path_value:
        return pd.DataFrame(columns=columns)
    path = Path(path_value)
    if not path.exists():
        return pd.DataFrame(columns=columns)
    df = pd.read_csv(path, sep="\t")
    for column in columns:
        if column not in df.columns:
            df[column] = ""
    return df[columns]


def summarize_sample_report_qc(df):
    if df.empty or "sample_report_qc_status" not in df.columns:
        return {
            "sample_report_qc_status": "not_configured",
            "sample_report_qc_sample_count": 0,
            "pass_reportable_count": 0,
            "sample_quality_review_count": 0,
            "no_call_recommended_count": 0,
        }
    status = df["sample_report_qc_status"].astype(str)
    return {
        "sample_report_qc_status": "completed",
        "sample_report_qc_sample_count": int(df["sample_id"].nunique()) if "sample_id" in df.columns else int(len(df)),
        "pass_reportable_count": int((status == "PASS_REPORTABLE").sum()),
        "sample_quality_review_count": int((status == "SAMPLE_QUALITY_REVIEW").sum()),
        "no_call_recommended_count": int((status == "NO_CALL_RECOMMENDED").sum()),
    }


def format_sample_report_qc_status(row):
    status = text_or_empty(row.get("sample_report_qc_status", ""))
    if not status:
        return "not_available"
    reasons = text_or_empty(row.get("sample_report_qc_reasons", ""))
    recommendation = text_or_empty(row.get("rebuild_or_resample_recommendation", ""))
    library = text_or_empty(row.get("library_qc_status", ""))
    details = [status]
    if library:
        details.append(f"library={library}")
    if reasons and reasons != "PASS":
        details.append(f"reasons={reasons}")
    if recommendation:
        details.append(f"recommendation={recommendation}")
    return "; ".join(details)


def format_a_branch_event(row, include_z=False):
    chrom = str(row["chr"])
    if not chrom.startswith("chr"):
        chrom = f"chr{chrom}"
    event = f"{chrom}:{int(row['start'])}-{int(row['end'])} {row['type']}"
    if include_z and pd.notna(row.get("abs_zscore")):
        event = f"{event} z={float(row['abs_zscore']):.2f}"
    return event


def merge_a_branch_events(df):
    required_columns = {"chr", "start", "end", "type", "zscore"}
    if df.empty or not required_columns.issubset(df.columns):
        return pd.DataFrame(columns=["chr", "start", "end", "type", "abs_zscore", "source_event_count"])

    normalized = df.copy()
    normalized["start"] = pd.to_numeric(normalized["start"], errors="coerce")
    normalized["end"] = pd.to_numeric(normalized["end"], errors="coerce")
    normalized["zscore"] = pd.to_numeric(normalized["zscore"], errors="coerce")
    normalized = normalized.dropna(subset=["start", "end", "zscore", "chr", "type"]).copy()
    if normalized.empty:
        return pd.DataFrame(columns=["chr", "start", "end", "type", "abs_zscore", "source_event_count"])

    normalized["abs_zscore"] = normalized["zscore"].abs()
    normalized = normalized.sort_values(["chr", "type", "start", "end"])
    merged_rows = []
    current = None
    for row in normalized.itertuples(index=False):
        row_chr = str(getattr(row, "chr"))
        row_type = str(getattr(row, "type"))
        row_start = int(getattr(row, "start"))
        row_end = int(getattr(row, "end"))
        row_abs_z = float(getattr(row, "abs_zscore"))
        if (
            current is None
            or current["chr"] != row_chr
            or current["type"] != row_type
            or row_start > current["end"] + 1
        ):
            if current is not None:
                merged_rows.append(current)
            current = {
                "chr": row_chr,
                "start": row_start,
                "end": row_end,
                "type": row_type,
                "abs_zscore": row_abs_z,
                "source_event_count": 1,
            }
            continue
        current["end"] = max(current["end"], row_end)
        current["abs_zscore"] = max(current["abs_zscore"], row_abs_z)
        current["source_event_count"] += 1
    if current is not None:
        merged_rows.append(current)

    merged_df = pd.DataFrame(merged_rows)
    if merged_df.empty:
        return merged_df
    merged_df["span"] = merged_df["end"] - merged_df["start"] + 1
    return merged_df.sort_values(["abs_zscore", "span"], ascending=[False, False]).drop(columns=["span"])


def load_a_branch(paths):
    rows = []
    for path_value in paths:
        path = Path(path_value)
        if not path.exists():
            continue
        df = pd.read_csv(path, sep="\t")
        sample_id = path.name.split("_")[0]
        if df.empty:
            rows.append(
                {
                    "sample_id": sample_id,
                    "a_branch_event_count": 0,
                    "a_branch_top_event": "",
                    "a_branch_review_candidate_count": 0,
                    "a_branch_review_shortlist": "",
                    "a_branch_strong_signal_count": 0,
                }
            )
            continue
        merged_df = merge_a_branch_events(df)
        if merged_df.empty:
            rows.append(
                {
                    "sample_id": sample_id,
                    "a_branch_event_count": 0,
                    "a_branch_top_event": "",
                    "a_branch_review_candidate_count": 0,
                    "a_branch_review_shortlist": "",
                    "a_branch_strong_signal_count": 0,
                }
            )
            continue
        shortlist = merged_df.head(A_BRANCH_REVIEW_SHORTLIST_SIZE).copy()
        top = shortlist.iloc[0]
        rows.append(
            {
                "sample_id": sample_id,
                "a_branch_event_count": int(len(merged_df)),
                "a_branch_top_event": format_a_branch_event(top),
                "a_branch_review_candidate_count": int(len(shortlist)),
                "a_branch_review_shortlist": "; ".join(
                    shortlist.apply(lambda row: format_a_branch_event(row, include_z=True), axis=1).tolist()
                ),
                "a_branch_strong_signal_count": int(
                    (merged_df["abs_zscore"].astype(float) >= A_BRANCH_STRONG_SIGNAL_Z).sum()
                ),
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


def _sample_from_summary_path(path):
    name = path.name
    for suffix in [".summary.json", ".json"]:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def _format_count_dict(value):
    if not isinstance(value, dict) or not value:
        return "not_recorded"
    return ";".join(f"{key}={int(count)}" for key, count in sorted(value.items()))


def load_branch_b_evidence_summaries(paths):
    rows = []
    for path_value in paths:
        path = Path(path_value)
        if not path.exists():
            continue
        summary = read_optional_json(path)
        if not summary:
            continue
        p3_counts = summary.get("p3_disposition_counts", {})
        review_required_count = int(p3_counts.get("REVIEW_REQUIRED", summary.get("review_burden_count", 0)) or 0)
        rows.append(
            {
                "sample_id": str(summary.get("sample_id") or _sample_from_summary_path(path)),
                "branch_b_evidence_candidate_count": int(summary.get("candidate_count", 0) or 0),
                "branch_b_evidence_review_required_count": review_required_count,
                "branch_b_evidence_missing_count": int(summary.get("missing_evidence_candidate_count", 0) or 0),
                "branch_b_evidence_background_source": _format_count_dict(summary.get("background_source_counts", {})),
                "branch_b_evidence_background_status": _format_count_dict(summary.get("background_status_counts", {})),
                "branch_b_evidence_disposition_counts": _format_count_dict(p3_counts),
                "branch_b_evidence_final_report_impact": str(summary.get("final_report_impact", "not_recorded")),
            }
        )
    return pd.DataFrame(rows)


def load_branch_b_v2_burden_summaries(sample_summary_path="", benchmark_summary_path="", report_reference_id=""):
    summary = read_optional_json(benchmark_summary_path)
    if not summary:
        return pd.DataFrame(), {
            "branch_b_v2_burden_status": "not_configured",
            "branch_b_v2_background_unknown_review_count": 0,
            "branch_b_v2_branch_s_review_count": 0,
            "branch_b_v2_technical_risk_review_count": 0,
            "branch_b_v2_report_candidate_count": 0,
            "branch_b_v2_report_event_count": 0,
            "branch_b_v2_report_strong_event_count": 0,
            "branch_b_v2_report_weak_event_count": 0,
            "branch_b_v2_internal_review_event_count": 0,
            "branch_b_v2_filtered_event_count": 0,
            "branch_b_v2_branch_s_event_count": 0,
            "branch_b_v2_legacy_fields_used": False,
            "branch_b_v2_final_impact": "development_review_only",
            "branch_b_v2_same_reference_config_status": "not_configured",
            "branch_b_v2_sample_summary_count": 0,
            "branch_b_v2_sample_report_burden_flag_count": 0,
            "branch_b_v2_sample_report_burden_threshold": 0,
            "branch_b_v2_note": "Branch B V2 burden stratification is not configured for this report.",
        }

    summary_reference = str(summary.get("reference_id", ""))
    report_reference = str(report_reference_id or "")
    same_reference = bool(summary_reference and report_reference and summary_reference == report_reference)
    same_reference_status = "matched" if same_reference else "mismatch"
    legacy_used = bool(summary.get("legacy_branch_b_decision_fields_used", False))
    contract = {
        "branch_b_v2_burden_status": str(summary.get("status", "unknown")),
        "branch_b_v2_background_unknown_review_count": int(
            summary.get("v2_background_unknown_review_burden_count", 0) or 0
        ),
        "branch_b_v2_branch_s_review_count": int(summary.get("v2_branch_s_review_burden_count", 0) or 0),
        "branch_b_v2_technical_risk_review_count": int(summary.get("v2_technical_risk_burden_count", 0) or 0),
        "branch_b_v2_report_candidate_count": int(summary.get("v2_report_candidate_burden_count", 0) or 0),
        "branch_b_v2_report_event_count": int(summary.get("v2_report_event_count", 0) or 0),
        "branch_b_v2_report_strong_event_count": int(summary.get("v2_report_strong_event_count", 0) or 0),
        "branch_b_v2_report_weak_event_count": int(summary.get("v2_report_weak_event_count", 0) or 0),
        "branch_b_v2_internal_review_event_count": int(summary.get("v2_internal_review_event_count", 0) or 0),
        "branch_b_v2_filtered_event_count": int(summary.get("v2_filtered_event_count", 0) or 0),
        "branch_b_v2_branch_s_event_count": int(summary.get("v2_branch_s_event_count", 0) or 0),
        "branch_b_v2_legacy_fields_used": legacy_used,
        "branch_b_v2_final_impact": "development_review_only",
        "branch_b_v2_same_reference_config_status": same_reference_status,
        "branch_b_v2_sample_summary_count": int(summary.get("sample_count", 0) or 0),
        "branch_b_v2_sample_report_burden_flag_count": int(
            summary.get("sample_report_burden_flag_count", 0) or 0
        ),
        "branch_b_v2_sample_report_burden_threshold": int(
            summary.get("sample_report_burden_threshold", 0) or 0
        ),
        "branch_b_v2_note": (
            "Branch B V2 burden stratification is displayed as development_review_only evidence; "
            "it is not FP-reduction proof and does not promote Branch B V2 or Branch S."
        ),
    }

    sample_path = Path(sample_summary_path) if sample_summary_path else None
    if sample_path is None or not sample_path.exists():
        return pd.DataFrame(), contract

    sample_summary = pd.read_csv(sample_path, sep="\t")
    if sample_summary.empty or "sample_id" not in sample_summary.columns:
        return pd.DataFrame(), contract

    rows = []
    for row in sample_summary.to_dict(orient="records"):
        rows.append(
            {
                "sample_id": str(row.get("sample_id", "")),
                "branch_b_v2_burden_status": contract["branch_b_v2_burden_status"],
                "branch_b_v2_background_unknown_review_count": int(
                    row.get("v2_background_unknown_review_burden_count", 0) or 0
                ),
                "branch_b_v2_branch_s_review_count": int(row.get("v2_branch_s_review_burden_count", 0) or 0),
                "branch_b_v2_technical_risk_review_count": int(row.get("v2_technical_risk_burden_count", 0) or 0),
                "branch_b_v2_report_candidate_count": int(row.get("v2_report_candidate_burden_count", 0) or 0),
                "branch_b_v2_review_candidate_count": int(row.get("v2_review_candidate_burden_count", 0) or 0),
                "branch_b_v2_report_event_count": int(row.get("v2_report_event_count", 0) or 0),
                "branch_b_v2_report_strong_event_count": int(row.get("v2_report_strong_event_count", 0) or 0),
                "branch_b_v2_report_weak_event_count": int(row.get("v2_report_weak_event_count", 0) or 0),
                "branch_b_v2_internal_review_event_count": int(row.get("v2_internal_review_event_count", 0) or 0),
                "branch_b_v2_filtered_event_count": int(row.get("v2_filtered_event_count", 0) or 0),
                "branch_b_v2_branch_s_event_count": int(row.get("v2_branch_s_event_count", 0) or 0),
                "branch_b_v2_sample_report_burden_flag": int(row.get("sample_report_burden_flag", 0) or 0),
                "branch_b_v2_sample_report_burden_reason": str(
                    row.get("sample_report_burden_reason", "") or ""
                ),
                "branch_b_v2_legacy_fields_used": legacy_used,
                "branch_b_v2_final_impact": "development_review_only",
                "branch_b_v2_same_reference_config_status": same_reference_status,
            }
        )
    return pd.DataFrame(rows), contract


def load_branch_a_validation_summaries(paths):
    rows = []
    for path_value in paths:
        path = Path(path_value)
        if not path.exists():
            continue
        summary = read_optional_json(path)
        if not summary:
            continue
        rows.append(
            {
                "path": str(path),
                "reference_id": str(summary.get("reference_id", "")),
                "sample_count": int(summary.get("sample_count", 0) or 0),
                "truth_event_count": int(summary.get("truth_event_count", 0) or 0),
                "truth_detected_count": int(summary.get("truth_detected_count", 0) or 0),
                "FN_count": int(summary.get("FN_count", 0) or 0),
                "H6_chr21_status": str(summary.get("H6_chr21_status", "not_applicable")),
                "status": str(summary.get("status", "unknown")),
            }
        )
    return rows


def summarize_branch_a_validation_gate(reference_id, summaries):
    if not summaries:
        return {
            "status": "not_configured",
            "summary_count": 0,
            "truth_event_count": 0,
            "truth_detected_count": 0,
            "FN_count": 0,
            "H6_chr21_status": "not_configured",
            "same_reference_config_status": "not_configured",
            "summaries": [],
        }
    reference_values = {str(item.get("reference_id", "")) for item in summaries}
    same_reference = reference_values == {str(reference_id)}
    fn_count = int(sum(int(item.get("FN_count", 0) or 0) for item in summaries))
    truth_event_count = int(sum(int(item.get("truth_event_count", 0) or 0) for item in summaries))
    truth_detected_count = int(sum(int(item.get("truth_detected_count", 0) or 0) for item in summaries))
    h6_statuses = [str(item.get("H6_chr21_status", "")) for item in summaries]
    if "detected" in h6_statuses:
        h6_status = "detected"
    elif h6_statuses:
        h6_status = ";".join(sorted(set(status for status in h6_statuses if status)))
    else:
        h6_status = "not_recorded"
    if not same_reference:
        status = "reference_mismatch"
        same_reference_status = "mismatch"
    elif truth_event_count <= 0:
        status = "no_truth_events"
        same_reference_status = "matched"
    elif fn_count == 0 and truth_detected_count >= truth_event_count:
        status = "passed_no_fn"
        same_reference_status = "matched"
    else:
        status = "failed_or_incomplete"
        same_reference_status = "matched"
    return {
        "status": status,
        "summary_count": int(len(summaries)),
        "truth_event_count": truth_event_count,
        "truth_detected_count": truth_detected_count,
        "FN_count": fn_count,
        "H6_chr21_status": h6_status,
        "same_reference_config_status": same_reference_status,
        "summaries": summaries,
    }


def load_branch_s_summaries(paths):
    rows = []
    for path_value in paths:
        path = Path(path_value)
        if not path.exists():
            continue
        summary = read_optional_json(path)
        if not summary:
            continue
        rows.append(
            {
                "sample_id": str(summary.get("sample_id") or summary.get("sample") or _sample_from_summary_path(path)),
                "sca_candidate_state": str(summary.get("sca_candidate_state", "not_available")),
                "sca_confidence_tier": str(summary.get("sca_confidence_tier", "not_available")),
                "sca_report_layer_class": str(summary.get("sca_report_layer_class", "not_available")),
                "sca_report_layer_reason": str(summary.get("sca_report_layer_reason", "")),
                "sca_output_mode": str(summary.get("sca_output_mode", "not_available")),
                "sca_uncertainty_reason": str(summary.get("sca_uncertainty_reason", "")),
                "report_text_status": str(summary.get("report_text_status", "not_available")),
            }
        )
    return pd.DataFrame(rows)


def text_or_empty(value):
    if pd.isna(value):
        return ""
    return str(value)


def build_plot_lookup(paths):
    lookup = {}
    for path_value in paths:
        if not path_value:
            continue
        path = Path(path_value)
        name = path.name
        if name.endswith(".final_cnv_cn.svg"):
            sample_id = name[: -len(".final_cnv_cn.svg")]
        elif name.endswith(".final_cnv.svg"):
            sample_id = name[: -len(".final_cnv.svg")]
        elif name.endswith(".branch_ab.svg"):
            sample_id = name[: -len(".branch_ab.svg")]
        else:
            sample_id = path.stem
        if sample_id:
            lookup[sample_id] = str(path)
    return lookup


def format_plot_link(sample_id, output_path, plot_lookup, html=False, label="plot"):
    plot_path = plot_lookup.get(str(sample_id), "")
    if not plot_path:
        return ""
    relative = os.path.relpath(plot_path, start=str(Path(output_path).parent))
    if html:
        escaped_href = html_lib.escape(relative)
        return f"<a href='{escaped_href}'>{html_lib.escape(label)}</a>"
    return f"[{label}]({relative})"


def is_suppressed_sex_review_event(row):
    artifact_status = str(row.get("artifact_status", "") or "").strip().lower()
    chrom = str(row.get("chrom", "") or "").strip()
    keep_event = int(row.get("keep_event", 0) or 0)
    return keep_event == 1 and artifact_status == "review" and chrom in SEX_CHROMS


def cnvseq_tier_rank(value):
    return CNVSEQ_TIER_RANK.get(str(value or "").strip(), 9)


def format_branch_b_top_event(row):
    tier = str(row.get("cnvseq_report_tier", "") or "").strip()
    suffix_items = [tier] if tier else []
    cytoband = text_or_empty(row.get("cytoband", ""))
    genes = text_or_empty(row.get("genes", ""))
    if cytoband:
        suffix_items.append(cytoband)
    if genes:
        suffix_items.append(genes)
    suffix = f"; {'; '.join(suffix_items)}" if suffix_items else ""
    return (
        f"{row['chrom']}:{int(row['start'])}-{int(row['end'])} {row['state']} "
        f"[{row['artifact_status']}/{row['technical_confidence']}{suffix}]"
    )


def format_branch_b_evidence_status(row):
    count = row.get("branch_b_evidence_candidate_count", "")
    if pd.isna(count) or count == "":
        return "not_available"
    background = text_or_empty(row.get("branch_b_evidence_background_status", "")) or "not_recorded"
    impact = text_or_empty(row.get("branch_b_evidence_final_report_impact", "")) or "not_recorded"
    review = int(row.get("branch_b_evidence_review_required_count", 0) or 0)
    return f"candidates={int(count)}; review_required={review}; background={background}; impact={impact}"


def format_branch_b_v2_burden_status(row):
    status = row.get("branch_b_v2_burden_status", "")
    if pd.isna(status) or status == "":
        return "not_available"
    return (
        f"status={status}; "
        f"background_unknown_review={int(row.get('branch_b_v2_background_unknown_review_count', 0) or 0)}; "
        f"branch_s_review={int(row.get('branch_b_v2_branch_s_review_count', 0) or 0)}; "
        f"technical_risk_review={int(row.get('branch_b_v2_technical_risk_review_count', 0) or 0)}; "
        f"report_candidate={int(row.get('branch_b_v2_report_candidate_count', 0) or 0)}; "
        f"report_event={int(row.get('branch_b_v2_report_event_count', 0) or 0)}; "
        f"report_strong={int(row.get('branch_b_v2_report_strong_event_count', 0) or 0)}; "
        f"report_weak={int(row.get('branch_b_v2_report_weak_event_count', 0) or 0)}; "
        f"internal_review_event={int(row.get('branch_b_v2_internal_review_event_count', 0) or 0)}; "
        f"filtered_event={int(row.get('branch_b_v2_filtered_event_count', 0) or 0)}; "
        f"branch_s_event={int(row.get('branch_b_v2_branch_s_event_count', 0) or 0)}; "
        f"sample_report_burden_flag={int(row.get('branch_b_v2_sample_report_burden_flag', 0) or 0)}; "
        f"impact={text_or_empty(row.get('branch_b_v2_final_impact', 'development_review_only')) or 'development_review_only'}"
    )


def cnvseq_tier_from_length_tier(value):
    text = str(value or "").strip()
    if text == "large_ge10mb":
        return "whole_chromosome"
    if text in {"broad_review_ge4mb", "reportable_candidate_ge2mb"}:
        return "reportable"
    if text == "review_only_1_2mb":
        return "review_1_2mb"
    if text == "subreportable_lt1mb":
        return "subreportable_lt1mb"
    return ""


def normalize_branch_b_report_events(events_df):
    ranked_df = events_df.copy()
    has_v2_report_contract = "v2_report_layer_class" in ranked_df.columns or "v2_report_visibility" in ranked_df.columns
    if not has_v2_report_contract:
        return ranked_df

    report_class = ranked_df.get("v2_report_layer_class", pd.Series("", index=ranked_df.index)).fillna("").astype(str)
    visibility = ranked_df.get("v2_report_visibility", pd.Series("", index=ranked_df.index)).fillna("").astype(str)
    report_mask = report_class.eq("report_event") | visibility.isin(
        {"report_strong_event", "report_weak_event", "final_report"}
    )
    ranked_df = ranked_df.loc[report_mask].copy()

    if "event_id" not in ranked_df.columns:
        ranked_df["event_id"] = ranked_df.get(
            "candidate_id",
            pd.Series([f"v2_report_event_{idx}" for idx in ranked_df.index], index=ranked_df.index),
        )
    if "keep_event" not in ranked_df.columns:
        ranked_df["keep_event"] = 1
    if "artifact_status" not in ranked_df.columns:
        ranked_df["artifact_status"] = visibility.map(
            {
                "report_strong_event": "pass",
                "report_weak_event": "review",
                "final_report": "review",
            }
        ).fillna("review")
    if "technical_confidence" not in ranked_df.columns:
        ranked_df["technical_confidence"] = visibility.map(
            {
                "report_strong_event": "strong",
                "report_weak_event": "weak",
                "final_report": "report",
            }
        ).fillna("report")
    if "priority_score" not in ranked_df.columns:
        ranked_df["priority_score"] = ranked_df.get(
            "a_abs_zscore",
            ranked_df.get("abs_zscore", pd.Series(0.0, index=ranked_df.index)),
        )
    if "n_bins" not in ranked_df.columns:
        ranked_df["n_bins"] = ranked_df.get("event_bin_count", pd.Series(0, index=ranked_df.index))
    if "cnvseq_report_tier" not in ranked_df.columns:
        ranked_df["cnvseq_report_tier"] = ranked_df.get(
            "v2_length_tier",
            pd.Series("", index=ranked_df.index),
        ).map(cnvseq_tier_from_length_tier)
    if "cnvseq_reportable" not in ranked_df.columns:
        ranked_df["cnvseq_reportable"] = ranked_df["cnvseq_report_tier"].isin(
            {"whole_chromosome", "reportable"}
        ).astype(int)
    for column in ("artifact_flags", "downgrade_reason", "filter_reason", "retain_reason"):
        if column not in ranked_df.columns:
            ranked_df[column] = ""
    return ranked_df


def format_sca_report_status(row):
    mode = text_or_empty(row.get("sca_output_mode", ""))
    if not mode:
        return "not_available"
    state = text_or_empty(row.get("sca_candidate_state", "")) or "not_available"
    tier = text_or_empty(row.get("sca_confidence_tier", "")) or "not_available"
    layer_class = text_or_empty(row.get("sca_report_layer_class", "")) or "not_available"
    text_status = text_or_empty(row.get("report_text_status", "")) or "not_available"
    return f"mode={mode}; state={state}; tier={tier}; report_layer={layer_class}; status={text_status}"


def summarize_branch_b_events(events_df):
    ranked_df = normalize_branch_b_report_events(events_df)
    if ranked_df.empty:
        return pd.DataFrame(columns=["sample_id"]), pd.DataFrame(columns=["sample_id", "branch_b_top_event"])
    if "cnvseq_report_tier" not in ranked_df.columns:
        ranked_df["cnvseq_report_tier"] = ""
    if "cnvseq_reportable" not in ranked_df.columns:
        ranked_df["cnvseq_reportable"] = 0
    ranked_df["cnvseq_reportable"] = pd.to_numeric(ranked_df["cnvseq_reportable"], errors="coerce").fillna(0).astype(int)
    ranked_df["cnvseq_tier_rank"] = ranked_df["cnvseq_report_tier"].map(cnvseq_tier_rank)
    kept_mask = ranked_df["keep_event"].astype(int).eq(1)
    chrom_series = ranked_df["chrom"].fillna("").astype(str) if "chrom" in ranked_df.columns else pd.Series("", index=ranked_df.index)
    ranked_df["_sex_review_top_candidate"] = ranked_df.apply(is_suppressed_sex_review_event, axis=1)
    ranked_df["_sample_has_nonsex_kept_event"] = (
        (kept_mask & ~chrom_series.isin(SEX_CHROMS))
        .groupby(ranked_df["sample_id"], dropna=False)
        .transform("max")
        .astype(bool)
    )
    ranked_df["branch_b_top_display_suppressed"] = (
        ranked_df["_sex_review_top_candidate"] & ranked_df["_sample_has_nonsex_kept_event"]
    )

    sample_df = (
        ranked_df.groupby("sample_id", dropna=False)
        .agg(
            branch_b_total_events=("event_id", "size"),
            branch_b_kept_events=("keep_event", "sum"),
            branch_b_pass_events=("artifact_status", lambda values: int((values == "pass").sum())),
            branch_b_review_events=("artifact_status", lambda values: int((values == "review").sum())),
            branch_b_reportable_events=("cnvseq_reportable", lambda values: int(values.loc[kept_mask.reindex(values.index, fill_value=False)].sum())),
            branch_b_review_tier_events=(
                "cnvseq_report_tier",
                lambda values: int(
                    (
                        values.loc[kept_mask.reindex(values.index, fill_value=False)]
                        .fillna("")
                        .astype(str)
                        .eq("review_1_2mb")
                    ).sum()
                ),
            ),
            branch_b_subreportable_events=(
                "cnvseq_report_tier",
                lambda values: int(
                    (
                        values.loc[kept_mask.reindex(values.index, fill_value=False)]
                        .fillna("")
                        .astype(str)
                        .eq("subreportable_lt1mb")
                    ).sum()
                ),
            ),
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
            ["sample_id", "cnvseq_tier_rank", "priority_score", "n_bins", "end"],
            ascending=[True, True, False, False, False],
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
                "cnvseq_report_tier",
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
                *[
                    column
                    for column in [
                        "cytoband",
                        "genes",
                        "gene_number",
                        "gene_location",
                        "omim_genes",
                        "omim_phenotypes",
                        "hpo_terms",
                        "region_context",
                        "annotation_status",
                    ]
                    if column in display_events_df.columns
                ],
            ]
        ]
        .copy()
    )
    top_branch_b["branch_b_top_event"] = top_branch_b.apply(format_branch_b_top_event, axis=1)
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


def ensure_branch_b_v2_sample_universe(sample_df, branch_b_v2_df):
    """Keep samples with zero final report events visible in the sample report."""
    if branch_b_v2_df.empty or "sample_id" not in branch_b_v2_df.columns:
        return sample_df

    sample_ids = branch_b_v2_df[["sample_id"]].copy()
    if not sample_df.empty and "sample_id" in sample_df.columns:
        sample_ids = pd.concat([sample_ids, sample_df[["sample_id"]]], ignore_index=True)
    sample_ids["sample_id"] = sample_ids["sample_id"].fillna("").astype(str)
    sample_ids = sample_ids.loc[sample_ids["sample_id"].ne("")].drop_duplicates().reset_index(drop=True)
    if sample_ids.empty:
        return sample_df

    if sample_df.empty or "sample_id" not in sample_df.columns:
        expanded = sample_ids.copy()
    else:
        expanded = sample_ids.merge(sample_df, on="sample_id", how="left")

    for column, default in BRANCH_B_SAMPLE_DEFAULTS.items():
        if column not in expanded.columns:
            expanded[column] = default
        expanded[column] = expanded[column].fillna(default)
    return expanded


def format_technical_conclusion(row):
    suppressed_count = int(row.get("branch_b_suppressed_sex_review_events", 0) or 0)
    top_event = text_or_empty(row.get("branch_b_top_event", "")) or ("suppressed_sex_chromosome_review" if suppressed_count > 0 else "none")
    top_fraction = text_or_empty(row.get("branch_b_top_fraction", "")) or "not_estimated"
    top_flags = text_or_empty(row.get("branch_b_top_flags", "")) or ("sex_chromosome_event" if suppressed_count > 0 else "none")
    top_downgrade = text_or_empty(row.get("branch_b_top_downgrade_reason", "")) or ("sex_chromosome_review_suppressed" if suppressed_count > 0 else "none")
    a_branch_count = int(row.get("a_branch_event_count", 0) or 0)
    a_branch_top_event = text_or_empty(row.get("a_branch_top_event", "")) or "none"
    a_branch_review_count = int(row.get("a_branch_review_candidate_count", 0) or 0)
    a_branch_strong_count = int(row.get("a_branch_strong_signal_count", 0) or 0)
    a_branch_shortlist = text_or_empty(row.get("a_branch_review_shortlist", "")) or "none"
    parts = [
        f"Legacy/current-code Branch B kept {int(row['branch_b_kept_events'])} events",
        f"legacy top event: {top_event}",
        f"fraction={top_fraction}",
        f"flags={top_flags}",
        f"downgrade={top_downgrade}",
        f"sex_review_suppressed={suppressed_count}",
        f"A-branch_sensitive_candidates={a_branch_count}",
        f"A-branch_review_shortlist_top{A_BRANCH_REVIEW_SHORTLIST_SIZE}={a_branch_review_count}",
        f"A-branch_strong_signals_z{A_BRANCH_STRONG_SIGNAL_Z:g}={a_branch_strong_count}",
        f"A-branch_top_signal={a_branch_top_event}",
        f"A-branch_review_shortlist={a_branch_shortlist}",
    ]
    sample_qc_status = text_or_empty(row.get("sample_report_qc_status", ""))
    if sample_qc_status:
        parts.extend(
            [
                f"sample_report_qc={sample_qc_status}",
                f"sample_report_qc_reasons={text_or_empty(row.get('sample_report_qc_reasons', 'PASS')) or 'PASS'}",
                f"sample_report_qc_recommendation={text_or_empty(row.get('rebuild_or_resample_recommendation', 'not_recorded')) or 'not_recorded'}",
            ]
        )
    evidence_count = row.get("branch_b_evidence_candidate_count", "")
    if not pd.isna(evidence_count) and evidence_count != "":
        parts.extend(
            [
                f"P3_branch_b_evidence_candidates={int(evidence_count)}",
                f"P3_review_required={int(row.get('branch_b_evidence_review_required_count', 0) or 0)}",
                f"P3_background={text_or_empty(row.get('branch_b_evidence_background_status', 'not_recorded')) or 'not_recorded'}",
                f"P3_report_impact={text_or_empty(row.get('branch_b_evidence_final_report_impact', 'not_recorded')) or 'not_recorded'}",
            ]
        )
    v2_status = row.get("branch_b_v2_burden_status", "")
    if not pd.isna(v2_status) and v2_status != "":
        parts.extend(
            [
                f"Bv2_burden_status={text_or_empty(v2_status)}",
                f"Bv2_background_unknown_review={int(row.get('branch_b_v2_background_unknown_review_count', 0) or 0)}",
                f"Bv2_branch_s_review={int(row.get('branch_b_v2_branch_s_review_count', 0) or 0)}",
                f"Bv2_technical_risk_review={int(row.get('branch_b_v2_technical_risk_review_count', 0) or 0)}",
                f"Bv2_report_candidate={int(row.get('branch_b_v2_report_candidate_count', 0) or 0)}",
                f"Bv2_report_event={int(row.get('branch_b_v2_report_event_count', 0) or 0)}",
                f"Bv2_report_strong={int(row.get('branch_b_v2_report_strong_event_count', 0) or 0)}",
                f"Bv2_report_weak={int(row.get('branch_b_v2_report_weak_event_count', 0) or 0)}",
                f"Bv2_internal_review_event={int(row.get('branch_b_v2_internal_review_event_count', 0) or 0)}",
                f"Bv2_filtered_event={int(row.get('branch_b_v2_filtered_event_count', 0) or 0)}",
                f"Bv2_branch_s_event={int(row.get('branch_b_v2_branch_s_event_count', 0) or 0)}",
                f"Bv2_sample_report_burden_flag={int(row.get('branch_b_v2_sample_report_burden_flag', 0) or 0)}",
                f"Bv2_final_impact={text_or_empty(row.get('branch_b_v2_final_impact', 'development_review_only')) or 'development_review_only'}",
                f"Bv2_legacy_fields_used={bool(row.get('branch_b_v2_legacy_fields_used', False))}",
            ]
        )
    sca_mode = text_or_empty(row.get("sca_output_mode", ""))
    if sca_mode:
        parts.extend(
            [
                f"SCA_mode={sca_mode}",
                f"SCA_state={text_or_empty(row.get('sca_candidate_state', 'not_available')) or 'not_available'}",
                f"SCA_tier={text_or_empty(row.get('sca_confidence_tier', 'not_available')) or 'not_available'}",
                f"SCA_report_layer={text_or_empty(row.get('sca_report_layer_class', 'not_available')) or 'not_available'}",
                f"SCA_status={text_or_empty(row.get('report_text_status', 'not_available')) or 'not_available'}",
            ]
        )
    parts.append(f"QC={text_or_empty(row.get('qc_status', 'NA')) or 'NA'}")
    parts.append(f"sex={text_or_empty(row.get('sex_call', 'NA')) or 'NA'}")
    return "; ".join(parts)


def format_biological_candidate_conclusion(row):
    v2_status = text_or_empty(row.get("branch_b_v2_burden_status", ""))
    if v2_status and v2_status != "not_configured":
        report_count = int(row.get("branch_b_v2_report_event_count", 0) or 0)
        internal_count = int(row.get("branch_b_v2_internal_review_event_count", 0) or 0)
        filtered_count = int(row.get("branch_b_v2_filtered_event_count", 0) or 0)
        branch_s_count = int(row.get("branch_b_v2_branch_s_event_count", 0) or 0)
        impact = text_or_empty(row.get("branch_b_v2_final_impact", "development_review_only")) or "development_review_only"
        if report_count > 0:
            return (
                "Branch B V2 report-layer: "
                f"report_events={report_count}; "
                f"internal_review_events={internal_count}; "
                f"filtered_events={filtered_count}; "
                f"branch_s_events={branch_s_count}; "
                f"impact={impact}"
            )
        if internal_count > 0 or filtered_count > 0 or branch_s_count > 0:
            return (
                "No Branch B V2 report-layer autosomal CNV event; "
                f"internal_review_events={internal_count}; "
                f"filtered_events={filtered_count}; "
                f"branch_s_events={branch_s_count}; "
                f"impact={impact}"
            )

    branch_b_top_event = text_or_empty(row.get("branch_b_top_event", ""))
    if branch_b_top_event:
        return branch_b_top_event
    a_branch_top_event = text_or_empty(row.get("a_branch_top_event", ""))
    if a_branch_top_event:
        if int(row.get("a_branch_strong_signal_count", 0) or 0) > 0:
            return f"A-branch strong sensitive signal only: {a_branch_top_event}; requires Branch B review"
        return f"A-branch sensitive signal only: {a_branch_top_event}; requires Branch B review"
    if int(row.get("branch_b_suppressed_sex_review_events", 0) or 0) > 0:
        return "Branch B top event suppressed (sex-chromosome review only)"
    return "No A-branch event"


def main():
    args = parse_args()
    logger = setup_logger("cnv_report", args.log or None)
    events_df = load_events(args.event_tsv)
    events_df = apply_event_annotations(events_df, args.event_annotation_tsv)
    gender_df = load_one_row_tables(args.gender_tsv)
    qc_df = load_one_row_tables(args.qc_tsv)
    a_branch_df = load_a_branch(args.a_branch_bed)
    branch_a_validation_summaries = load_branch_a_validation_summaries(args.branch_a_validation_summary)
    branch_b_evidence_df = load_branch_b_evidence_summaries(args.branch_b_evidence_summary)
    branch_s_df = load_branch_s_summaries(args.branch_s_summary)
    branch_b_v2_df, branch_b_v2_contract = load_branch_b_v2_burden_summaries(
        sample_summary_path=args.branch_b_v2_sample_summary,
        benchmark_summary_path=args.branch_b_v2_benchmark_summary,
        report_reference_id=args.reference_id,
    )
    sample_report_qc_df = load_sample_report_qc(args.sample_report_qc)
    sample_report_qc_contract = summarize_sample_report_qc(sample_report_qc_df)
    event_annotation_contract = summarize_event_annotation_status(args.event_annotation_tsv)
    evaluation_summary = read_optional_json(args.evaluation_summary)
    ml_summary = read_optional_json(args.ml_summary)
    benchmark_summary = read_optional_json(args.benchmark_summary)
    truth_validation_summary = read_optional_json(args.truth_validation_summary)
    plot_lookup = build_plot_lookup(args.plot_svg)
    copy_number_plot_lookup = build_plot_lookup(args.copy_number_plot_svg)
    fraction_benchmark = benchmark_summary.get("fraction_estimation", {}) if isinstance(benchmark_summary, dict) else {}
    low_fraction_detection = benchmark_summary.get("low_fraction_detection", []) if isinstance(benchmark_summary, dict) else []

    if events_df.empty and branch_b_v2_df.empty and sample_report_qc_df.empty:
        empty = pd.DataFrame(columns=["sample_id"])
        ensure_parent(args.output_tsv)
        empty.to_csv(args.output_tsv, sep="\t", index=False)
        ensure_parent(args.output_json).write_text(json.dumps({"status": "empty"}, indent=2), encoding="utf-8")
        ensure_parent(args.output_md).write_text("# CNV Report\n\nNo events available.\n", encoding="utf-8")
        ensure_parent(args.output_html).write_text("<html><body><h1>CNV Report</h1><p>No events available.</p></body></html>", encoding="utf-8")
        logger.info("report skipped: no events")
        return

    sample_df, top_branch_b = summarize_branch_b_events(events_df)
    sample_df = ensure_branch_b_v2_sample_universe(sample_df, branch_b_v2_df)
    sample_df = ensure_branch_b_v2_sample_universe(sample_df, sample_report_qc_df)

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
    if not branch_b_evidence_df.empty:
        sample_df = sample_df.merge(branch_b_evidence_df, on="sample_id", how="left")
    if not branch_s_df.empty:
        sample_df = sample_df.merge(branch_s_df, on="sample_id", how="left")
    if not branch_b_v2_df.empty:
        sample_df = sample_df.merge(branch_b_v2_df, on="sample_id", how="left")
    if not sample_report_qc_df.empty:
        sample_df = sample_df.merge(sample_report_qc_df, on="sample_id", how="left")

    sample_df["technical_conclusion"] = sample_df.apply(format_technical_conclusion, axis=1)
    sample_df["biological_candidate_conclusion"] = sample_df.apply(format_biological_candidate_conclusion, axis=1)
    sample_df["branch_b_evidence_status"] = sample_df.apply(format_branch_b_evidence_status, axis=1)
    sample_df["branch_b_v2_burden_display"] = sample_df.apply(format_branch_b_v2_burden_status, axis=1)
    sample_df["sample_report_qc_display"] = sample_df.apply(format_sample_report_qc_status, axis=1)
    sample_df["sca_report_status"] = sample_df.apply(format_sca_report_status, axis=1)
    sample_df["plot_svg"] = sample_df["sample_id"].map(lambda sample_id: plot_lookup.get(str(sample_id), ""))
    sample_df["copy_number_plot_svg"] = sample_df["sample_id"].map(lambda sample_id: copy_number_plot_lookup.get(str(sample_id), ""))
    sample_df = sample_df.drop(columns=[column for column in A_BRANCH_INTERNAL_COLUMNS if column in sample_df.columns])

    ensure_parent(args.output_tsv)
    sample_df.to_csv(args.output_tsv, sep="\t", index=False)
    branch_a_validation_gate = summarize_branch_a_validation_gate(args.reference_id, branch_a_validation_summaries)
    report_contract = {
        "status": "development_only_not_final_release",
        "reference_id": str(args.reference_id or "not_configured"),
        "wisecondorx_predict_command": str(args.wisecondorx_predict_command or "not_configured"),
        "branch_a_validation_summary_count": int(len(branch_a_validation_summaries)),
        "branch_a_no_fn_status": branch_a_validation_gate["status"],
        "same_reference_config_status": branch_a_validation_gate["same_reference_config_status"],
        "branch_b_evidence_summary_count": int(len(branch_b_evidence_df)),
        "branch_s_summary_count": int(len(branch_s_df)),
        **branch_b_v2_contract,
        "branch_b_raw_ledger_used": False,
        "branch_s_raw_evidence_used": False,
        **event_annotation_contract,
        **sample_report_qc_contract,
        "note": (
            "P6 report package carries P3/P5 review summaries and Branch B V2 burden display only; "
            "it does not promote Branch B V2, Branch S, or SCA."
        ),
    }
    payload = {
        "status": "completed",
        "report_contract": report_contract,
        "reference_id": report_contract["reference_id"],
        "wisecondorx_predict_command": report_contract["wisecondorx_predict_command"],
        "branch_a_validation_gate": branch_a_validation_gate,
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
        f"- Report release status: `{report_contract['status']}`",
        f"- Reference ID: `{report_contract['reference_id']}`",
        f"- Branch A no-FN status: `{branch_a_validation_gate['status']}`",
        f"- Branch A no-FN truth: `{branch_a_validation_gate['truth_detected_count']}/{branch_a_validation_gate['truth_event_count']}`",
        f"- Branch A no-FN FN count: `{branch_a_validation_gate['FN_count']}`",
        f"- H6 chr21 status: `{branch_a_validation_gate['H6_chr21_status']}`",
        f"- Branch B evidence summaries: `{report_contract['branch_b_evidence_summary_count']}`",
        f"- Branch S summaries: `{report_contract['branch_s_summary_count']}`",
        f"- Branch B V2 burden status: `{report_contract['branch_b_v2_burden_status']}`",
        f"- Branch B V2 final impact: `{report_contract['branch_b_v2_final_impact']}`",
        f"- Branch B V2 legacy fields used: `{report_contract['branch_b_v2_legacy_fields_used']}`",
        f"- Branch B V2 background-unknown review burden: `{report_contract['branch_b_v2_background_unknown_review_count']}`",
        f"- Branch B V2 Branch S review burden: `{report_contract['branch_b_v2_branch_s_review_count']}`",
        f"- Branch B V2 technical-risk review burden: `{report_contract['branch_b_v2_technical_risk_review_count']}`",
        f"- Branch B V2 report-candidate burden: `{report_contract['branch_b_v2_report_candidate_count']}`",
        f"- Branch B V2 final report events: `{report_contract['branch_b_v2_report_event_count']}`",
        f"- Branch B V2 final report strong events: `{report_contract['branch_b_v2_report_strong_event_count']}`",
        f"- Branch B V2 final report weak events: `{report_contract['branch_b_v2_report_weak_event_count']}`",
        f"- Branch B V2 internal review events: `{report_contract['branch_b_v2_internal_review_event_count']}`",
        f"- Branch B V2 filtered events: `{report_contract['branch_b_v2_filtered_event_count']}`",
        f"- Branch B V2 Branch S events: `{report_contract['branch_b_v2_branch_s_event_count']}`",
        f"- Branch B V2 sample report-burden flags: `{report_contract['branch_b_v2_sample_report_burden_flag_count']}`"
        f" (threshold `{report_contract['branch_b_v2_sample_report_burden_threshold']}` report events/sample; audit-only)",
        "- Branch B V2 limitation: `development_review_only display; not FP-reduction proof; not final promotion`",
        f"- Sample reportability QC status: `{report_contract['sample_report_qc_status']}`",
        f"- Sample reportability QC: pass `{report_contract['pass_reportable_count']}`, review `{report_contract['sample_quality_review_count']}`, no-call `{report_contract['no_call_recommended_count']}`",
        f"- Event annotation status: `{report_contract['event_annotation_status']}`",
        f"- Event annotation rows: `{report_contract['event_annotation_row_count']}`",
        f"- Event annotation backend: `{report_contract['event_annotation_backend']}`",
        f"- Event annotation gene source: `{report_contract['event_annotation_gene_source_status']}`",
        f"- Event annotation OMIM source: `{report_contract['event_annotation_omim_source_status']}`",
        f"- Event annotation HPO source: `{report_contract['event_annotation_hpo_source_status']}`",
        "",
        "## Sample Table",
        "",
        "| Sample | QC | Reportability QC | Sex | Plot | CN Plot | Legacy Branch B Top Event | P3 Evidence | Branch B V2 Burden | SCA Status | Technical Conclusion | Biological Candidate Conclusion |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
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
        plot_link = format_plot_link(row.sample_id, args.output_md, plot_lookup, html=False, label="z plot")
        cn_plot_link = format_plot_link(row.sample_id, args.output_md, copy_number_plot_lookup, html=False, label="CN plot")
        md_lines.append(
            f"| `{row.sample_id}` | `{getattr(row, 'qc_status', 'NA')}` | `{getattr(row, 'sample_report_qc_display', 'not_available')}` | `{getattr(row, 'sex_call', 'NA')}` | "
            f"{plot_link or ''} | {cn_plot_link or ''} | "
            f"{getattr(row, 'branch_b_top_event', '') or 'none'} | "
            f"{getattr(row, 'branch_b_evidence_status', 'not_available')} | "
            f"{getattr(row, 'branch_b_v2_burden_display', 'not_available')} | "
            f"{getattr(row, 'sca_report_status', 'not_available')} | "
            f"{row.technical_conclusion} | {row.biological_candidate_conclusion} |"
        )
    ensure_parent(args.output_md).write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    html_rows = []
    for row in sample_df.itertuples(index=False):
        plot_link = format_plot_link(row.sample_id, args.output_html, plot_lookup, html=True, label="z plot")
        cn_plot_link = format_plot_link(row.sample_id, args.output_html, copy_number_plot_lookup, html=True, label="CN plot")
        html_rows.append(
            "<tr>"
            f"<td>{html_lib.escape(str(row.sample_id))}</td>"
            f"<td>{html_lib.escape(str(getattr(row, 'qc_status', 'NA')))}</td>"
            f"<td>{html_lib.escape(str(getattr(row, 'sample_report_qc_display', 'not_available')))}</td>"
            f"<td>{html_lib.escape(str(getattr(row, 'sex_call', 'NA')))}</td>"
            f"<td>{plot_link}</td>"
            f"<td>{cn_plot_link}</td>"
            f"<td>{html_lib.escape(str(getattr(row, 'branch_b_top_event', '') or 'none'))}</td>"
            f"<td>{html_lib.escape(str(getattr(row, 'branch_b_evidence_status', 'not_available')))}</td>"
            f"<td>{html_lib.escape(str(getattr(row, 'branch_b_v2_burden_display', 'not_available')))}</td>"
            f"<td>{html_lib.escape(str(getattr(row, 'sca_report_status', 'not_available')))}</td>"
            f"<td>{html_lib.escape(str(row.technical_conclusion))}</td>"
            f"<td>{html_lib.escape(str(row.biological_candidate_conclusion))}</td>"
            "</tr>"
        )
    html = (
        "<html><head><meta charset='utf-8'><title>CNV Report</title>"
        "<style>body{font-family:Arial,sans-serif;margin:24px;}table{border-collapse:collapse;width:100%;}"
        "th,td{border:1px solid #ccc;padding:8px;vertical-align:top;}th{background:#f3f3f3;text-align:left;}</style>"
        "</head><body><h1>CNV Report</h1>"
        f"<p>Samples: {sample_df['sample_id'].nunique()} | Branch B kept events: {int(sample_df['branch_b_kept_events'].sum())} | "
        f"Report release status: {html_lib.escape(report_contract['status'])}</p>"
        f"<p>Branch B V2 report-layer events: report={int(report_contract['branch_b_v2_report_event_count'])}; "
        f"report_strong={int(report_contract['branch_b_v2_report_strong_event_count'])}; "
        f"report_weak={int(report_contract['branch_b_v2_report_weak_event_count'])}; "
        f"internal_review={int(report_contract['branch_b_v2_internal_review_event_count'])}; "
        f"filtered={int(report_contract['branch_b_v2_filtered_event_count'])}; "
        f"branch_s={int(report_contract['branch_b_v2_branch_s_event_count'])}; "
        f"sample_report_burden_flags={int(report_contract['branch_b_v2_sample_report_burden_flag_count'])}</p>"
        f"<p>Sample reportability QC: status={html_lib.escape(str(report_contract['sample_report_qc_status']))}; "
        f"pass={int(report_contract['pass_reportable_count'])}; "
        f"review={int(report_contract['sample_quality_review_count'])}; "
        f"no_call={int(report_contract['no_call_recommended_count'])}</p>"
        "<p>Branch B V2 burden display is development_review_only evidence; it is not FP-reduction proof and is not final promotion.</p>"
        "<table><thead><tr><th>Sample</th><th>QC</th><th>Reportability QC</th><th>Sex</th><th>Plot</th><th>CN Plot</th><th>Legacy Branch B Top Event</th>"
        "<th>P3 Evidence</th><th>Branch B V2 Burden</th><th>SCA Status</th><th>Technical Conclusion</th><th>Biological Candidate Conclusion</th></tr></thead><tbody>"
        + "".join(html_rows)
        + "</tbody></table></body></html>"
    )
    ensure_parent(args.output_html).write_text(html, encoding="utf-8")
    logger.info("report completed: samples=%d", sample_df["sample_id"].nunique())


if __name__ == "__main__":
    main()
