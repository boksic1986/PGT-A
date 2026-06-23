#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import pandas as pd


REPORT_QC_COLUMNS = [
    "sample_id",
    "sample_report_qc_status",
    "sample_report_qc_reasons",
    "rebuild_or_resample_recommendation",
    "library_qc_status",
    "global_wave_context",
    "multi_chromosome_cnv_context",
    "possible_contamination_or_mixture_context",
    "cnv_qc_status",
    "mad_log1p",
    "cohort_mad_log1p_median",
    "cohort_mad_log1p_robust_z",
    "report_event_count",
    "internal_review_event_count",
    "filtered_event_count",
    "branch_s_event_count",
    "multi_chromosome_report_event_count",
    "autosomal_cn_outside_1p7_2p3_fraction",
]


def _safe_float(value, default=math.nan):
    try:
        if value is None:
            return default
        number = float(value)
    except (TypeError, ValueError):
        return default
    return number if math.isfinite(number) else default


def _safe_int(value, default=0):
    number = _safe_float(value)
    return int(number) if math.isfinite(number) else default


def _is_truthy_text(value) -> bool:
    return str(value or "").strip().lower() not in {"", "na", "nan", "none", "not_available"}


def classify_sample_reportability(row: dict) -> dict:
    reasons: list[str] = []
    no_call = False

    bam_status = str(row.get("bam_qc_status") or row.get("library_qc_status") or "BAM_QC_NOT_AVAILABLE")
    cnv_status = str(row.get("cnv_qc_status") or row.get("status") or "UNKNOWN")
    cohort_size = _safe_int(row.get("cohort_size"), 0)
    mad_z = _safe_float(row.get("cohort_mad_log1p_robust_z"))
    cn_outside = _safe_float(row.get("autosomal_cn_outside_1p7_2p3_fraction"), 0.0)
    report_count = _safe_int(row.get("report_event_count", row.get("v2_report_event_count")), 0)
    internal_count = _safe_int(row.get("internal_review_event_count", row.get("v2_internal_review_event_count")), 0)
    filtered_count = _safe_int(row.get("filtered_event_count", row.get("v2_filtered_event_count")), 0)
    branch_s_count = _safe_int(row.get("branch_s_event_count", row.get("v2_branch_s_event_count")), 0)
    multi_chrom_count = _safe_int(row.get("multi_chromosome_report_event_count"), 0)

    if bam_status == "BAM_QC_FAIL":
        no_call = True
        reasons.append("BAM_QC_FAIL")
    elif bam_status == "BAM_QC_REVIEW":
        reasons.append("BAM_QC_REVIEW")
    elif bam_status == "BAM_QC_NOT_AVAILABLE":
        reasons.append("BAM_QC_NOT_AVAILABLE_REVIEW")

    if cnv_status.upper() == "FAIL":
        no_call = True
        reasons.append("CNV_INPUT_QC_FAIL")
    elif cnv_status.upper() not in {"PASS", "UNKNOWN", "NA", "NAN", ""}:
        reasons.append(f"CNV_INPUT_QC_{cnv_status.upper()}")

    if math.isfinite(mad_z) and mad_z >= 5.0:
        reasons.append("HIGH_MAD_LOG1P_RELATIVE_OUTLIER")
    if cn_outside >= 0.20:
        reasons.append("AUTOSOMAL_CN_OUTSIDE_1P7_2P3_REVIEW")
    if report_count >= 10 or multi_chrom_count >= 5:
        reasons.append("HIGH_MULTI_CHROMOSOME_REPORT_BURDEN_REVIEW")
    elif report_count >= 5 and multi_chrom_count >= 3:
        reasons.append("MULTI_CHROMOSOME_CNV_REVIEW")

    severe_global = (
        cn_outside >= 0.45
        and report_count >= 20
        and multi_chrom_count >= 10
        and bam_status != "BAM_QC_PASS"
    )
    if severe_global and cohort_size >= 10:
        no_call = True
        reasons.append("GLOBAL_UNINTERPRETABLE_SIGNAL_NO_CALL")

    if no_call:
        status = "NO_CALL_RECOMMENDED"
    elif reasons:
        status = "SAMPLE_QUALITY_REVIEW"
    else:
        status = "PASS_REPORTABLE"

    if status == "NO_CALL_RECOMMENDED":
        recommendation = "do_not_release_formal_cnv_report_rebuild_or_resample"
    elif status == "SAMPLE_QUALITY_REVIEW":
        recommendation = "review_library_or_resequence_before_release"
    else:
        recommendation = "reportable_no_sample_level_qc_action"

    global_wave_context = "GLOBAL_WAVE_REVIEW" if any(
        reason in reasons for reason in ["HIGH_MAD_LOG1P_RELATIVE_OUTLIER", "AUTOSOMAL_CN_OUTSIDE_1P7_2P3_REVIEW"]
    ) else "global_wave_not_flagged"
    multi_chrom_context = (
        "MULTI_CHROMOSOME_CNV_REVIEW"
        if any(reason in reasons for reason in ["HIGH_MULTI_CHROMOSOME_REPORT_BURDEN_REVIEW", "MULTI_CHROMOSOME_CNV_REVIEW"])
        else "multi_chromosome_burden_not_flagged"
    )
    mixture_context = (
        "possible_contamination_or_mixture_review"
        if global_wave_context == "GLOBAL_WAVE_REVIEW" and multi_chrom_context == "MULTI_CHROMOSOME_CNV_REVIEW"
        else "no_mixture_evidence"
    )
    if _is_truthy_text(row.get("chrY_cn_artifact_context")):
        mixture_context = f"{mixture_context}; chrY_context={row.get('chrY_cn_artifact_context')}"

    out = {column: row.get(column, "") for column in REPORT_QC_COLUMNS}
    out.update(
        {
            "sample_id": str(row.get("sample_id", "")),
            "sample_report_qc_status": status,
            "sample_report_qc_reasons": ";".join(reasons) if reasons else "PASS",
            "rebuild_or_resample_recommendation": recommendation,
            "library_qc_status": bam_status,
            "global_wave_context": global_wave_context,
            "multi_chromosome_cnv_context": multi_chrom_context,
            "possible_contamination_or_mixture_context": mixture_context,
            "cnv_qc_status": cnv_status,
            "report_event_count": report_count,
            "internal_review_event_count": internal_count,
            "filtered_event_count": filtered_count,
            "branch_s_event_count": branch_s_count,
            "multi_chromosome_report_event_count": multi_chrom_count,
            "autosomal_cn_outside_1p7_2p3_fraction": cn_outside,
        }
    )
    return out


def _robust_z(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    median = numeric.median()
    mad = (numeric - median).abs().median()
    if not math.isfinite(float(mad or math.nan)) or mad == 0:
        return pd.Series([0.0] * len(numeric), index=numeric.index)
    return (numeric - median) / (1.4826 * mad)


def _load_tsv(path_value: str, columns: list[str] | None = None) -> pd.DataFrame:
    path = Path(path_value)
    if not path.exists():
        return pd.DataFrame(columns=columns or [])
    df = pd.read_csv(path, sep="\t")
    return df


def _load_many_tsv(paths: list[str]) -> pd.DataFrame:
    frames = []
    for path_value in paths or []:
        path = Path(path_value)
        if path.exists():
            df = pd.read_csv(path, sep="\t")
            if not df.empty:
                if "sample_id" not in df.columns:
                    df["sample_id"] = path.name.split(".")[0]
                frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _summarize_cn_bins(paths: list[str]) -> pd.DataFrame:
    frames = []
    for path_value in paths or []:
        path = Path(path_value)
        if not path.exists():
            continue
        df = pd.read_csv(path, sep="\t")
        if df.empty or "copy_number" not in df.columns:
            continue
        sample_id = df["sample_id"].iloc[0] if "sample_id" in df.columns and not df["sample_id"].isna().all() else path.name.split(".")[0]
        chrom = df.get("chrom", pd.Series(dtype=str)).astype(str)
        auto = df[chrom.str.match(r"^(chr)?([1-9]|1[0-9]|2[0-2])$")].copy()
        if auto.empty:
            outside_fraction = 0.0
        else:
            cn = pd.to_numeric(auto["copy_number"], errors="coerce")
            valid = cn.dropna()
            outside_fraction = float(((valid < 1.7) | (valid > 2.3)).mean()) if len(valid) else 0.0
        frames.append(
            {
                "sample_id": str(sample_id),
                "autosomal_cn_outside_1p7_2p3_fraction": outside_fraction,
            }
        )
    return pd.DataFrame(frames)


def _summarize_report_events(report_events: pd.DataFrame) -> pd.DataFrame:
    if report_events.empty or "sample_id" not in report_events.columns:
        return pd.DataFrame(columns=["sample_id", "report_event_count", "multi_chromosome_report_event_count"])
    df = report_events.copy()
    chrom_col = "chrom" if "chrom" in df.columns else None
    grouped = df.groupby("sample_id")
    rows = []
    for sample_id, group in grouped:
        chrom_count = int(group[chrom_col].astype(str).nunique()) if chrom_col else 0
        rows.append(
            {
                "sample_id": str(sample_id),
                "report_event_count": int(len(group)),
                "multi_chromosome_report_event_count": chrom_count,
            }
        )
    return pd.DataFrame(rows)


def build_sample_report_qc_table(
    *,
    bam_qc: pd.DataFrame,
    cnv_qc: pd.DataFrame,
    v2_sample_summary: pd.DataFrame | None = None,
    report_events: pd.DataFrame | None = None,
    cn_bins_summary: pd.DataFrame | None = None,
) -> pd.DataFrame:
    sample_ids = []
    for df in [bam_qc, cnv_qc, v2_sample_summary, report_events, cn_bins_summary]:
        if df is not None and not df.empty and "sample_id" in df.columns:
            sample_ids.extend([str(value) for value in df["sample_id"].dropna().unique()])
    base = pd.DataFrame({"sample_id": sorted(set(sample_ids))})
    if base.empty:
        return pd.DataFrame(columns=REPORT_QC_COLUMNS)
    merged = base
    if not bam_qc.empty:
        merged = merged.merge(bam_qc, on="sample_id", how="left")
    if not cnv_qc.empty:
        keep = [column for column in ["sample_id", "status", "mad_log1p", "nonzero_fraction", "total_counts"] if column in cnv_qc.columns]
        merged = merged.merge(cnv_qc[keep].rename(columns={"status": "cnv_qc_status"}), on="sample_id", how="left")
    if v2_sample_summary is not None and not v2_sample_summary.empty:
        rename_map = {
            "v2_report_event_count": "report_event_count",
            "v2_internal_review_event_count": "internal_review_event_count",
            "v2_filtered_event_count": "filtered_event_count",
            "v2_branch_s_event_count": "branch_s_event_count",
        }
        keep = [column for column in ["sample_id", *rename_map.keys(), "report_event_count", "internal_review_event_count", "filtered_event_count", "branch_s_event_count"] if column in v2_sample_summary.columns]
        v2 = v2_sample_summary[keep].rename(columns=rename_map)
        merged = merged.merge(v2, on="sample_id", how="left")
    if report_events is not None and not report_events.empty:
        merged = merged.merge(_summarize_report_events(report_events), on="sample_id", how="left", suffixes=("", "_events"))
        for column in ["report_event_count", "multi_chromosome_report_event_count"]:
            events_column = f"{column}_events"
            if events_column in merged.columns:
                merged[column] = merged[column].fillna(merged[events_column]) if column in merged.columns else merged[events_column]
                merged = merged.drop(columns=[events_column])
    if cn_bins_summary is not None and not cn_bins_summary.empty:
        merged = merged.merge(cn_bins_summary, on="sample_id", how="left")

    if "mad_log1p" in merged.columns:
        merged["cohort_mad_log1p_median"] = pd.to_numeric(merged["mad_log1p"], errors="coerce").median()
        merged["cohort_mad_log1p_robust_z"] = _robust_z(merged["mad_log1p"]).fillna(0.0)
    else:
        merged["cohort_mad_log1p_median"] = math.nan
        merged["cohort_mad_log1p_robust_z"] = 0.0
    merged["cohort_size"] = int(len(merged))

    rows = []
    for row in merged.to_dict(orient="records"):
        rows.append(classify_sample_reportability(row))
    return pd.DataFrame(rows, columns=REPORT_QC_COLUMNS)


def parse_args():
    parser = argparse.ArgumentParser(description="Post-predict sample reportability QC.")
    parser.add_argument("--bam-qc-summary", required=True)
    parser.add_argument("--cnv-qc-tsv", action="append", default=[])
    parser.add_argument("--branch-b-v2-sample-summary", default="")
    parser.add_argument("--report-events", default="")
    parser.add_argument("--copy-number-bin-tsv", action="append", default=[])
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--output-json", required=True)
    parser.add_argument("--log", default="")
    return parser.parse_args()


def main():
    args = parse_args()
    bam_qc = _load_tsv(args.bam_qc_summary)
    cnv_qc = _load_many_tsv(args.cnv_qc_tsv)
    v2_sample_summary = _load_tsv(args.branch_b_v2_sample_summary) if args.branch_b_v2_sample_summary else pd.DataFrame()
    report_events = _load_tsv(args.report_events) if args.report_events else pd.DataFrame()
    cn_bins_summary = _summarize_cn_bins(args.copy_number_bin_tsv)
    table = build_sample_report_qc_table(
        bam_qc=bam_qc,
        cnv_qc=cnv_qc,
        v2_sample_summary=v2_sample_summary,
        report_events=report_events,
        cn_bins_summary=cn_bins_summary,
    )
    output_tsv = Path(args.output_tsv)
    output_json = Path(args.output_json)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_tsv, sep="\t", index=False)
    payload = {
        "status": "completed",
        "sample_count": int(table["sample_id"].nunique()) if "sample_id" in table.columns else 0,
        "pass_reportable_count": int((table.get("sample_report_qc_status", pd.Series(dtype=str)) == "PASS_REPORTABLE").sum()),
        "sample_quality_review_count": int((table.get("sample_report_qc_status", pd.Series(dtype=str)) == "SAMPLE_QUALITY_REVIEW").sum()),
        "no_call_recommended_count": int((table.get("sample_report_qc_status", pd.Series(dtype=str)) == "NO_CALL_RECOMMENDED").sum()),
        "samples": table.to_dict(orient="records"),
    }
    output_json.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


if __name__ == "__main__":
    main()
