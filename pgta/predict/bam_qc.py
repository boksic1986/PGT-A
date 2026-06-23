#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import pandas as pd

from pgta.qc.report import parse_fastp_json, parse_flagstat, parse_idxstats


BAM_QC_COLUMNS = [
    "sample_id",
    "bam_qc_status",
    "bam_qc_reasons",
    "clean_reads",
    "clean_gc_content",
    "duplication_rate",
    "insert_peak",
    "mapping_rate_bam",
    "proper_pair_rate_bam",
    "autosome_mapped_reads",
    "autosome_bias_cv",
    "chrX_mapped_reads",
    "chrY_mapped_reads",
    "library_qc_note",
]


def _safe_float(value, default=math.nan):
    try:
        if value is None:
            return default
        number = float(value)
    except (TypeError, ValueError):
        return default
    return number if math.isfinite(number) else default


def _is_missing(value) -> bool:
    return not math.isfinite(_safe_float(value))


def _chrom_token(chrom: str) -> str:
    text = str(chrom)
    return text[3:] if text.startswith("chr") else text


def _is_autosome(chrom: str) -> bool:
    token = _chrom_token(chrom)
    return token.isdigit() and 1 <= int(token) <= 22


def summarize_idxstats_autosomes(rows) -> dict[str, float]:
    df = pd.DataFrame(rows)
    if df.empty:
        return {
            "autosome_mapped_reads": 0,
            "autosome_bias_cv": math.nan,
            "chrX_mapped_reads": 0,
            "chrY_mapped_reads": 0,
        }
    for column in ["length", "mapped_reads"]:
        if column not in df.columns:
            df[column] = 0
        df[column] = pd.to_numeric(df[column], errors="coerce").fillna(0)
    df["chrom"] = df["chrom"].astype(str)
    auto = df[df["chrom"].map(_is_autosome)].copy()
    if auto.empty:
        auto_depth = pd.Series(dtype=float)
    else:
        auto = auto[auto["length"] > 0].copy()
        auto_depth = auto["mapped_reads"] / auto["length"]
    mean_depth = float(auto_depth.mean()) if len(auto_depth) else math.nan
    bias_cv = float(auto_depth.std(ddof=0) / mean_depth) if mean_depth and math.isfinite(mean_depth) else math.nan
    chr_x = int(df.loc[df["chrom"].isin(["chrX", "X"]), "mapped_reads"].sum())
    chr_y = int(df.loc[df["chrom"].isin(["chrY", "Y"]), "mapped_reads"].sum())
    return {
        "autosome_mapped_reads": int(auto["mapped_reads"].sum()) if not auto.empty else 0,
        "autosome_bias_cv": bias_cv,
        "chrX_mapped_reads": chr_x,
        "chrY_mapped_reads": chr_y,
    }


def classify_predict_bam_qc(metrics: dict, thresholds: dict | None = None) -> dict:
    thresholds = thresholds or {}
    sample_id = str(metrics.get("sample_id", ""))
    reasons: list[str] = []
    fail = False

    clean_reads = _safe_float(metrics.get("clean_reads"))
    clean_gc = _safe_float(metrics.get("clean_gc_content"))
    duplication = _safe_float(metrics.get("duplication_rate"))
    insert_peak = _safe_float(metrics.get("insert_peak"))
    mapping_rate = _safe_float(metrics.get("mapping_rate_bam"))
    proper_pair_rate = _safe_float(metrics.get("proper_pair_rate_bam"))
    autosome_bias_cv = _safe_float(metrics.get("autosome_bias_cv"))

    if _is_missing(mapping_rate):
        fail = True
        reasons.append("FLAGSTAT_MAPPING_RATE_MISSING")
    elif mapping_rate < float(thresholds.get("mapping_rate_fail", 0.85)):
        fail = True
        reasons.append("LOW_MAPPING_RATE_FAIL")
    elif mapping_rate < float(thresholds.get("mapping_rate_review", 0.95)):
        reasons.append("LOW_MAPPING_RATE_REVIEW")

    if _is_missing(proper_pair_rate):
        fail = True
        reasons.append("FLAGSTAT_PROPER_PAIR_RATE_MISSING")
    elif proper_pair_rate < float(thresholds.get("proper_pair_rate_fail", 0.70)):
        fail = True
        reasons.append("LOW_PROPER_PAIR_RATE_FAIL")
    elif proper_pair_rate < float(thresholds.get("proper_pair_rate_review", 0.90)):
        reasons.append("LOW_PROPER_PAIR_RATE_REVIEW")

    if _is_missing(clean_reads) or _is_missing(clean_gc):
        reasons.append("FASTP_METRICS_MISSING")
    else:
        if clean_reads < float(thresholds.get("clean_reads_fail", 500_000)):
            fail = True
            reasons.append("LOW_CLEAN_READS_FAIL")
        elif clean_reads < float(thresholds.get("clean_reads_review", 1_000_000)):
            reasons.append("LOW_CLEAN_READS_REVIEW")
        if clean_gc < float(thresholds.get("gc_fail_low", 0.20)) or clean_gc > float(thresholds.get("gc_fail_high", 0.70)):
            fail = True
            reasons.append("GC_CONTENT_FAIL")
        elif clean_gc < float(thresholds.get("gc_review_low", 0.30)) or clean_gc > float(thresholds.get("gc_review_high", 0.60)):
            reasons.append("GC_CONTENT_REVIEW")

    if not _is_missing(duplication):
        if duplication > float(thresholds.get("duplication_fail", 0.95)):
            fail = True
            reasons.append("DUPLICATION_RATE_FAIL")
        elif duplication > float(thresholds.get("duplication_review", 0.80)):
            reasons.append("DUPLICATION_RATE_REVIEW")

    if not _is_missing(insert_peak) and insert_peak <= 0:
        reasons.append("INSERT_SIZE_UNINFORMATIVE_REVIEW")

    if _is_missing(autosome_bias_cv):
        reasons.append("IDXSTATS_AUTOSOME_BIAS_MISSING")
    elif autosome_bias_cv > float(thresholds.get("autosome_bias_cv_fail", 0.50)):
        fail = True
        reasons.append("AUTOSOME_CHROMOSOME_BIAS_FAIL")
    elif autosome_bias_cv > float(thresholds.get("autosome_bias_cv_review", 0.25)):
        reasons.append("AUTOSOME_CHROMOSOME_BIAS_REVIEW")

    if fail:
        status = "BAM_QC_FAIL"
    elif reasons:
        status = "BAM_QC_REVIEW"
    else:
        status = "BAM_QC_PASS"

    row = {column: metrics.get(column, "") for column in BAM_QC_COLUMNS}
    row.update(
        {
            "sample_id": sample_id,
            "bam_qc_status": status,
            "bam_qc_reasons": ";".join(reasons) if reasons else "PASS",
            "library_qc_note": (
                "library_qc_blocks_formal_report"
                if status == "BAM_QC_FAIL"
                else ("library_qc_review_recommended" if status == "BAM_QC_REVIEW" else "library_qc_pass")
            ),
        }
    )
    return row


def build_predict_bam_qc_row(
    *,
    sample_id: str,
    flagstat_path: str,
    idxstats_path: str,
    fastp_json_path: str = "",
    thresholds: dict | None = None,
) -> dict:
    flagstat = parse_flagstat(flagstat_path)
    flagstat["sample_id"] = sample_id
    total_reads = _safe_float(flagstat.get("total_reads_bam"))
    mapped_reads = _safe_float(flagstat.get("mapped_reads_bam"))
    paired_reads = _safe_float(flagstat.get("paired_reads_bam"))
    proper_pair_reads = _safe_float(flagstat.get("properly_paired_reads_bam"))
    mapping_rate = mapped_reads / total_reads if total_reads > 0 else math.nan
    proper_pair_rate = proper_pair_reads / paired_reads if paired_reads > 0 else math.nan

    idx_df = parse_idxstats(idxstats_path)
    idx_summary = summarize_idxstats_autosomes(idx_df.to_dict(orient="records"))

    metrics = {
        "sample_id": sample_id,
        "mapping_rate_bam": mapping_rate,
        "proper_pair_rate_bam": proper_pair_rate,
        **flagstat,
        **idx_summary,
    }
    if fastp_json_path and Path(fastp_json_path).exists():
        fastp = parse_fastp_json(fastp_json_path)
        fastp["sample_id"] = sample_id
        metrics.update(fastp)
    return classify_predict_bam_qc(metrics, thresholds=thresholds)


def write_single_sample(args) -> None:
    row = build_predict_bam_qc_row(
        sample_id=args.sample_id,
        flagstat_path=args.flagstat,
        idxstats_path=args.idxstats,
        fastp_json_path=args.fastp_json or "",
    )
    output_tsv = Path(args.output_tsv)
    output_json = Path(args.output_json)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([row], columns=BAM_QC_COLUMNS).to_csv(output_tsv, sep="\t", index=False)
    output_json.write_text(json.dumps(row, indent=2, ensure_ascii=False), encoding="utf-8")


def write_summary(args) -> None:
    frames = []
    for path_value in args.input_tsv or []:
        path = Path(path_value)
        if path.exists():
            frames.append(pd.read_csv(path, sep="\t"))
    summary = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=BAM_QC_COLUMNS)
    output_tsv = Path(args.summary_tsv)
    output_json = Path(args.summary_json)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output_tsv, sep="\t", index=False)
    payload = {
        "status": "completed",
        "sample_count": int(summary["sample_id"].nunique()) if "sample_id" in summary.columns else 0,
        "bam_qc_pass_count": int((summary.get("bam_qc_status", pd.Series(dtype=str)) == "BAM_QC_PASS").sum()),
        "bam_qc_review_count": int((summary.get("bam_qc_status", pd.Series(dtype=str)) == "BAM_QC_REVIEW").sum()),
        "bam_qc_fail_count": int((summary.get("bam_qc_status", pd.Series(dtype=str)) == "BAM_QC_FAIL").sum()),
        "samples": summary.to_dict(orient="records"),
    }
    output_json.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def parse_args():
    parser = argparse.ArgumentParser(description="Predict-side BAM/library QC.")
    parser.add_argument("--sample-id", default="")
    parser.add_argument("--flagstat", default="")
    parser.add_argument("--idxstats", default="")
    parser.add_argument("--fastp-json", default="")
    parser.add_argument("--output-tsv", default="")
    parser.add_argument("--output-json", default="")
    parser.add_argument("--input-tsv", action="append", default=[])
    parser.add_argument("--summary-tsv", default="")
    parser.add_argument("--summary-json", default="")
    parser.add_argument("--log", default="")
    return parser.parse_args()


def main():
    args = parse_args()
    if args.summary_tsv:
        write_summary(args)
    else:
        missing = [
            name
            for name, value in [("sample-id", args.sample_id), ("flagstat", args.flagstat), ("idxstats", args.idxstats), ("output-tsv", args.output_tsv), ("output-json", args.output_json)]
            if not value
        ]
        if missing:
            raise SystemExit(f"missing required arguments for per-sample BAM QC: {','.join(missing)}")
        write_single_sample(args)


if __name__ == "__main__":
    main()
