#!/biosoftware/miniconda/envs/snakemake_env/bin/python
from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd

from pgta.core.logging import setup_logger
from pgta.predict.branch_b.common import read_table, write_json, write_table


LEDGER_COLUMNS = [
    "sample_id",
    "candidate_id",
    "reference_id",
    "negative_bank_version",
    "chrom",
    "start",
    "end",
    "state",
    "a_zscore",
    "a_abs_zscore",
    "a_ratio",
    "a_support_level",
    "branch_b_event_id",
    "branch_b_report_class",
    "branch_b_artifact_status",
    "branch_b_keep_event",
    "final_disposition",
    "raw_amplitude",
    "corrected_amplitude",
    "attenuation_ratio",
    "same_direction_fraction",
    "flank_contrast",
    "hard_region_fraction",
    "sample_noise_mad",
    "refmap_status",
    "calibration_null_status",
    "evidence_missing_reason",
]


def parse_args():
    parser = argparse.ArgumentParser(description="Build Branch B V2 shadow evidence rows for every Branch A candidate.")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--a-candidates", required=True)
    parser.add_argument("--branch-b-events", default="")
    parser.add_argument("--input-bins", default="")
    parser.add_argument("--output-ledger", required=True)
    parser.add_argument("--output-summary", required=True)
    parser.add_argument("--reference-id", default="UNKNOWN_REFERENCE")
    parser.add_argument("--negative-bank-version", default="UNKNOWN_NEGATIVE_BANK")
    parser.add_argument("--log", default="")
    return parser.parse_args()


def _empty_table(columns=None):
    return pd.DataFrame(columns=list(columns or []))


def read_optional_table(path_value, columns=None):
    if not path_value:
        return _empty_table(columns)
    path = Path(path_value)
    if not path.exists() or path.stat().st_size == 0:
        return _empty_table(columns)
    return read_table(path, empty_ok=True)


def safe_float(value, default=np.nan):
    try:
        result = float(value)
    except (TypeError, ValueError):
        return float(default)
    return result if np.isfinite(result) else float(default)


def state_direction(state):
    text = str(state).strip().lower()
    if text == "gain":
        return 1.0
    if text == "loss":
        return -1.0
    return 0.0


def finite_median_mad(values):
    vector = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(dtype=np.float64)
    vector = vector[np.isfinite(vector)]
    if vector.size == 0:
        return np.nan
    median = float(np.median(vector))
    return float(np.median(np.abs(vector - median)))


def first_numeric(row, columns, default=np.nan):
    for column in columns:
        if column not in row.index:
            continue
        value = safe_float(row[column], default=np.nan)
        if np.isfinite(value):
            return value
    return float(default)


def first_text(row, columns, default=""):
    for column in columns:
        if column not in row.index:
            continue
        value = str(row[column]).strip()
        if value and value.lower() != "nan":
            return value
    return default


def append_reason(reasons, reason):
    if reason and reason not in reasons:
        reasons.append(reason)


def event_overlap_bins(candidate_row, bins_df):
    if bins_df is None or bins_df.empty:
        return pd.DataFrame()
    required = {"chrom", "start", "end"}
    if not required.issubset(bins_df.columns):
        return pd.DataFrame()
    chrom = str(candidate_row["chrom"])
    start = int(candidate_row["start"])
    end = int(candidate_row["end"])
    return bins_df[
        bins_df["chrom"].astype(str).eq(chrom)
        & (pd.to_numeric(bins_df["start"], errors="coerce") < end)
        & (pd.to_numeric(bins_df["end"], errors="coerce") > start)
    ].copy()


def flank_bins(candidate_row, bins_df):
    if bins_df is None or bins_df.empty or not {"chrom", "start", "end"}.issubset(bins_df.columns):
        return pd.DataFrame()
    chrom = str(candidate_row["chrom"])
    start = int(candidate_row["start"])
    end = int(candidate_row["end"])
    span = max(end - start, 1)
    left = max(start - span, 0)
    right = end + span
    chrom_df = bins_df[bins_df["chrom"].astype(str).eq(chrom)].copy()
    return chrom_df[
        ((pd.to_numeric(chrom_df["start"], errors="coerce") >= left) & (pd.to_numeric(chrom_df["end"], errors="coerce") <= start))
        | ((pd.to_numeric(chrom_df["start"], errors="coerce") >= end) & (pd.to_numeric(chrom_df["end"], errors="coerce") <= right))
    ].copy()


def choose_branch_b_event(candidate_row, branch_b_events):
    if branch_b_events is None or branch_b_events.empty:
        return None
    candidate_id = str(candidate_row["candidate_id"])
    matches = pd.DataFrame()
    if "a_candidate_id" in branch_b_events.columns:
        matches = branch_b_events[branch_b_events["a_candidate_id"].astype(str).eq(candidate_id)].copy()
    if matches.empty and "candidate_id" in branch_b_events.columns:
        matches = branch_b_events[branch_b_events["candidate_id"].astype(str).eq(candidate_id)].copy()
    if matches.empty:
        return None
    score = pd.Series(0.0, index=matches.index)
    if "keep_event" in matches.columns:
        score += pd.to_numeric(matches["keep_event"], errors="coerce").fillna(0.0) * 1000.0
    for column in ["event_corr_adjusted_z", "calibrated_mean_z", "segment_mean_robust_z", "a_abs_zscore"]:
        if column in matches.columns:
            score += pd.to_numeric(matches[column], errors="coerce").abs().fillna(0.0)
    return matches.loc[score.sort_values(ascending=False).index[0]]


def signal_column(frame):
    for column in ["calibrated_z", "robust_z", "signal_for_calling", "normalized_signal"]:
        if column in frame.columns:
            return column
    return ""


def signed_mean(frame, direction, preferred_columns):
    column = next((item for item in preferred_columns if item in frame.columns), "")
    if not column:
        return np.nan
    values = pd.to_numeric(frame[column], errors="coerce").to_numpy(dtype=np.float64)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.nan
    return float(direction * np.mean(values))


def calculate_same_direction_fraction(candidate_row, event_bins):
    if event_bins.empty:
        return np.nan
    direction = state_direction(candidate_row["state"])
    column = signal_column(event_bins)
    if direction == 0.0 or not column:
        return np.nan
    values = pd.to_numeric(event_bins[column], errors="coerce").to_numpy(dtype=np.float64)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.nan
    return float(np.mean((values * direction) > 0.0))


def calculate_flank_contrast(candidate_row, event_bins, flanks):
    if event_bins.empty or flanks.empty:
        return np.nan
    direction = state_direction(candidate_row["state"])
    if direction == 0.0:
        return np.nan
    event_mean = signed_mean(event_bins, direction, ["calibrated_z", "robust_z"])
    flank_mean = signed_mean(flanks, direction, ["calibrated_z", "robust_z"])
    if not np.isfinite(event_mean) or not np.isfinite(flank_mean):
        return np.nan
    return float(event_mean - flank_mean)


def calculate_hard_region_fraction(event_row, event_bins):
    if event_row is not None:
        value = first_numeric(event_row, ["hard_region_fraction", "high_risk_bin_fraction"], default=np.nan)
        if np.isfinite(value):
            return float(np.clip(value, 0.0, 1.0))
    if event_bins.empty:
        return np.nan
    hard_parts = []
    if "mask_label" in event_bins.columns:
        hard_parts.append(event_bins["mask_label"].fillna("").astype(str).isin({"hard", "dynamic"}).to_numpy(dtype=bool))
    if "region_risk_class" in event_bins.columns:
        hard_parts.append(event_bins["region_risk_class"].fillna("").astype(str).eq("high").to_numpy(dtype=bool))
    if not hard_parts:
        return np.nan
    combined = np.logical_or.reduce(hard_parts)
    return float(np.mean(combined))


def classify_refmap_status(event_row, event_bins):
    values = {}
    source_available = False
    for column in ["wisecondorx_ref_bin_count", "ref_size_after_cutoff", "reference_bin_count", "ref_bin_count"]:
        if event_row is not None and column in event_row.index and np.isfinite(safe_float(event_row[column], default=np.nan)):
            values["ref_bin_count"] = safe_float(event_row[column])
            source_available = True
            break
        if not event_bins.empty and column in event_bins.columns:
            series = pd.to_numeric(event_bins[column], errors="coerce")
            if series.notna().any():
                values["ref_bin_count"] = float(series.min())
                source_available = True
                break
    for column in ["low_refbin_fraction", "wisecondorx_low_refbin_component"]:
        if event_row is not None and column in event_row.index and np.isfinite(safe_float(event_row[column], default=np.nan)):
            values[column] = safe_float(event_row[column])
            source_available = True
        elif not event_bins.empty and column in event_bins.columns:
            series = pd.to_numeric(event_bins[column], errors="coerce")
            if series.notna().any():
                values[column] = float(series.max())
                source_available = True
    if event_row is not None and "same_chrom_ref_bin_count" in event_row.index and np.isfinite(safe_float(event_row["same_chrom_ref_bin_count"], default=np.nan)):
        values["same_chrom_ref_bin_count"] = safe_float(event_row["same_chrom_ref_bin_count"])
        source_available = True
    elif not event_bins.empty and "same_chrom_ref_bin_count" in event_bins.columns:
        series = pd.to_numeric(event_bins["same_chrom_ref_bin_count"], errors="coerce")
        if series.notna().any():
            values["same_chrom_ref_bin_count"] = float(series.max())
            source_available = True

    if not source_available:
        return "UNKNOWN"
    if values.get("same_chrom_ref_bin_count", 0.0) > 0.0:
        return "SAME_CHROM_REF_LEAKAGE"
    ref_count = values.get("ref_bin_count", np.nan)
    if np.isfinite(ref_count) and ref_count < 50.0:
        return "LOW_REF_BINS_HARD_RISK"
    if np.isfinite(ref_count) and ref_count < 150.0:
        return "LOW_REF_BINS_REVIEW"
    if max(values.get("low_refbin_fraction", 0.0), values.get("wisecondorx_low_refbin_component", 0.0)) >= 0.50:
        return "LOW_REFBIN_BURDEN"
    return "OK"


def classify_calibration_null_status(event_bins):
    if event_bins.empty or "calibration_null_eligible" not in event_bins.columns:
        return "UNKNOWN"
    values = pd.to_numeric(event_bins["calibration_null_eligible"], errors="coerce")
    if values.notna().sum() == 0:
        return "UNKNOWN"
    eligible_fraction = float(values.fillna(0.0).astype(float).mean())
    if eligible_fraction <= 0.0:
        return "NO_NULL_SUPPORT"
    if eligible_fraction < 0.50:
        return "LIMITED_NULL_SUPPORT"
    return "OK"


def derive_disposition(event_row):
    if event_row is None:
        return "REVIEW_REQUIRED"
    report_class = first_text(event_row, ["report_class"], default="").lower()
    artifact_status = first_text(event_row, ["artifact_status"], default="").lower()
    keep_event = first_numeric(event_row, ["keep_event"], default=np.nan)
    mosaic_status = first_text(event_row, ["mosaic_fraction_status", "mosaic_status"], default="").lower()
    if "mosaic" in mosaic_status and "artifact" not in artifact_status:
        return "MOSAIC_SUSPECT"
    if report_class == "candidate_pass" or artifact_status in {"pass", "confirmed"}:
        return "CONFIRMED"
    if report_class == "candidate_suppressed" or artifact_status == "artifact" or keep_event == 0.0:
        return "LIKELY_ARTIFACT"
    if artifact_status in {"no_call", "nocall"}:
        return "NO_CALL"
    return "REVIEW_REQUIRED"


def build_candidate_evidence_ledger(
    sample_id,
    a_candidates,
    branch_b_events=None,
    bins_df=None,
    reference_id="UNKNOWN_REFERENCE",
    negative_bank_version="UNKNOWN_NEGATIVE_BANK",
):
    if a_candidates is None or a_candidates.empty:
        return pd.DataFrame(columns=LEDGER_COLUMNS)
    required = {"candidate_id", "sample_id", "chrom", "start", "end", "state"}
    missing = sorted(required.difference(a_candidates.columns))
    if missing:
        raise ValueError(f"Branch A candidate table missing columns: {','.join(missing)}")
    branch_b_events = branch_b_events if branch_b_events is not None else pd.DataFrame()
    bins_df = bins_df if bins_df is not None else pd.DataFrame()
    sample_noise_column = signal_column(bins_df) if not bins_df.empty else ""
    sample_noise_mad = finite_median_mad(bins_df[sample_noise_column]) if sample_noise_column else np.nan

    rows = []
    for candidate in a_candidates.sort_values(["chrom", "start", "end", "candidate_id"]).iterrows():
        _, candidate_row = candidate
        event_row = choose_branch_b_event(candidate_row, branch_b_events)
        event_bins = event_overlap_bins(candidate_row, bins_df)
        flanks = flank_bins(candidate_row, bins_df)
        raw_amplitude = (
            first_numeric(event_row, ["raw_amplitude", "segment_mean_robust_z"], default=np.nan)
            if event_row is not None
            else np.nan
        )
        if not np.isfinite(raw_amplitude) and not event_bins.empty:
            raw_amplitude = signed_mean(event_bins, 1.0, ["robust_z", "normalized_signal"])
        corrected_amplitude = (
            first_numeric(event_row, ["corrected_amplitude", "event_corr_adjusted_z", "calibrated_mean_z"], default=np.nan)
            if event_row is not None
            else np.nan
        )
        if not np.isfinite(corrected_amplitude) and not event_bins.empty:
            corrected_amplitude = signed_mean(event_bins, 1.0, ["calibrated_z", "robust_z"])
        attenuation_ratio = (
            abs(corrected_amplitude) / abs(raw_amplitude)
            if np.isfinite(corrected_amplitude) and np.isfinite(raw_amplitude) and abs(raw_amplitude) > 1e-12
            else np.nan
        )
        refmap_status = classify_refmap_status(event_row, event_bins)
        calibration_null_status = classify_calibration_null_status(event_bins)
        reasons = []
        if event_row is None:
            append_reason(reasons, "branch_b_event_missing")
        if refmap_status == "UNKNOWN":
            append_reason(reasons, "refmap_missing")
        if calibration_null_status == "UNKNOWN":
            append_reason(reasons, "calibration_null_missing")
        if event_bins.empty:
            append_reason(reasons, "event_bins_missing")

        payload = {
            "sample_id": str(sample_id or candidate_row["sample_id"]),
            "candidate_id": str(candidate_row["candidate_id"]),
            "reference_id": str(reference_id),
            "negative_bank_version": str(negative_bank_version),
            "chrom": str(candidate_row["chrom"]),
            "start": int(candidate_row["start"]),
            "end": int(candidate_row["end"]),
            "state": str(candidate_row["state"]),
            "a_zscore": safe_float(candidate_row.get("a_zscore", np.nan), default=np.nan),
            "a_abs_zscore": safe_float(candidate_row.get("a_abs_zscore", np.nan), default=np.nan),
            "a_ratio": safe_float(candidate_row.get("a_ratio", np.nan), default=np.nan),
            "a_support_level": str(candidate_row.get("a_support_level", "")),
            "branch_b_event_id": first_text(event_row, ["event_id", "candidate_id"], default="") if event_row is not None else "",
            "branch_b_report_class": first_text(event_row, ["report_class"], default="") if event_row is not None else "",
            "branch_b_artifact_status": first_text(event_row, ["artifact_status"], default="") if event_row is not None else "",
            "branch_b_keep_event": first_numeric(event_row, ["keep_event"], default=np.nan) if event_row is not None else np.nan,
            "final_disposition": derive_disposition(event_row),
            "raw_amplitude": raw_amplitude,
            "corrected_amplitude": corrected_amplitude,
            "attenuation_ratio": attenuation_ratio,
            "same_direction_fraction": calculate_same_direction_fraction(candidate_row, event_bins),
            "flank_contrast": calculate_flank_contrast(candidate_row, event_bins, flanks),
            "hard_region_fraction": calculate_hard_region_fraction(event_row, event_bins),
            "sample_noise_mad": sample_noise_mad,
            "refmap_status": refmap_status,
            "calibration_null_status": calibration_null_status,
            "evidence_missing_reason": ";".join(reasons),
        }
        rows.append(payload)
    return pd.DataFrame(rows, columns=LEDGER_COLUMNS)


def summarize_evidence_ledger(sample_id, ledger_df):
    disposition_counts = (
        ledger_df["final_disposition"].fillna("UNKNOWN").astype(str).value_counts().sort_index().to_dict()
        if "final_disposition" in ledger_df.columns
        else {}
    )
    missing_count = (
        int(ledger_df["evidence_missing_reason"].fillna("").astype(str).ne("").sum())
        if "evidence_missing_reason" in ledger_df.columns
        else 0
    )
    return {
        "sample_id": str(sample_id),
        "candidate_count": int(len(ledger_df)),
        "disposition_counts": {str(key): int(value) for key, value in disposition_counts.items()},
        "missing_evidence_candidate_count": missing_count,
        "final_report_impact": "none_shadow_only",
    }


def main():
    args = parse_args()
    logger = setup_logger("branch_b_evidence_ledger", args.log or None)
    a_candidates = read_table(
        args.a_candidates,
        required_columns=["candidate_id", "sample_id", "chrom", "start", "end", "state"],
        empty_ok=True,
    )
    branch_b_events = read_optional_table(args.branch_b_events)
    bins_df = read_optional_table(args.input_bins)
    ledger = build_candidate_evidence_ledger(
        sample_id=args.sample_id,
        a_candidates=a_candidates,
        branch_b_events=branch_b_events,
        bins_df=bins_df,
        reference_id=args.reference_id,
        negative_bank_version=args.negative_bank_version,
    )
    write_table(args.output_ledger, ledger)
    write_json(args.output_summary, summarize_evidence_ledger(args.sample_id, ledger))
    logger.info("wrote evidence ledger rows=%d to %s", len(ledger), args.output_ledger)


if __name__ == "__main__":
    main()
