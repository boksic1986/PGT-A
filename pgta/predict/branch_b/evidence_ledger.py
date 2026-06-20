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
    "cnvpro_like_gc_rc_background_status",
    "dynamic_reference_status",
    "matched_negative_source",
    "matched_negative_percentile",
    "copy_number_estimate",
    "sex_adjusted_copy_number",
    "mosaic_fraction_proxy",
    "mosaic_proxy_status",
    "event_arm_class",
    "event_par_class",
    "crosses_centromere",
    "crosses_par_boundary",
    "whole_chromosome_fraction",
    "cnvpro_large_segment_tier",
    "waviness",
    "sample_noise_status",
    "cnvpro_like_evidence_status",
    "evidence_missing_reason",
    "sample",
    "branch_b_direction_support",
    "copy_number_like_amplitude",
    "mosaic_proxy",
    "loh_evidence",
    "upd_evidence",
    "background_source",
    "background_status",
    "region_risk_context",
    "sample_noise_context",
    "cnvpro_consistency_status",
    "disposition",
    "disposition_reason",
    "report_impact",
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


def classify_cnvpro_like_gc_rc_background_status(raw_amplitude, corrected_amplitude, event_bins):
    if event_bins.empty:
        return "UNKNOWN"
    if np.isfinite(raw_amplitude) and np.isfinite(corrected_amplitude):
        return "GC_RC_AVAILABLE"
    if np.isfinite(raw_amplitude):
        return "RAW_ONLY"
    return "UNKNOWN"


def classify_dynamic_reference_status(refmap_status):
    status = str(refmap_status or "").strip().upper()
    if status == "OK":
        return "OK"
    if status == "UNKNOWN":
        return "UNKNOWN"
    if status == "SAME_CHROM_REF_LEAKAGE":
        return "INVALID_REF_LEAKAGE"
    if status.startswith("LOW_REF") or status.startswith("LOW_REFBIN"):
        return "REVIEW"
    return "REVIEW"


def estimate_copy_number(candidate_row, event_row):
    if event_row is not None:
        value = first_numeric(event_row, ["copy_number_estimate", "sex_adjusted_copy_number"], default=np.nan)
        if np.isfinite(value):
            return value
    ratio = safe_float(candidate_row.get("a_ratio", np.nan), default=np.nan)
    if not np.isfinite(ratio):
        return np.nan
    return float(2.0 * (1.0 + ratio))


def estimate_sex_adjusted_copy_number(candidate_row, event_row, copy_number):
    if event_row is not None:
        value = first_numeric(event_row, ["sex_adjusted_copy_number"], default=np.nan)
        if np.isfinite(value):
            return value
    return copy_number


def estimate_mosaic_fraction(candidate_row, copy_number):
    if not np.isfinite(copy_number):
        return np.nan, "UNKNOWN"
    direction = state_direction(candidate_row["state"])
    if direction > 0.0:
        value = copy_number - 2.0
    elif direction < 0.0:
        value = 2.0 - copy_number
    else:
        return np.nan, "UNKNOWN"
    return float(np.clip(value, 0.0, 1.0)), "AVAILABLE"


def classify_event_arm(candidate_row, event_row):
    if event_row is not None:
        value = first_text(event_row, ["event_arm_class", "arm_class"], default="")
        if value:
            return value
    chrom = str(candidate_row["chrom"])
    if chrom in {"chrX", "chrY", "X", "Y"}:
        return "sex_chromosome"
    return "autosome"


def _boolean_column_fraction(frame, columns):
    if frame.empty:
        return np.nan
    for column in columns:
        if column not in frame.columns:
            continue
        values = pd.to_numeric(frame[column], errors="coerce")
        if values.notna().any():
            return float(values.fillna(0.0).gt(0.0).mean())
        text = frame[column].fillna("").astype(str).str.lower()
        if text.ne("").any():
            return float(text.isin({"1", "true", "yes", "par"}).mean())
    return np.nan


def classify_event_par(candidate_row, event_bins):
    chrom = str(candidate_row["chrom"])
    if chrom not in {"chrX", "chrY", "X", "Y"}:
        return "autosome"
    par_fraction = _boolean_column_fraction(
        event_bins,
        ["is_PAR", "is_par", "is_par_region", "par_overlap_fraction"],
    )
    if not np.isfinite(par_fraction):
        return "UNKNOWN"
    if par_fraction >= 0.95:
        return "par_only"
    if par_fraction > 0.0:
        return "mixed_par_nonpar"
    return "nonpar"


def calculate_crosses_centromere(event_row, event_bins):
    if event_row is not None:
        value = first_numeric(event_row, ["crosses_centromere", "cnvseq_crosses_gap_or_centromere"], default=np.nan)
        if np.isfinite(value):
            return int(value > 0.0)
    if event_bins.empty:
        return 0
    for column in ["is_near_centromere", "near_centromere", "gap_centromere_telomere_overlap_fraction"]:
        if column not in event_bins.columns:
            continue
        values = pd.to_numeric(event_bins[column], errors="coerce")
        if values.notna().any() and bool(values.fillna(0.0).gt(0.0).any()):
            return 1
    return 0


def calculate_whole_chromosome_fraction(candidate_row, event_bins, bins_df):
    if event_bins.empty or bins_df is None or bins_df.empty or "chrom" not in bins_df.columns:
        return np.nan
    chrom = str(candidate_row["chrom"])
    chrom_bins = bins_df[bins_df["chrom"].astype(str).eq(chrom)]
    if chrom_bins.empty:
        return np.nan
    return float(len(event_bins) / len(chrom_bins))


def classify_large_segment_tier(candidate_row, whole_chromosome_fraction):
    length = max(int(candidate_row["end"]) - int(candidate_row["start"]), 0)
    if length >= 20_000_000 and np.isfinite(whole_chromosome_fraction) and whole_chromosome_fraction >= 0.90:
        return "whole_chromosome"
    if length >= 10_000_000:
        return "large_ge10mb"
    if length >= 4_000_000:
        return "large_ge4mb"
    if length >= 2_000_000:
        return "reportable_ge2mb"
    if length >= 1_000_000:
        return "review_ge1mb"
    return "subreportable_lt1mb"


def classify_sample_noise_status(waviness):
    if not np.isfinite(waviness):
        return "UNKNOWN"
    if waviness <= 1.0:
        return "OK"
    if waviness <= 2.0:
        return "MODERATE_WAVINESS"
    return "HIGH_WAVINESS"


def classify_cnvpro_like_evidence_status(event_bins):
    if event_bins.empty:
        return "UNKNOWN"
    return "SHADOW_EVIDENCE_ONLY"


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


def has_branch_a_positive_support(candidate_row):
    support_level = str(candidate_row.get("a_support_level", "")).strip().lower()
    a_abs_zscore = safe_float(candidate_row.get("a_abs_zscore", np.nan), default=np.nan)
    return "strong" in support_level or "sensitive" in support_level or (np.isfinite(a_abs_zscore) and a_abs_zscore >= 5.0)


def classify_branch_b_direction_support(candidate_row, same_direction_fraction, raw_amplitude, corrected_amplitude):
    if np.isfinite(same_direction_fraction):
        if same_direction_fraction >= 0.70:
            return "SUPPORTED"
        if same_direction_fraction >= 0.50:
            return "BORDERLINE_SUPPORTED"
        return "DISCORDANT"
    direction = state_direction(candidate_row.get("state", ""))
    for value in [corrected_amplitude, raw_amplitude]:
        if np.isfinite(value) and direction != 0.0:
            return "SUPPORTED" if value * direction > 0.0 else "DISCORDANT"
    return "UNKNOWN"


def p3_copy_number_like_amplitude(corrected_amplitude, raw_amplitude):
    if np.isfinite(corrected_amplitude):
        return corrected_amplitude
    if np.isfinite(raw_amplitude):
        return raw_amplitude
    return np.nan


def p3_mosaic_proxy(mosaic_fraction_proxy, mosaic_proxy_status):
    if str(mosaic_proxy_status).upper() == "AVAILABLE" and np.isfinite(mosaic_fraction_proxy):
        return mosaic_fraction_proxy
    return "not_available"


def p3_background_status(source):
    text = str(source or "").strip()
    return text if text else "UNKNOWN_BACKGROUND"


def p3_region_risk_context(refmap_status, hard_region_fraction):
    status = str(refmap_status or "").strip().upper()
    if status and status not in {"OK", "UNKNOWN"}:
        return status
    if np.isfinite(hard_region_fraction) and hard_region_fraction >= 0.50:
        return "HIGH_REGION_RISK"
    if status == "OK":
        return "OK"
    return "UNKNOWN"


def derive_p3_disposition(candidate_row, final_disposition, refmap_status):
    if str(refmap_status or "").strip().upper() == "SAME_CHROM_REF_LEAKAGE":
        return "NO_CALL_CONTRACT_RISK"
    # P3 is evidence refinement, not final confirmation/suppression. Legacy
    # dispositions are retained in final_disposition but do not become P3 hard
    # report decisions.
    return "REVIEW_REQUIRED"


def derive_p3_disposition_reason(candidate_row, final_disposition, evidence_missing_reason, refmap_status):
    reasons = []
    if evidence_missing_reason:
        reasons.extend([item for item in str(evidence_missing_reason).split(";") if item])
    legacy = str(final_disposition or "").strip().upper()
    if legacy:
        append_reason(reasons, f"legacy_final_disposition={legacy}")
    if legacy == "LIKELY_ARTIFACT" and has_branch_a_positive_support(candidate_row):
        append_reason(reasons, "legacy_artifact_not_hard_suppressed")
    if str(refmap_status or "").strip().upper() == "SAME_CHROM_REF_LEAKAGE":
        append_reason(reasons, "same_chrom_ref_leakage_contract_risk")
    return ";".join(reasons) if reasons else "p3_review_required"


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
        cnvpro_like_gc_rc_background_status = classify_cnvpro_like_gc_rc_background_status(
            raw_amplitude,
            corrected_amplitude,
            event_bins,
        )
        dynamic_reference_status = classify_dynamic_reference_status(refmap_status)
        copy_number_estimate = estimate_copy_number(candidate_row, event_row)
        sex_adjusted_copy_number = estimate_sex_adjusted_copy_number(candidate_row, event_row, copy_number_estimate)
        mosaic_fraction_proxy, mosaic_proxy_status = estimate_mosaic_fraction(candidate_row, sex_adjusted_copy_number)
        event_par_class = classify_event_par(candidate_row, event_bins)
        crosses_centromere = calculate_crosses_centromere(event_row, event_bins)
        crosses_par_boundary = int(event_par_class == "mixed_par_nonpar")
        whole_chromosome_fraction = calculate_whole_chromosome_fraction(candidate_row, event_bins, bins_df)
        cnvpro_large_segment_tier = classify_large_segment_tier(candidate_row, whole_chromosome_fraction)
        waviness = sample_noise_mad
        sample_noise_status = classify_sample_noise_status(waviness)
        cnvpro_like_evidence_status = classify_cnvpro_like_evidence_status(event_bins)
        same_direction_fraction = calculate_same_direction_fraction(candidate_row, event_bins)
        flank_contrast = calculate_flank_contrast(candidate_row, event_bins, flanks)
        hard_region_fraction = calculate_hard_region_fraction(event_row, event_bins)
        reasons = []
        if event_row is None:
            append_reason(reasons, "branch_b_event_missing")
        if refmap_status == "UNKNOWN":
            append_reason(reasons, "refmap_missing")
        if calibration_null_status == "UNKNOWN":
            append_reason(reasons, "calibration_null_missing")
        if event_bins.empty:
            append_reason(reasons, "event_bins_missing")
        matched_negative_source = "UNKNOWN_BACKGROUND"
        final_disposition = derive_disposition(event_row)
        evidence_missing_reason = ";".join(reasons)
        copy_number_like_amplitude = p3_copy_number_like_amplitude(corrected_amplitude, raw_amplitude)
        mosaic_proxy = p3_mosaic_proxy(mosaic_fraction_proxy, mosaic_proxy_status)
        background_status = p3_background_status(matched_negative_source)
        region_risk_context = p3_region_risk_context(refmap_status, hard_region_fraction)
        disposition = derive_p3_disposition(candidate_row, final_disposition, refmap_status)
        disposition_reason = derive_p3_disposition_reason(
            candidate_row,
            final_disposition,
            evidence_missing_reason,
            refmap_status,
        )

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
            "final_disposition": final_disposition,
            "raw_amplitude": raw_amplitude,
            "corrected_amplitude": corrected_amplitude,
            "attenuation_ratio": attenuation_ratio,
            "same_direction_fraction": same_direction_fraction,
            "flank_contrast": flank_contrast,
            "hard_region_fraction": hard_region_fraction,
            "sample_noise_mad": sample_noise_mad,
            "refmap_status": refmap_status,
            "calibration_null_status": calibration_null_status,
            "cnvpro_like_gc_rc_background_status": cnvpro_like_gc_rc_background_status,
            "dynamic_reference_status": dynamic_reference_status,
            "matched_negative_source": matched_negative_source,
            "matched_negative_percentile": np.nan,
            "copy_number_estimate": copy_number_estimate,
            "sex_adjusted_copy_number": sex_adjusted_copy_number,
            "mosaic_fraction_proxy": mosaic_fraction_proxy,
            "mosaic_proxy_status": mosaic_proxy_status,
            "event_arm_class": classify_event_arm(candidate_row, event_row),
            "event_par_class": event_par_class,
            "crosses_centromere": crosses_centromere,
            "crosses_par_boundary": crosses_par_boundary,
            "whole_chromosome_fraction": whole_chromosome_fraction,
            "cnvpro_large_segment_tier": cnvpro_large_segment_tier,
            "waviness": waviness,
            "sample_noise_status": sample_noise_status,
            "cnvpro_like_evidence_status": cnvpro_like_evidence_status,
            "evidence_missing_reason": evidence_missing_reason,
            "sample": str(sample_id or candidate_row["sample_id"]),
            "branch_b_direction_support": classify_branch_b_direction_support(
                candidate_row,
                same_direction_fraction,
                raw_amplitude,
                corrected_amplitude,
            ),
            "copy_number_like_amplitude": copy_number_like_amplitude,
            "mosaic_proxy": mosaic_proxy,
            "loh_evidence": "not_available",
            "upd_evidence": "not_available",
            "background_source": matched_negative_source,
            "background_status": background_status,
            "region_risk_context": region_risk_context,
            "sample_noise_context": sample_noise_status,
            "cnvpro_consistency_status": cnvpro_like_evidence_status,
            "disposition": disposition,
            "disposition_reason": disposition_reason,
            "report_impact": "none_shadow_only",
        }
        rows.append(payload)
    return pd.DataFrame(rows, columns=LEDGER_COLUMNS)


def summarize_evidence_ledger(sample_id, ledger_df):
    def _value_counts(column):
        if column not in ledger_df.columns:
            return {}
        return ledger_df[column].fillna("UNKNOWN").astype(str).value_counts().sort_index().to_dict()

    disposition_counts = (
        ledger_df["final_disposition"].fillna("UNKNOWN").astype(str).value_counts().sort_index().to_dict()
        if "final_disposition" in ledger_df.columns
        else {}
    )
    p3_disposition_counts = (
        ledger_df["disposition"].fillna("UNKNOWN").astype(str).value_counts().sort_index().to_dict()
        if "disposition" in ledger_df.columns
        else {}
    )
    missing_count = (
        int(ledger_df["evidence_missing_reason"].fillna("").astype(str).ne("").sum())
        if "evidence_missing_reason" in ledger_df.columns
        else 0
    )
    review_burden_count = (
        int(ledger_df["disposition"].fillna("").astype(str).str.contains("REVIEW", case=False, regex=False).sum())
        if "disposition" in ledger_df.columns
        else 0
    )
    return {
        "sample_id": str(sample_id),
        "candidate_count": int(len(ledger_df)),
        "disposition_counts": {str(key): int(value) for key, value in disposition_counts.items()},
        "p3_disposition_counts": {str(key): int(value) for key, value in p3_disposition_counts.items()},
        "background_source_counts": {str(key): int(value) for key, value in _value_counts("background_source").items()},
        "background_status_counts": {str(key): int(value) for key, value in _value_counts("background_status").items()},
        "review_burden_count": review_burden_count,
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
