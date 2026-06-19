#!/biosoftware/miniconda/envs/snakemake_env/bin/python
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from pgta.core.logging import setup_logger
from pgta.predict.branch_b.common import read_table, write_json, write_table


MATCHED_NEGATIVE_COLUMNS = [
    "matched_negative_version",
    "matched_negative_feature",
    "matched_negative_query_abs",
    "matched_negative_background_status",
    "matched_negative_scope",
    "matched_negative_n",
    "matched_negative_abs_percentile",
    "matched_negative_abs_median",
    "matched_negative_abs_p95",
    "matched_negative_action",
]


def parse_args():
    parser = argparse.ArgumentParser(description="Add Phase 3 matched-negative percentile shadow evidence.")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--input-ledger", required=True)
    parser.add_argument("--negative-bank-labels", required=True)
    parser.add_argument("--background-ledger", action="append", default=[])
    parser.add_argument("--output-ledger", required=True)
    parser.add_argument("--output-summary", required=True)
    parser.add_argument("--version", default="branch_ab_v2")
    parser.add_argument("--feature-column", default="corrected_amplitude")
    parser.add_argument("--min-background", type=int, default=5)
    parser.add_argument("--similar-length-fold", type=float, default=2.0)
    parser.add_argument("--shadow-background-label", action="append", default=[])
    parser.add_argument("--log", default="")
    return parser.parse_args()


def safe_float(value, default=np.nan):
    try:
        result = float(value)
    except (TypeError, ValueError):
        return float(default)
    return result if np.isfinite(result) else float(default)


def read_optional_tables(paths):
    frames = []
    for raw_path in paths or []:
        if not raw_path:
            continue
        path = Path(raw_path)
        if not path.exists() or path.stat().st_size == 0:
            continue
        frames.append(read_table(path, empty_ok=True))
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True, sort=False)


def normalize_sample_ids(values):
    return {str(value).strip() for value in values if str(value).strip() and str(value).strip().lower() != "nan"}


def eligible_n0_samples(labels_df):
    if labels_df is None or labels_df.empty or "sample_id" not in labels_df.columns:
        return set()
    frame = labels_df.copy()
    eligible = pd.Series(False, index=frame.index)
    if "matched_negative_eligible" in frame.columns:
        eligible |= pd.to_numeric(frame["matched_negative_eligible"], errors="coerce").fillna(0).astype(float).gt(0)
    if "negative_bank_label" in frame.columns:
        eligible |= frame["negative_bank_label"].fillna("").astype(str).str.upper().eq("N0")
    return normalize_sample_ids(frame.loc[eligible, "sample_id"])


def eligible_label_samples(labels_df, allowed_labels):
    labels = {str(label).strip().upper() for label in allowed_labels or [] if str(label).strip()}
    if labels_df is None or labels_df.empty or "sample_id" not in labels_df.columns or not labels:
        return set()
    if "negative_bank_label" not in labels_df.columns:
        return set()
    frame = labels_df.copy()
    eligible = frame["negative_bank_label"].fillna("").astype(str).str.upper().isin(labels)
    return normalize_sample_ids(frame.loc[eligible, "sample_id"])


def normalize_background(background_df, labels_df, query_sample_id="", allowed_labels=None):
    background_samples = (
        eligible_label_samples(labels_df, allowed_labels)
        if allowed_labels
        else eligible_n0_samples(labels_df)
    )
    if background_df is None or background_df.empty or not background_samples:
        return pd.DataFrame(columns=list(background_df.columns) if background_df is not None else [])
    if "sample_id" not in background_df.columns:
        return pd.DataFrame(columns=list(background_df.columns))
    frame = background_df.copy()
    sample_ids = frame["sample_id"].fillna("").astype(str).str.strip()
    keep = sample_ids.isin(background_samples)
    if query_sample_id:
        keep &= ~sample_ids.eq(str(query_sample_id))
    frame = frame.loc[keep].copy()
    if frame.empty:
        return frame
    for column in ["start", "end"]:
        if column in frame.columns:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
    if "length_bp" not in frame.columns and {"start", "end"}.issubset(frame.columns):
        frame["length_bp"] = frame["end"] - frame["start"]
    return frame


def is_autosome(chrom):
    text = str(chrom).strip()
    if text.lower().startswith("chr"):
        text = text[3:]
    try:
        value = int(text)
    except ValueError:
        return False
    return 1 <= value <= 22


def similar_length_mask(frame, length_bp, fold):
    if frame.empty or "length_bp" not in frame.columns or not np.isfinite(length_bp) or length_bp <= 0:
        return pd.Series(False, index=frame.index)
    fold = max(float(fold), 1.0)
    lengths = pd.to_numeric(frame["length_bp"], errors="coerce")
    return lengths.between(length_bp / fold, length_bp * fold, inclusive="both")


def candidate_length(row):
    if "length_bp" in row.index:
        value = safe_float(row["length_bp"], default=np.nan)
        if np.isfinite(value) and value > 0:
            return value
    start = safe_float(row.get("start", np.nan), default=np.nan)
    end = safe_float(row.get("end", np.nan), default=np.nan)
    if np.isfinite(start) and np.isfinite(end) and end > start:
        return float(end - start)
    return np.nan


def select_background_rows(candidate, background_df, feature_column, min_background=5, similar_length_fold=2.0):
    empty = background_df.iloc[0:0].copy() if background_df is not None and not background_df.empty else pd.DataFrame()
    if background_df is None or background_df.empty:
        return "UNKNOWN_BACKGROUND", "UNKNOWN_BACKGROUND", empty
    required = {"chrom", "start", "end", feature_column}
    if not required.issubset(background_df.columns):
        return "UNKNOWN_BACKGROUND", "UNKNOWN_BACKGROUND", empty

    chrom = str(candidate.get("chrom", ""))
    start = safe_float(candidate.get("start", np.nan), default=np.nan)
    end = safe_float(candidate.get("end", np.nan), default=np.nan)
    length_bp = candidate_length(candidate)
    state = str(candidate.get("state", "")).strip().lower()

    state_mask = pd.Series(True, index=background_df.index)
    if state and "state" in background_df.columns:
        state_mask = background_df["state"].fillna("").astype(str).str.lower().eq(state)

    same_region = background_df[
        background_df["chrom"].astype(str).eq(chrom)
        & pd.to_numeric(background_df["start"], errors="coerce").eq(start)
        & pd.to_numeric(background_df["end"], errors="coerce").eq(end)
        & state_mask
    ].copy()
    same_region = same_region[pd.to_numeric(same_region[feature_column], errors="coerce").notna()]
    if len(same_region) >= min_background:
        return "OK", "same_region", same_region

    similar_length = similar_length_mask(background_df, length_bp, similar_length_fold)
    same_chrom = background_df[
        background_df["chrom"].astype(str).eq(chrom)
        & similar_length
        & state_mask
        & pd.to_numeric(background_df[feature_column], errors="coerce").notna()
    ].copy()
    if len(same_chrom) >= min_background:
        return "OK", "same_chrom_similar_length", same_chrom

    autosome = background_df[
        background_df["chrom"].map(is_autosome)
        & similar_length
        & state_mask
        & pd.to_numeric(background_df[feature_column], errors="coerce").notna()
    ].copy()
    if len(autosome) >= min_background:
        return "OK", "autosome_similar_length", autosome

    return "UNKNOWN_BACKGROUND", "UNKNOWN_BACKGROUND", empty


def percentile_leq(value, background_values):
    vector = pd.to_numeric(pd.Series(background_values), errors="coerce").abs().to_numpy(dtype=np.float64)
    vector = vector[np.isfinite(vector)]
    if vector.size == 0 or not np.isfinite(value):
        return np.nan
    return float(np.mean(vector <= abs(float(value))))


def matched_negative_payload(candidate, background_rows, status, scope, feature_column, version):
    feature_value = safe_float(candidate.get(feature_column, np.nan), default=np.nan)
    background_values = (
        pd.to_numeric(background_rows[feature_column], errors="coerce").abs().dropna().to_numpy(dtype=np.float64)
        if background_rows is not None and not background_rows.empty and feature_column in background_rows.columns
        else np.array([], dtype=np.float64)
    )
    if status not in {"OK", "SHADOW_BACKGROUND"} or background_values.size == 0 or not np.isfinite(feature_value):
        return {
            "matched_negative_version": str(version),
            "matched_negative_feature": str(feature_column),
            "matched_negative_query_abs": abs(feature_value) if np.isfinite(feature_value) else np.nan,
            "matched_negative_background_status": "UNKNOWN_BACKGROUND" if status == "UNKNOWN_BACKGROUND" else status,
            "matched_negative_scope": scope,
            "matched_negative_n": int(background_values.size),
            "matched_negative_abs_percentile": np.nan,
            "matched_negative_abs_median": np.nan,
            "matched_negative_abs_p95": np.nan,
            "matched_negative_action": "REVIEW_NO_CALL",
        }
    return {
        "matched_negative_version": str(version),
        "matched_negative_feature": str(feature_column),
        "matched_negative_query_abs": abs(feature_value),
        "matched_negative_background_status": str(status),
        "matched_negative_scope": scope,
        "matched_negative_n": int(background_values.size),
        "matched_negative_abs_percentile": percentile_leq(feature_value, background_values),
        "matched_negative_abs_median": float(np.median(background_values)),
        "matched_negative_abs_p95": float(np.percentile(background_values, 95)),
        "matched_negative_action": "BACKGROUND_SUPPORTED" if status == "OK" else "SHADOW_CONTEXT_ONLY",
    }


def build_matched_negative_percentiles(
    query_ledger,
    background_ledgers,
    negative_bank_labels,
    min_background=5,
    similar_length_fold=2.0,
    feature_column="corrected_amplitude",
    shadow_background_labels=None,
    version="branch_ab_v2",
):
    if query_ledger is None or query_ledger.empty:
        return pd.DataFrame(columns=list(query_ledger.columns if query_ledger is not None else []) + MATCHED_NEGATIVE_COLUMNS)
    if feature_column not in query_ledger.columns:
        raise ValueError(f"query ledger missing feature column: {feature_column}")
    query = query_ledger.copy()
    query_sample_id = str(query["sample_id"].iloc[0]) if "sample_id" in query.columns and len(query) else ""
    background = normalize_background(background_ledgers, negative_bank_labels, query_sample_id=query_sample_id)
    shadow_background = normalize_background(
        background_ledgers,
        negative_bank_labels,
        query_sample_id=query_sample_id,
        allowed_labels=shadow_background_labels,
    )

    payloads = []
    for _, candidate in query.iterrows():
        status, scope, background_rows = select_background_rows(
            candidate,
            background,
            feature_column,
            min_background=max(int(min_background), 1),
            similar_length_fold=similar_length_fold,
        )
        if status == "UNKNOWN_BACKGROUND" and shadow_background_labels:
            shadow_status, shadow_scope, shadow_rows = select_background_rows(
                candidate,
                shadow_background,
                feature_column,
                min_background=max(int(min_background), 1),
                similar_length_fold=similar_length_fold,
            )
            if shadow_status == "OK":
                status = "SHADOW_BACKGROUND"
                scope = shadow_scope
                background_rows = shadow_rows
        payloads.append(matched_negative_payload(candidate, background_rows, status, scope, feature_column, version))
    annotated = pd.concat([query.reset_index(drop=True), pd.DataFrame(payloads)], axis=1)
    return annotated


def summarize_matched_negative(sample_id, annotated_df, version="branch_ab_v2"):
    status_counts = (
        annotated_df["matched_negative_background_status"].fillna("UNKNOWN").astype(str).value_counts().sort_index().to_dict()
        if annotated_df is not None and "matched_negative_background_status" in annotated_df.columns
        else {}
    )
    action_counts = (
        annotated_df["matched_negative_action"].fillna("UNKNOWN").astype(str).value_counts().sort_index().to_dict()
        if annotated_df is not None and "matched_negative_action" in annotated_df.columns
        else {}
    )
    return {
        "sample_id": str(sample_id),
        "version": str(version),
        "candidate_count": int(len(annotated_df)) if annotated_df is not None else 0,
        "background_status_counts": {str(key): int(value) for key, value in status_counts.items()},
        "action_counts": {str(key): int(value) for key, value in action_counts.items()},
        "final_report_impact": "none_shadow_only",
    }


def main():
    args = parse_args()
    logger = setup_logger("branch_b_matched_negative", args.log or None)
    query = read_table(args.input_ledger, required_columns=["sample_id", "candidate_id", "chrom", "start", "end"], empty_ok=True)
    labels = read_table(args.negative_bank_labels, required_columns=["sample_id"], empty_ok=True)
    background = read_optional_tables(args.background_ledger)
    annotated = build_matched_negative_percentiles(
        query,
        background,
        labels,
        min_background=args.min_background,
        similar_length_fold=args.similar_length_fold,
        feature_column=args.feature_column,
        shadow_background_labels=set(args.shadow_background_label or []),
        version=args.version,
    )
    write_table(args.output_ledger, annotated)
    write_json(args.output_summary, summarize_matched_negative(args.sample_id, annotated, version=args.version))
    logger.info("wrote matched-negative rows=%d to %s", len(annotated), args.output_ledger)


if __name__ == "__main__":
    main()
