#!/biosoftware/miniconda/envs/snakemake_env/bin/python
from __future__ import annotations

import argparse
import math
from pathlib import Path

import pandas as pd

from pgta.predict.branch_b.common import interval_overlap, read_table, write_json, write_table


LOWRES_VERSION = "branch_b_v2_lowres_evidence_v1"
LOWRES_2MB_INFORMATIVE_MIN_BP = 3_000_000
LOWRES_3MB_INFORMATIVE_MIN_BP = 4_000_000
LOWRES_MIN_OVERLAP_FRACTION = 0.50


def parse_args():
    parser = argparse.ArgumentParser(description="Annotate Branch A candidates with 2Mb/3Mb low-resolution evidence.")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--input-ledger", required=True)
    parser.add_argument("--output-ledger", required=True)
    parser.add_argument("--output-summary", required=True)
    parser.add_argument("--lowres-2mb-events", default="")
    parser.add_argument("--lowres-3mb-events", default="")
    parser.add_argument("--ref-stability-events", default="")
    return parser.parse_args()


def _clean_chrom(value):
    text = str(value if value is not None else "").strip()
    if not text:
        return ""
    return text if text.lower().startswith("chr") else f"chr{text}"


def _safe_float(value, default=math.nan):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return float(default)
    return number if math.isfinite(number) else float(default)


def _event_length(row):
    start = _safe_float(row.get("start", math.nan))
    end = _safe_float(row.get("end", math.nan))
    if not math.isfinite(start) or not math.isfinite(end) or end <= start:
        return math.nan
    return float(end - start)


def _is_sex_chrom(row):
    chrom = _clean_chrom(row.get("chrom", "")).lower().removeprefix("chr")
    return chrom in {"x", "y"}


def _candidate_id(row, index):
    value = row.get("candidate_id", "")
    try:
        if pd.isna(value):
            value = ""
    except (TypeError, ValueError):
        pass
    text = str(value).strip()
    return text if text else f"candidate_{index}"


def _score_from_event(row):
    for column in [
        "a_abs_zscore",
        "segment_abs_max_robust_z",
        "event_abs_z",
        "abs_zscore",
        "zscore",
        "score",
    ]:
        if column not in row:
            continue
        value = _safe_float(row.get(column, math.nan))
        if math.isfinite(value):
            return abs(value)
    return math.nan


def _normalize_events(events, sample_id):
    if events is None:
        return None
    frame = events.copy()
    if frame.empty:
        return frame
    for column in ["chrom", "start", "end", "state"]:
        if column not in frame.columns:
            frame[column] = ""
    if "sample_id" in frame.columns:
        frame = frame.loc[frame["sample_id"].astype(str).eq(str(sample_id))].copy()
    frame["chrom"] = frame["chrom"].map(_clean_chrom)
    frame["start"] = pd.to_numeric(frame["start"], errors="coerce")
    frame["end"] = pd.to_numeric(frame["end"], errors="coerce")
    frame["state"] = frame["state"].fillna("").astype(str).str.lower()
    return frame.loc[frame["start"].notna() & frame["end"].notna() & (frame["end"] > frame["start"])].copy()


def _annotate_resolution(candidates, events, sample_id, prefix, informative_min_bp):
    frame = candidates.copy()
    events = _normalize_events(events, sample_id)
    configured = events is not None
    rows = []

    for index, row in frame.iterrows():
        candidate_length = _event_length(row)
        base = {
            f"{prefix}_support_label": "LOWRES_NOT_CONFIGURED",
            f"{prefix}_same_direction": 0,
            f"{prefix}_overlap_fraction": 0.0,
            f"{prefix}_z_or_score": math.nan,
        }
        if _is_sex_chrom(row):
            base[f"{prefix}_support_label"] = "LOWRES_NOT_APPLICABLE_SEX_CHROM"
            rows.append(base)
            continue
        if not configured:
            rows.append(base)
            continue
        if not math.isfinite(candidate_length) or candidate_length < informative_min_bp:
            base[f"{prefix}_support_label"] = "LOWRES_NOT_INFORMATIVE_SHORT_OR_BOUNDARY_EVENT"
            rows.append(base)
            continue

        chrom = _clean_chrom(row.get("chrom", ""))
        state = str(row.get("state", "")).strip().lower()
        start = int(_safe_float(row.get("start", 0), default=0))
        end = int(_safe_float(row.get("end", 0), default=0))
        overlaps = []
        if events is not None and not events.empty:
            subset = events.loc[events["chrom"].eq(chrom) & events["state"].eq(state)]
            for _, event in subset.iterrows():
                overlap_bp = interval_overlap(start, end, int(event["start"]), int(event["end"]))
                if overlap_bp <= 0:
                    continue
                overlaps.append((float(overlap_bp) / candidate_length, _score_from_event(event)))

        if overlaps:
            best_overlap = max(value[0] for value in overlaps)
            best_score = max((value[1] for value in overlaps if math.isfinite(value[1])), default=math.nan)
            base[f"{prefix}_overlap_fraction"] = round(float(best_overlap), 6)
            base[f"{prefix}_z_or_score"] = best_score
            if best_overlap >= LOWRES_MIN_OVERLAP_FRACTION:
                base[f"{prefix}_support_label"] = "LOWRES_SAME_DIRECTION_SUPPORT"
                base[f"{prefix}_same_direction"] = 1
            else:
                base[f"{prefix}_support_label"] = "LOWRES_PARTIAL_SAME_DIRECTION_CONTEXT"
        else:
            base[f"{prefix}_support_label"] = "LOWRES_NO_SAME_DIRECTION_SUPPORT"
        rows.append(base)

    annotations = pd.DataFrame(rows, index=frame.index)
    return pd.concat([frame, annotations], axis=1)


def _consensus_label(row):
    if _is_sex_chrom(row):
        return "LOWRES_NOT_APPLICABLE_SEX_CHROM"
    length = _event_length(row)
    support_2mb = row.get("lowres_2mb_support_label") == "LOWRES_SAME_DIRECTION_SUPPORT"
    support_3mb = row.get("lowres_3mb_support_label") == "LOWRES_SAME_DIRECTION_SUPPORT"
    no_2mb = row.get("lowres_2mb_support_label") == "LOWRES_NO_SAME_DIRECTION_SUPPORT"
    no_3mb = row.get("lowres_3mb_support_label") == "LOWRES_NO_SAME_DIRECTION_SUPPORT"
    if support_2mb and support_3mb:
        return "LOWRES_2MB_3MB_SAME_DIRECTION_SUPPORT"
    if support_2mb and math.isfinite(length) and length < LOWRES_3MB_INFORMATIVE_MIN_BP:
        return "LOWRES_2MB_SUPPORT_FOR_3_4MB_EVENT"
    if support_2mb:
        return "LOWRES_2MB_SAME_DIRECTION_SUPPORT"
    if support_3mb:
        return "LOWRES_3MB_SAME_DIRECTION_SUPPORT"
    if math.isfinite(length) and length < LOWRES_2MB_INFORMATIVE_MIN_BP:
        return "LOWRES_CONTEXT_ONLY_SHORT_OR_BOUNDARY_EVENT"
    if no_2mb and row.get("lowres_3mb_support_label") in {
        "LOWRES_NO_SAME_DIRECTION_SUPPORT",
        "LOWRES_NOT_INFORMATIVE_SHORT_OR_BOUNDARY_EVENT",
        "LOWRES_NOT_CONFIGURED",
    }:
        return "LOWRES_NO_SUPPORT_INFORMATIVE_BUT_NOT_FILTER"
    if row.get("lowres_2mb_support_label") == "LOWRES_NOT_CONFIGURED" and row.get("lowres_3mb_support_label") == "LOWRES_NOT_CONFIGURED":
        return "LOWRES_NOT_CONFIGURED"
    return "LOWRES_CONTEXT_ONLY"


def _merge_ref_stability(candidates, ref_stability_events):
    if ref_stability_events is None or ref_stability_events.empty:
        return candidates
    frame = candidates.copy()
    stability = ref_stability_events.copy()
    if "candidate_id" not in stability.columns:
        return frame
    keep_columns = [
        column
        for column in [
            "candidate_id",
            "event_ref_mad_median",
            "event_ref_mad_p90",
            "high_ref_mad_bin_fraction",
            "ref_stability_context",
        ]
        if column in stability.columns
    ]
    if len(keep_columns) <= 1:
        return frame
    return frame.merge(stability[keep_columns], on="candidate_id", how="left", suffixes=("", "_ref_stability"))


def summarize_lowres_evidence(sample_id, annotated):
    def counts(column):
        if column not in annotated.columns:
            return {}
        return {
            str(key): int(value)
            for key, value in annotated[column].fillna("UNKNOWN").astype(str).value_counts().sort_index().to_dict().items()
        }

    return {
        "sample_id": str(sample_id),
        "version": LOWRES_VERSION,
        "candidate_count": int(len(annotated)),
        "lowres_2mb_support_label_counts": counts("lowres_2mb_support_label"),
        "lowres_3mb_support_label_counts": counts("lowres_3mb_support_label"),
        "lowres_consensus_label_counts": counts("lowres_consensus_label"),
        "ref_stability_context_counts": counts("ref_stability_context"),
        "final_report_impact": "development_review_only",
        "lowres_absence_is_filter_evidence": False,
    }


def annotate_lowres_evidence(
    candidates,
    sample_id,
    lowres_2mb_events=None,
    lowres_3mb_events=None,
    ref_stability_events=None,
):
    frame = candidates.copy()
    if frame.empty:
        annotated = frame.copy()
        for column in [
            "lowres_2mb_support_label",
            "lowres_2mb_same_direction",
            "lowres_2mb_overlap_fraction",
            "lowres_2mb_z_or_score",
            "lowres_3mb_support_label",
            "lowres_3mb_same_direction",
            "lowres_3mb_overlap_fraction",
            "lowres_3mb_z_or_score",
            "lowres_consensus_label",
        ]:
            annotated[column] = []
        return annotated, summarize_lowres_evidence(sample_id, annotated)

    frame = _annotate_resolution(
        frame,
        lowres_2mb_events,
        sample_id,
        "lowres_2mb",
        LOWRES_2MB_INFORMATIVE_MIN_BP,
    )
    frame = _annotate_resolution(
        frame,
        lowres_3mb_events,
        sample_id,
        "lowres_3mb",
        LOWRES_3MB_INFORMATIVE_MIN_BP,
    )
    frame["lowres_consensus_label"] = frame.apply(_consensus_label, axis=1)
    frame = _merge_ref_stability(frame, ref_stability_events)
    return frame, summarize_lowres_evidence(sample_id, frame)


def _read_optional_table(path_value):
    if not path_value:
        return None
    path = Path(path_value)
    if not path.exists() or path.stat().st_size == 0:
        return None
    return read_table(path, empty_ok=True)


def main():
    args = parse_args()
    candidates = read_table(args.input_ledger, required_columns=["candidate_id"], empty_ok=True)
    lowres_2mb = _read_optional_table(args.lowres_2mb_events)
    lowres_3mb = _read_optional_table(args.lowres_3mb_events)
    ref_stability = _read_optional_table(args.ref_stability_events)
    annotated, summary = annotate_lowres_evidence(
        candidates,
        sample_id=args.sample_id,
        lowres_2mb_events=lowres_2mb,
        lowres_3mb_events=lowres_3mb,
        ref_stability_events=ref_stability,
    )
    write_table(args.output_ledger, annotated)
    write_json(args.output_summary, summary)


if __name__ == "__main__":
    main()
