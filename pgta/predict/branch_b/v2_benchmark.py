#!/biosoftware/miniconda/envs/snakemake_env/bin/python
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import pandas as pd

from pgta.predict.truth import (
    filter_truth_to_samples,
    normalize_chrom,
    normalize_state,
    overlap_fraction,
)


TRUTH_METRIC_COLUMNS = [
    "sample_id",
    "truth_index",
    "chrom",
    "start",
    "end",
    "state",
    "v2_overlap_count",
    "v2_preserved_count",
    "v2_positive_support_count",
    "v2_background_informative_count",
    "v2_hard_suppressed_count",
    "top_candidate_id",
    "top_a_abs_zscore",
    "top_v2_candidate_class",
    "top_v2_evidence_tier",
    "top_v2_classifier_action",
    "top_v2_disposition",
    "top_v2_filter_action",
    "top_v2_length_tier",
    "top_v2_clean_support_label",
    "top_v2_gc_rc_context_label",
    "top_v2_b_signal_context_label",
    "top_v2_burden_reduction_tier",
    "top_v2_burden_reduction_action",
    "top_v2_burden_reduction_reason",
    "top_v2_burden_evidence_tags",
    "top_v2_report_layer_class",
    "top_v2_report_visibility",
    "top_v2_report_filter_reason",
    "top_v2_report_filter_rule_tags",
    "top_attenuation_ratio",
]


SAMPLE_SUMMARY_COLUMNS = [
    "sample_id",
    "truth_event_count",
    "truth_preserved_count",
    "FN_count",
    "candidate_count",
    "v2_positive_support_count",
    "v2_background_compatible_count",
    "v2_technical_review_count",
    "v2_contract_risk_count",
    "v2_hard_suppressed_count",
    "v2_filter_suppressed_count",
    "v2_report_candidate_burden_count",
    "v2_review_candidate_burden_count",
    "v2_background_unknown_review_burden_count",
    "v2_technical_risk_burden_count",
    "v2_branch_s_review_burden_count",
    "v2_report_event_count",
    "sample_report_burden_flag",
    "sample_report_burden_reason",
    "v2_internal_review_event_count",
    "v2_filtered_event_count",
    "v2_branch_s_event_count",
]

SAMPLE_REPORT_BURDEN_THRESHOLD = 3


def parse_args():
    parser = argparse.ArgumentParser(
        description="Benchmark Branch B V2-only candidate classifications against locked truth."
    )
    parser.add_argument("--sample-id", action="append", default=[])
    parser.add_argument("--classification-tsv", action="append", default=[])
    parser.add_argument("--truth-tsv", default="")
    parser.add_argument("--reference-id", default="UNKNOWN_REFERENCE")
    parser.add_argument("--output-truth-metrics", required=True)
    parser.add_argument("--output-sample-summary", required=True)
    parser.add_argument("--output-filtered-events", default="")
    parser.add_argument("--output-filtered-events-json", default="")
    parser.add_argument("--output-report-events", default="")
    parser.add_argument("--output-report-events-json", default="")
    parser.add_argument("--output-summary", required=True)
    return parser.parse_args()


def _read_table(path_value: str, columns: list[str] | None = None) -> pd.DataFrame:
    if not path_value:
        return pd.DataFrame(columns=columns or [])
    path = Path(path_value)
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame(columns=columns or [])
    return pd.read_csv(path, sep="\t")


def _sample_from_classification_path(path_value: str) -> str:
    name = Path(path_value).name
    suffix = ".candidate_classification.tsv"
    return name[: -len(suffix)] if name.endswith(suffix) else Path(path_value).stem


def load_classifications(paths: list[str]) -> pd.DataFrame:
    frames = []
    for path_value in paths:
        df = _read_table(path_value)
        if df.empty:
            sample_id = _sample_from_classification_path(path_value)
            frames.append(pd.DataFrame([{"sample_id": sample_id}]).iloc[0:0])
            continue
        frames.append(df)
    if not frames:
        return pd.DataFrame(columns=["sample_id", "candidate_id", "chrom", "start", "end", "state"])
    out = pd.concat(frames, ignore_index=True, sort=False)
    if "sample" in out.columns and "sample_id" not in out.columns:
        out = out.rename(columns={"sample": "sample_id"})
    for column in ["sample_id", "candidate_id", "chrom", "state"]:
        if column not in out.columns:
            out[column] = ""
    for column in ["start", "end", "a_abs_zscore"]:
        if column not in out.columns:
            out[column] = math.nan
    out["sample_id"] = out["sample_id"].astype(str)
    out["chrom"] = out["chrom"].map(normalize_chrom)
    out["state"] = out["state"].map(normalize_state)
    out["start"] = pd.to_numeric(out["start"], errors="coerce")
    out["end"] = pd.to_numeric(out["end"], errors="coerce")
    out["a_abs_zscore"] = pd.to_numeric(out["a_abs_zscore"], errors="coerce")
    return out


def load_truth(truth_tsv: str, sample_ids: list[str]) -> pd.DataFrame:
    if not truth_tsv:
        return pd.DataFrame(columns=["sample_id", "chrom", "start", "end", "state"])
    truth = _read_table(truth_tsv)
    if truth.empty:
        return pd.DataFrame(columns=["sample_id", "chrom", "start", "end", "state"])
    if "sample" in truth.columns and "sample_id" not in truth.columns:
        truth = truth.rename(columns={"sample": "sample_id"})
    if "expected_state" in truth.columns and "state" not in truth.columns:
        truth = truth.rename(columns={"expected_state": "state"})
    if "svtype" in truth.columns and "state" not in truth.columns:
        truth = truth.rename(columns={"svtype": "state"})
    missing = sorted({"sample_id", "chrom", "state"}.difference(truth.columns))
    if missing:
        raise ValueError(f"truth table missing columns: {','.join(missing)}")
    for column in ["start", "end"]:
        if column not in truth.columns:
            truth[column] = math.nan
    truth = filter_truth_to_samples(truth, sample_ids)
    truth["sample_id"] = truth["sample_id"].astype(str)
    truth["chrom"] = truth["chrom"].map(normalize_chrom)
    truth["state"] = truth["state"].map(normalize_state)
    truth["start"] = pd.to_numeric(truth["start"], errors="coerce")
    truth["end"] = pd.to_numeric(truth["end"], errors="coerce")
    return truth.reset_index(drop=True)


def hard_suppression_mask(frame: pd.DataFrame) -> pd.Series:
    if frame.empty:
        return pd.Series(False, index=frame.index, dtype=bool)
    action = frame.get("v2_classifier_action", pd.Series("", index=frame.index)).fillna("").astype(str).str.upper()
    filter_action = frame.get("v2_filter_action", pd.Series("", index=frame.index)).fillna("").astype(str).str.lower()
    klass = frame.get("v2_candidate_class", pd.Series("", index=frame.index)).fillna("").astype(str).str.upper()
    hard_action = (
        action.str.startswith("V2_SUPPRESS_")
        | action.str.startswith("SUPPRESS_")
        | action.str.startswith("V2_HARD_SUPPRESS")
        | action.str.startswith("HARD_SUPPRESS")
        | action.str.contains("ARTIFACT", regex=False)
    )
    hard_filter_action = filter_action.eq("suppress_workflow_contract_risk") | filter_action.str.startswith("suppress_")
    return hard_action | hard_filter_action | klass.str.contains(
        "ARTIFACT", regex=False
    )


def report_layer_filtered_mask(frame: pd.DataFrame) -> pd.Series:
    if frame.empty:
        return pd.Series(False, index=frame.index, dtype=bool)
    report_class = frame.get("v2_report_layer_class", pd.Series("", index=frame.index)).fillna("").astype(str)
    visibility = frame.get("v2_report_visibility", pd.Series("", index=frame.index)).fillna("").astype(str)
    return report_class.eq("filtered_event") | visibility.eq("audit_only")


def positive_support_mask(frame: pd.DataFrame) -> pd.Series:
    if frame.empty:
        return pd.Series(False, index=frame.index, dtype=bool)
    tier = frame.get("v2_evidence_tier", pd.Series("", index=frame.index)).fillna("").astype(str).str.upper()
    klass = frame.get("v2_candidate_class", pd.Series("", index=frame.index)).fillna("").astype(str).str.upper()
    return tier.str.contains("POSITIVE_SUPPORT", regex=False) | klass.str.contains("POSITIVE_SUPPORT", regex=False)


def background_informative_mask(frame: pd.DataFrame) -> pd.Series:
    if frame.empty:
        return pd.Series(False, index=frame.index, dtype=bool)
    gate = frame.get("v2_evidence_gate", pd.Series("", index=frame.index)).fillna("").astype(str).str.upper()
    return gate.isin({"BACKGROUND_INFORMATIVE", "SHADOW_BACKGROUND_CONTEXT"})


def _candidate_matches_truth(candidates: pd.DataFrame, truth_row) -> pd.DataFrame:
    if candidates.empty:
        return candidates
    subset = candidates[
        candidates["sample_id"].astype(str).eq(str(truth_row.sample_id))
        & candidates["chrom"].astype(str).eq(str(truth_row.chrom))
        & candidates["state"].astype(str).eq(str(truth_row.state))
    ].copy()
    if subset.empty or pd.isna(truth_row.start) or pd.isna(truth_row.end):
        return subset
    overlaps = subset.apply(
        lambda row: overlap_fraction(row["start"], row["end"], truth_row.start, truth_row.end),
        axis=1,
    )
    return subset[overlaps.fillna(0.0) > 0.0].copy()


def _count_contains(frame: pd.DataFrame, column: str, token: str) -> int:
    if frame.empty or column not in frame.columns:
        return 0
    return int(frame[column].fillna("").astype(str).str.upper().str.contains(token, regex=False).sum())


def _count_equals(frame: pd.DataFrame, column: str, value: str) -> int:
    if frame.empty or column not in frame.columns:
        return 0
    return int(frame[column].fillna("").astype(str).str.lower().eq(value.lower()).sum())


def _value_counts(frame: pd.DataFrame, column: str) -> dict[str, int]:
    if frame.empty or column not in frame.columns:
        return {}
    counts = frame[column].fillna("UNKNOWN").astype(str).value_counts().sort_index().to_dict()
    return {str(key): int(value) for key, value in counts.items()}


def _rule_tag_counts(frame: pd.DataFrame, column: str) -> dict[str, int]:
    if frame.empty or column not in frame.columns:
        return {}
    counts: dict[str, int] = {}
    for tag_text in frame[column].fillna("").astype(str):
        for tag in [part.strip() for part in tag_text.split(";") if part.strip()]:
            counts[tag] = counts.get(tag, 0) + 1
    return dict(sorted(counts.items()))


def build_v2_benchmark(
    classification_paths: list[str],
    truth_tsv: str,
    sample_ids: list[str],
    reference_id: str,
    output_truth_metrics: str,
    output_sample_summary: str,
    output_summary: str,
    output_filtered_events: str = "",
    output_filtered_events_json: str = "",
    output_report_events: str = "",
    output_report_events_json: str = "",
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    inferred_samples = [_sample_from_classification_path(path) for path in classification_paths]
    ordered_samples = list(dict.fromkeys([str(sample) for sample in sample_ids if str(sample)] + inferred_samples))
    classifications = load_classifications(classification_paths)
    truth = load_truth(truth_tsv, ordered_samples)

    truth_rows = []
    for truth_index, row in truth.iterrows():
        matches = _candidate_matches_truth(classifications, row)
        hard_suppressed = hard_suppression_mask(matches)
        report_layer_filtered = report_layer_filtered_mask(matches)
        preserved = matches.loc[~hard_suppressed & ~report_layer_filtered].copy()
        top = preserved.sort_values("a_abs_zscore", ascending=False).iloc[0] if not preserved.empty else None
        truth_rows.append(
            {
                "sample_id": str(row.sample_id),
                "truth_index": int(truth_index),
                "chrom": str(row.chrom),
                "start": row.start,
                "end": row.end,
                "state": str(row.state),
                "v2_overlap_count": int(len(matches)),
                "v2_preserved_count": int(len(preserved)),
                "v2_positive_support_count": int(positive_support_mask(matches).sum()) if not matches.empty else 0,
                "v2_background_informative_count": int(background_informative_mask(matches).sum()) if not matches.empty else 0,
                "v2_hard_suppressed_count": int(hard_suppressed.sum()) if not matches.empty else 0,
                "top_candidate_id": "" if top is None else str(top.get("candidate_id", "")),
                "top_a_abs_zscore": math.nan if top is None else top.get("a_abs_zscore", math.nan),
                "top_v2_candidate_class": "" if top is None else str(top.get("v2_candidate_class", "")),
                "top_v2_evidence_tier": "" if top is None else str(top.get("v2_evidence_tier", "")),
                "top_v2_classifier_action": "" if top is None else str(top.get("v2_classifier_action", "")),
                "top_v2_disposition": "" if top is None else str(top.get("v2_disposition", "")),
                "top_v2_filter_action": "" if top is None else str(top.get("v2_filter_action", "")),
                "top_v2_length_tier": "" if top is None else str(top.get("v2_length_tier", "")),
                "top_v2_clean_support_label": "" if top is None else str(top.get("v2_clean_support_label", "")),
                "top_v2_gc_rc_context_label": "" if top is None else str(top.get("v2_gc_rc_context_label", "")),
                "top_v2_b_signal_context_label": "" if top is None else str(top.get("v2_b_signal_context_label", "")),
                "top_v2_burden_reduction_tier": "" if top is None else str(top.get("v2_burden_reduction_tier", "")),
                "top_v2_burden_reduction_action": "" if top is None else str(top.get("v2_burden_reduction_action", "")),
                "top_v2_burden_reduction_reason": "" if top is None else str(top.get("v2_burden_reduction_reason", "")),
                "top_v2_burden_evidence_tags": "" if top is None else str(top.get("v2_burden_evidence_tags", "")),
                "top_v2_report_layer_class": "" if top is None else str(top.get("v2_report_layer_class", "")),
                "top_v2_report_visibility": "" if top is None else str(top.get("v2_report_visibility", "")),
                "top_v2_report_filter_reason": "" if top is None else str(top.get("v2_report_filter_reason", "")),
                "top_v2_report_filter_rule_tags": "" if top is None else str(top.get("v2_report_filter_rule_tags", "")),
                "top_attenuation_ratio": math.nan if top is None else top.get("attenuation_ratio", math.nan),
            }
        )
    truth_metrics = pd.DataFrame(truth_rows, columns=TRUTH_METRIC_COLUMNS)

    sample_rows = []
    for sample_id in ordered_samples:
        sample_class = classifications[classifications["sample_id"].astype(str).eq(str(sample_id))].copy()
        sample_truth = truth_metrics[truth_metrics["sample_id"].astype(str).eq(str(sample_id))].copy()
        truth_event_count = int(len(sample_truth))
        truth_preserved = int(sample_truth["v2_preserved_count"].gt(0).sum()) if not sample_truth.empty else 0
        report_event_count = _count_equals(sample_class, "v2_report_layer_class", "report_event")
        sample_report_burden_flag = int(report_event_count >= SAMPLE_REPORT_BURDEN_THRESHOLD)
        sample_report_burden_reason = (
            f"report_event_count_ge_{SAMPLE_REPORT_BURDEN_THRESHOLD}" if sample_report_burden_flag else ""
        )
        sample_rows.append(
            {
                "sample_id": str(sample_id),
                "truth_event_count": truth_event_count,
                "truth_preserved_count": truth_preserved,
                "FN_count": max(truth_event_count - truth_preserved, 0),
                "candidate_count": int(len(sample_class)),
                "v2_positive_support_count": int(positive_support_mask(sample_class).sum()) if not sample_class.empty else 0,
                "v2_background_compatible_count": _count_contains(sample_class, "v2_candidate_class", "BACKGROUND_COMPATIBLE"),
                "v2_technical_review_count": _count_contains(sample_class, "v2_candidate_class", "TECHNICAL_REVIEW"),
                "v2_contract_risk_count": _count_contains(sample_class, "v2_candidate_class", "CONTRACT_RISK"),
                "v2_hard_suppressed_count": int(hard_suppression_mask(sample_class).sum()) if not sample_class.empty else 0,
                "v2_filter_suppressed_count": _count_contains(sample_class, "v2_filter_action", "suppress_"),
                "v2_report_candidate_burden_count": _count_equals(
                    sample_class, "v2_burden_reduction_tier", "report_candidate"
                ),
                "v2_review_candidate_burden_count": (
                    _count_equals(sample_class, "v2_burden_reduction_tier", "review_candidate")
                    + _count_equals(sample_class, "v2_burden_reduction_tier", "background_unknown_review")
                ),
                "v2_background_unknown_review_burden_count": _count_equals(
                    sample_class, "v2_burden_reduction_tier", "background_unknown_review"
                ),
                "v2_technical_risk_burden_count": _count_equals(
                    sample_class, "v2_burden_reduction_tier", "technical_risk_review"
                ),
                "v2_branch_s_review_burden_count": _count_equals(
                    sample_class, "v2_burden_reduction_tier", "branch_s_review"
                ),
                "v2_report_event_count": report_event_count,
                "sample_report_burden_flag": sample_report_burden_flag,
                "sample_report_burden_reason": sample_report_burden_reason,
                "v2_internal_review_event_count": _count_equals(
                    sample_class, "v2_report_layer_class", "internal_review_event"
                ),
                "v2_filtered_event_count": _count_equals(sample_class, "v2_report_layer_class", "filtered_event"),
                "v2_branch_s_event_count": _count_equals(sample_class, "v2_report_layer_class", "branch_s_event"),
            }
        )
    sample_summary = pd.DataFrame(sample_rows, columns=SAMPLE_SUMMARY_COLUMNS)
    truth_event_count = int(len(truth_metrics))
    truth_preserved_count = int(truth_metrics["v2_preserved_count"].gt(0).sum()) if not truth_metrics.empty else 0
    summary = {
        "status": "ready" if truth_tsv else "skipped_no_truth",
        "reference_id": str(reference_id),
        "sample_count": int(len(sample_summary)),
        "candidate_count": int(sample_summary["candidate_count"].sum()) if not sample_summary.empty else 0,
        "truth_event_count": truth_event_count,
        "truth_preserved_count": truth_preserved_count,
        "FN_count": max(truth_event_count - truth_preserved_count, 0),
        "truth_recall_by_v2_preservation": (truth_preserved_count / truth_event_count if truth_event_count else None),
        "truth_hard_suppressed_count": int(truth_metrics["v2_hard_suppressed_count"].sum()) if not truth_metrics.empty else 0,
        "truth_report_layer_filtered_count": 0,
        "v2_positive_support_candidate_count": int(sample_summary["v2_positive_support_count"].sum()) if not sample_summary.empty else 0,
        "v2_filter_suppressed_candidate_count": int(sample_summary["v2_filter_suppressed_count"].sum()) if not sample_summary.empty else 0,
        "v2_report_candidate_burden_count": int(sample_summary["v2_report_candidate_burden_count"].sum()) if not sample_summary.empty else 0,
        "v2_review_candidate_burden_count": int(sample_summary["v2_review_candidate_burden_count"].sum()) if not sample_summary.empty else 0,
        "v2_background_unknown_review_burden_count": int(sample_summary["v2_background_unknown_review_burden_count"].sum()) if not sample_summary.empty else 0,
        "v2_technical_risk_burden_count": int(sample_summary["v2_technical_risk_burden_count"].sum()) if not sample_summary.empty else 0,
        "v2_branch_s_review_burden_count": int(sample_summary["v2_branch_s_review_burden_count"].sum()) if not sample_summary.empty else 0,
        "v2_report_event_count": int(sample_summary["v2_report_event_count"].sum()) if not sample_summary.empty else 0,
        "sample_report_burden_threshold": SAMPLE_REPORT_BURDEN_THRESHOLD,
        "sample_report_burden_flag_count": (
            int(sample_summary["sample_report_burden_flag"].sum()) if not sample_summary.empty else 0
        ),
        "v2_internal_review_event_count": int(sample_summary["v2_internal_review_event_count"].sum()) if not sample_summary.empty else 0,
        "v2_filtered_event_count": int(sample_summary["v2_filtered_event_count"].sum()) if not sample_summary.empty else 0,
        "v2_branch_s_event_count": int(sample_summary["v2_branch_s_event_count"].sum()) if not sample_summary.empty else 0,
        "burden_stratification_fields": [
            "v2_filter_action",
            "v2_disposition",
            "v2_length_tier",
            "v2_clean_support_label",
            "v2_b_signal_context_label",
            "v2_gc_rc_context_label",
            "v2_burden_reduction_tier",
            "v2_burden_reduction_action",
            "v2_report_layer_class",
            "v2_report_visibility",
        ],
        "burden_stratification_counts": {
            "v2_filter_action": _value_counts(classifications, "v2_filter_action"),
            "v2_disposition": _value_counts(classifications, "v2_disposition"),
            "v2_length_tier": _value_counts(classifications, "v2_length_tier"),
            "v2_clean_support_label": _value_counts(classifications, "v2_clean_support_label"),
            "v2_b_signal_context_label": _value_counts(classifications, "v2_b_signal_context_label"),
            "v2_gc_rc_context_label": _value_counts(classifications, "v2_gc_rc_context_label"),
            "v2_burden_reduction_tier": _value_counts(classifications, "v2_burden_reduction_tier"),
            "v2_burden_reduction_action": _value_counts(classifications, "v2_burden_reduction_action"),
            "v2_report_layer_class": _value_counts(classifications, "v2_report_layer_class"),
            "v2_report_visibility": _value_counts(classifications, "v2_report_visibility"),
        },
        "legacy_branch_b_decision_fields_used": False,
        "ignored_legacy_decision_fields": [
            "final_disposition",
            "branch_b_keep_event",
            "branch_b_report_class",
            "branch_b_artifact_status",
        ],
        "final_report_impact": "none_shadow_only",
        "benchmark_scope": "v2_classifier_rows_only",
    }
    if not truth_metrics.empty:
        filtered_truth = 0
        for _, row in truth.iterrows():
            matches = _candidate_matches_truth(classifications, row)
            filtered_truth += int(report_layer_filtered_mask(matches).sum()) if not matches.empty else 0
        summary["truth_report_layer_filtered_count"] = int(filtered_truth)
    else:
        summary["truth_report_layer_filtered_count"] = 0

    filtered_classifications = classifications.loc[report_layer_filtered_mask(classifications)].copy()
    filtered_payload = {
        "filtered_event_count": int(len(filtered_classifications)),
        "filtered_event_rule_counts": _rule_tag_counts(filtered_classifications, "v2_report_filter_rule_tags"),
        "filtered_event_ids": (
            filtered_classifications.get("candidate_id", pd.Series([], dtype=object)).fillna("").astype(str).tolist()
            if not filtered_classifications.empty
            else []
        ),
    }
    report_classifications = classifications.loc[
        classifications.get("v2_report_layer_class", pd.Series("", index=classifications.index))
        .fillna("")
        .astype(str)
        .eq("report_event")
    ].copy()
    report_payload = {
        "report_event_count": int(len(report_classifications)),
        "report_event_rule_counts": _rule_tag_counts(report_classifications, "v2_report_filter_rule_tags"),
        "report_event_ids": (
            report_classifications.get("candidate_id", pd.Series([], dtype=object)).fillna("").astype(str).tolist()
            if not report_classifications.empty
            else []
        ),
    }

    for path_value, frame in [
        (output_truth_metrics, truth_metrics),
        (output_sample_summary, sample_summary),
        (output_filtered_events, filtered_classifications),
        (output_report_events, report_classifications),
    ]:
        if not path_value:
            continue
        path = Path(path_value)
        path.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(path, sep="\t", index=False)
    if output_filtered_events_json:
        path = Path(output_filtered_events_json)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(filtered_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if output_report_events_json:
        path = Path(output_report_events_json)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(report_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    summary_path = Path(output_summary)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return truth_metrics, sample_summary, summary


def main():
    args = parse_args()
    build_v2_benchmark(
        classification_paths=args.classification_tsv,
        truth_tsv=args.truth_tsv,
        sample_ids=args.sample_id,
        reference_id=args.reference_id,
        output_truth_metrics=args.output_truth_metrics,
        output_sample_summary=args.output_sample_summary,
        output_summary=args.output_summary,
        output_filtered_events=args.output_filtered_events,
        output_filtered_events_json=args.output_filtered_events_json,
        output_report_events=args.output_report_events,
        output_report_events_json=args.output_report_events_json,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
