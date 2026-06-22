#!/biosoftware/miniconda/envs/snakemake_env/bin/python
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import pandas as pd


AUDIT_COLUMNS = [
    "sample_id",
    "candidate_id",
    "event_id",
    "chrom",
    "start",
    "end",
    "state",
    "original_report_layer_class",
    "original_report_visibility",
    "plot_layer_class",
    "plot_visibility",
    "plot_visibility_reason",
    "plot_support_class",
    "support_interpretation_status",
    "cn_direction_consistency_status",
    "same_direction_median_display_ref_z",
    "cn_same_direction_fraction",
    "truth_guard_status",
    "protected_event_label",
    "proposed_report_table_action",
    "proposed_report_table_reason",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Audit plot-manifest low-confidence events before changing report-table classification."
    )
    parser.add_argument("--report-events-tsv", required=True)
    parser.add_argument("--plot-manifest-tsv", action="append", default=[])
    parser.add_argument("--plot-support-tsv", action="append", default=[])
    parser.add_argument("--truth-metrics-tsv", default="")
    parser.add_argument("--sample-summary-tsv", default="")
    parser.add_argument("--reference-id", default="UNKNOWN_REFERENCE")
    parser.add_argument("--output-audit-tsv", required=True)
    parser.add_argument("--output-summary-json", required=True)
    parser.add_argument("--output-report-md", required=True)
    return parser.parse_args()


def _read_table(path_value: str | Path, columns: list[str] | None = None) -> pd.DataFrame:
    if not path_value:
        return pd.DataFrame(columns=columns or [])
    path = Path(path_value)
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame(columns=columns or [])
    try:
        return pd.read_csv(path, sep="\t")
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=columns or [])


def _concat_tables(frames_or_paths, columns: list[str] | None = None) -> pd.DataFrame:
    frames = []
    for item in frames_or_paths or []:
        if isinstance(item, pd.DataFrame):
            frame = item.copy()
        else:
            frame = _read_table(item, columns=columns)
        if not frame.empty:
            frames.append(frame)
    if not frames:
        return pd.DataFrame(columns=columns or [])
    return pd.concat(frames, ignore_index=True, sort=False)


def _normalize_chrom(value) -> str:
    token = str(value or "").strip()
    if not token:
        return ""
    if token.lower().startswith("chr"):
        token = token[3:]
    token = token.upper()
    if token == "23":
        token = "X"
    elif token == "24":
        token = "Y"
    return f"chr{token}"


def _prepare_events(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    if out.empty:
        return pd.DataFrame(columns=["sample_id", "candidate_id", "chrom", "start", "end", "state"])
    if "sample" in out.columns and "sample_id" not in out.columns:
        out = out.rename(columns={"sample": "sample_id"})
    for column in ["sample_id", "candidate_id", "event_id", "chrom", "state"]:
        if column not in out.columns:
            out[column] = ""
    for column in ["start", "end"]:
        if column not in out.columns:
            out[column] = math.nan
        out[column] = pd.to_numeric(out[column], errors="coerce")
    out["sample_id"] = out["sample_id"].fillna("").astype(str)
    out["candidate_id"] = out["candidate_id"].fillna("").astype(str)
    out["event_id"] = out["event_id"].fillna("").astype(str)
    out["chrom"] = out["chrom"].map(_normalize_chrom)
    out["state"] = out["state"].fillna("").astype(str).str.lower().replace({"dup": "gain", "del": "loss"})
    return out


def _merge_by_candidate_id(left: pd.DataFrame, right: pd.DataFrame, suffix: str) -> pd.DataFrame:
    if left.empty or right.empty or "candidate_id" not in right.columns:
        return left
    usable = right.copy()
    usable["candidate_id"] = usable["candidate_id"].fillna("").astype(str)
    usable = usable[usable["candidate_id"].ne("")]
    if usable.empty:
        return left
    usable = usable.drop_duplicates("candidate_id", keep="first")
    return left.merge(usable, on="candidate_id", how="left", suffixes=("", suffix))


def _truth_candidate_ids(truth_metrics: pd.DataFrame) -> set[str]:
    if truth_metrics is None or truth_metrics.empty:
        return set()
    ids: set[str] = set()
    for column in ["top_candidate_id", "candidate_id"]:
        if column not in truth_metrics.columns:
            continue
        values = truth_metrics[column].fillna("").astype(str)
        ids.update(value for value in values if value)
    return ids


def _summary_truth_counts(truth_metrics: pd.DataFrame, sample_summary: pd.DataFrame) -> dict[str, int]:
    if sample_summary is not None and not sample_summary.empty:
        truth_event_count = int(pd.to_numeric(sample_summary.get("truth_event_count", 0), errors="coerce").fillna(0).sum())
        truth_preserved_count = int(
            pd.to_numeric(sample_summary.get("truth_preserved_count", 0), errors="coerce").fillna(0).sum()
        )
        fn_count = int(pd.to_numeric(sample_summary.get("FN_count", 0), errors="coerce").fillna(0).sum())
    elif truth_metrics is not None and not truth_metrics.empty:
        truth_event_count = int(len(truth_metrics))
        preserved = pd.to_numeric(truth_metrics.get("v2_preserved_count", 0), errors="coerce").fillna(0)
        truth_preserved_count = int(preserved.gt(0).sum())
        fn_count = max(truth_event_count - truth_preserved_count, 0)
    else:
        truth_event_count = 0
        truth_preserved_count = 0
        fn_count = 0
    hard_suppressed = 0
    if truth_metrics is not None and not truth_metrics.empty and "v2_hard_suppressed_count" in truth_metrics.columns:
        hard_suppressed = int(pd.to_numeric(truth_metrics["v2_hard_suppressed_count"], errors="coerce").fillna(0).sum())
    return {
        "truth_event_count": truth_event_count,
        "truth_preserved_count": truth_preserved_count,
        "FN_count": fn_count,
        "truth_hard_suppressed_count": hard_suppressed,
    }


def _classify_action(row: pd.Series, has_truth: bool) -> tuple[str, str, str, str]:
    candidate_id = str(row.get("candidate_id", ""))
    truth_guard = bool(row.get("_truth_guard", False))
    if not has_truth:
        truth_status = "context_only_no_truth"
        protected_label = ""
    elif truth_guard:
        truth_status = "truth_overlap_protected"
        protected_label = "locked_truth_top_candidate"
    else:
        truth_status = "no_truth_overlap"
        protected_label = ""

    support_class = str(row.get("plot_support_class", ""))
    if truth_guard:
        return (
            truth_status,
            protected_label,
            "retain_report_event_truth_guard",
            f"{candidate_id} overlaps locked truth and must not be downgraded by plot-support ablation",
        )
    if support_class == "Z_SUPPORTED_CN_NOT_SUPPORTED":
        return (
            truth_status,
            protected_label,
            "downgrade_to_internal_review_candidate",
            "autosomal report_event has z support but CN direction is not supported; downgrade candidate for audit only",
        )
    return (
        truth_status,
        protected_label,
        "retain_report_event",
        "no low-confidence report-table ablation trigger",
    )


def build_report_table_ablation_audit(
    report_events: pd.DataFrame,
    plot_manifests,
    plot_supports,
    truth_metrics: pd.DataFrame | None = None,
    sample_summary: pd.DataFrame | None = None,
    reference_id: str = "UNKNOWN_REFERENCE",
    output_audit_tsv: str | Path | None = None,
    output_summary_json: str | Path | None = None,
    output_report_md: str | Path | None = None,
):
    reports = _prepare_events(report_events)
    manifests = _prepare_events(_concat_tables(plot_manifests))
    supports = _prepare_events(_concat_tables(plot_supports))
    truth_metrics = truth_metrics.copy() if truth_metrics is not None else pd.DataFrame()
    sample_summary = sample_summary.copy() if sample_summary is not None else pd.DataFrame()

    audit = reports.copy()
    if "v2_report_layer_class" in audit.columns:
        audit = audit[audit["v2_report_layer_class"].fillna("").astype(str).eq("report_event")].copy()
    audit["original_report_layer_class"] = audit.get("v2_report_layer_class", "report_event")
    audit["original_report_visibility"] = audit.get("v2_report_visibility", "")

    manifest_columns = [
        "candidate_id",
        "event_id",
        "plot_layer_class",
        "plot_visibility",
        "plot_visibility_reason",
        "plot_support_class",
        "event_layer",
    ]
    if not manifests.empty:
        audit = _merge_by_candidate_id(audit, manifests[[c for c in manifest_columns if c in manifests.columns]], "_manifest")

    support_columns = [
        "candidate_id",
        "event_id",
        "support_interpretation_status",
        "cn_direction_consistency_status",
        "same_direction_median_display_ref_z",
        "cn_same_direction_fraction",
    ]
    if not supports.empty:
        audit = _merge_by_candidate_id(audit, supports[[c for c in support_columns if c in supports.columns]], "_support")

    for column in [
        "event_id",
        "plot_layer_class",
        "plot_visibility",
        "plot_visibility_reason",
        "plot_support_class",
        "support_interpretation_status",
        "cn_direction_consistency_status",
    ]:
        if column not in audit.columns:
            audit[column] = ""
        audit[column] = audit[column].fillna("").astype(str)
    for column in ["same_direction_median_display_ref_z", "cn_same_direction_fraction"]:
        if column not in audit.columns:
            audit[column] = math.nan
        audit[column] = pd.to_numeric(audit[column], errors="coerce")

    truth_ids = _truth_candidate_ids(truth_metrics)
    has_truth = bool(len(truth_metrics) > 0)
    audit["_truth_guard"] = audit["candidate_id"].fillna("").astype(str).isin(truth_ids)
    action_rows = audit.apply(lambda row: _classify_action(row, has_truth=has_truth), axis=1)
    audit["truth_guard_status"] = [item[0] for item in action_rows]
    audit["protected_event_label"] = [item[1] for item in action_rows]
    audit["proposed_report_table_action"] = [item[2] for item in action_rows]
    audit["proposed_report_table_reason"] = [item[3] for item in action_rows]
    audit = audit.drop(columns=["_truth_guard"], errors="ignore")

    for column in AUDIT_COLUMNS:
        if column not in audit.columns:
            audit[column] = ""
    audit = audit[AUDIT_COLUMNS].copy()

    truth_counts = _summary_truth_counts(truth_metrics, sample_summary)
    original_count = int(len(audit))
    demotion_count = int(audit["proposed_report_table_action"].eq("downgrade_to_internal_review_candidate").sum())
    final_plot_count = int(audit["plot_visibility"].eq("final_report_plot").sum())
    summary = {
        "status": "ready" if has_truth else "context_only_no_truth",
        "reference_id": str(reference_id),
        "original_report_event_count": original_count,
        "final_report_plot_event_count": final_plot_count,
        "proposed_report_event_count": max(original_count - demotion_count, 0),
        "proposed_internal_review_demotion_count": demotion_count,
        "truth_guarded_report_event_count": int(audit["truth_guard_status"].eq("truth_overlap_protected").sum()),
        "plot_support_class_counts": audit["plot_support_class"].fillna("").astype(str).value_counts().to_dict(),
        "proposed_action_counts": audit["proposed_report_table_action"].fillna("").astype(str).value_counts().to_dict(),
        "TP_FN_FP_status": "not_computed_no_truth" if not has_truth else "truth_preservation_only_no_fp_claim",
        "legacy_branch_b_decision_fields_used": False,
        "final_report_impact": "ablation_audit_only_not_production_filtering",
    }
    summary.update(truth_counts)

    if output_audit_tsv:
        path = Path(output_audit_tsv)
        path.parent.mkdir(parents=True, exist_ok=True)
        audit.to_csv(path, sep="\t", index=False)
    if output_summary_json:
        path = Path(output_summary_json)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if output_report_md:
        _write_markdown_report(output_report_md, summary)
    return audit, summary


def _write_markdown_report(path_value: str | Path, summary: dict):
    path = Path(path_value)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Branch B V2 Report-Table Ablation Audit",
        "",
        "This is an ablation audit, not production filtering.",
        "",
        f"- reference_id: `{summary.get('reference_id', '')}`",
        f"- original report events: {summary.get('original_report_event_count', 0)}",
        f"- final plot events: {summary.get('final_report_plot_event_count', 0)}",
        f"- proposed internal-review demotions: {summary.get('proposed_internal_review_demotion_count', 0)}",
        f"- truth events: {summary.get('truth_event_count', 0)}",
        f"- FN_count: {summary.get('FN_count', 0)}",
        f"- truth hard-suppressed: {summary.get('truth_hard_suppressed_count', 0)}",
        "",
        "0615/context cohorts remain burden-only; TP/FN/FP is not computed when truth is absent.",
        "",
        "Candidate rule under audit: autosomal `report_event` with "
        "`plot_support_class=Z_SUPPORTED_CN_NOT_SUPPORTED` may be downgraded to internal review, "
        "while locked-truth candidates remain protected.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    args = parse_args()
    build_report_table_ablation_audit(
        report_events=_read_table(args.report_events_tsv),
        plot_manifests=args.plot_manifest_tsv,
        plot_supports=args.plot_support_tsv,
        truth_metrics=_read_table(args.truth_metrics_tsv),
        sample_summary=_read_table(args.sample_summary_tsv),
        reference_id=args.reference_id,
        output_audit_tsv=args.output_audit_tsv,
        output_summary_json=args.output_summary_json,
        output_report_md=args.output_report_md,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
