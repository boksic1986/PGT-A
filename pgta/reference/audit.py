#!/biosoftware/miniconda/envs/snakemake_env/bin/python
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

AUDIT_COLUMNS = [
    "sample",
    "batch",
    "qc_status",
    "sex_call",
    "reference_membership",
    "branch_a_candidate_count",
    "branch_a_strong_count",
    "branch_a_sensitive_count",
    "broad_event_count",
    "altered_genome_fraction",
    "acrocentric_or_high_repeat_burden",
    "centromere_telomere_proximal_burden",
    "shared_batch_signal_flag",
    "sample_specific_signal_flag",
    "R_label",
    "R_label_reason",
    "recommended_action",
]

ACROCENTRIC_CHROMS = {"chr13", "chr14", "chr15", "chr21", "chr22"}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Aggregate Branch A/reference evidence for G1 reference cohort audit."
    )
    parser.add_argument("--sample-id", action="append", default=[])
    parser.add_argument("--sample-metadata", default="")
    parser.add_argument("--a-candidates-dir", required=True)
    parser.add_argument("--qc-dir", required=True)
    parser.add_argument("--gender-dir", default="")
    parser.add_argument("--evidence-ledger-dir", default="")
    parser.add_argument("--reference-samples-file", default="")
    parser.add_argument("--bin-annotations", default="")
    parser.add_argument("--reference-id", default="UNKNOWN_REFERENCE")
    parser.add_argument("--output-audit", required=True)
    parser.add_argument("--output-summary", required=True)
    parser.add_argument("--strong-z", type=float, default=10.0)
    parser.add_argument("--broad-event-min-bp", type=int, default=10_000_000)
    parser.add_argument("--sample-specific-fraction-threshold", type=float, default=0.05)
    parser.add_argument("--high-risk-burden-threshold", type=float, default=0.20)
    parser.add_argument("--shared-signal-min-samples", type=int, default=3)
    return parser.parse_args()


def _read_table(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t")
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def _safe_float(value, default=0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return float(default)
    if not np.isfinite(out):
        return float(default)
    return out


def _safe_int(value, default=0) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return int(default)


def normalize_chrom(value) -> str:
    text = str(value).strip()
    if not text:
        return text
    return text if text.startswith("chr") else f"chr{text}"


def load_sample_metadata(path_value: str, sample_ids: list[str]) -> pd.DataFrame:
    rows = []
    if path_value:
        path = Path(path_value)
        df = _read_table(path)
        if df.empty or "sample_id" not in df.columns:
            raise ValueError(f"sample metadata must contain sample_id: {path}")
        for _, row in df.iterrows():
            sample_id = str(row.get("sample_id", "")).strip()
            if not sample_id:
                continue
            rows.append(
                {
                    "sample_id": sample_id,
                    "batch": str(row.get("batch", row.get("batch_group", "unknown")) or "unknown"),
                    "reference_membership": str(row.get("reference_membership", "candidate_only") or "candidate_only"),
                    "known_status": str(row.get("known_status", "unknown") or "unknown").lower(),
                    "manual_r_label": str(row.get("manual_r_label", "") or "").upper(),
                    "manual_r_label_reason": str(row.get("manual_r_label_reason", "") or ""),
                }
            )
    else:
        for sample_id in sample_ids:
            rows.append(
                {
                    "sample_id": str(sample_id),
                    "batch": "unknown",
                    "reference_membership": "candidate_only",
                    "known_status": "unknown",
                    "manual_r_label": "",
                    "manual_r_label_reason": "",
                }
            )
    if sample_ids:
        allowed = {str(item) for item in sample_ids}
        rows = [row for row in rows if row["sample_id"] in allowed]
    if not rows:
        raise ValueError("reference audit has no samples to evaluate")
    return pd.DataFrame(rows).drop_duplicates(subset=["sample_id"], keep="first")


def load_reference_membership(path_value: str) -> set[str]:
    if not path_value:
        return set()
    path = Path(path_value)
    if not path.exists():
        raise FileNotFoundError(f"reference sample list not found: {path}")
    return {line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()}


def infer_genome_bp(path_value: str, fallback=3_000_000_000) -> float:
    if not path_value:
        return float(fallback)
    df = _read_table(Path(path_value))
    if df.empty:
        return float(fallback)
    chrom = df.get("chrom", pd.Series(dtype=str)).map(normalize_chrom)
    mask = ~chrom.isin({"chrM", "chrMT"})
    if "effective_size" in df.columns:
        total = pd.to_numeric(df.loc[mask, "effective_size"], errors="coerce").fillna(0).sum()
    else:
        start = pd.to_numeric(df.loc[mask, "start"], errors="coerce")
        end = pd.to_numeric(df.loc[mask, "end"], errors="coerce")
        total = (end - start).clip(lower=0).fillna(0).sum()
    return float(total) if total > 0 else float(fallback)


def sample_file(directory: str, sample_id: str, suffix: str) -> Path:
    return Path(directory) / f"{sample_id}{suffix}"


def read_qc_status(qc_dir: str, sample_id: str) -> tuple[str, str]:
    df = _read_table(sample_file(qc_dir, sample_id, ".qc.tsv"))
    if df.empty:
        return "UNKNOWN", "missing_qc"
    row = df.iloc[0]
    return str(row.get("status", "UNKNOWN") or "UNKNOWN").upper(), str(row.get("reason", "") or "")


def read_sex_call(gender_dir: str, sample_id: str) -> str:
    if not gender_dir:
        return "UNKNOWN"
    df = _read_table(sample_file(gender_dir, sample_id, ".gender.tsv"))
    if df.empty:
        return "UNKNOWN"
    return str(df.iloc[0].get("sex_call", "UNKNOWN") or "UNKNOWN").upper()


def load_candidates(a_candidates_dir: str, sample_id: str) -> pd.DataFrame:
    df = _read_table(sample_file(a_candidates_dir, sample_id, ".candidate_events.tsv"))
    if df.empty:
        return pd.DataFrame(columns=["chrom", "start", "end", "state", "a_abs_zscore", "a_support_level"])
    for column in ["start", "end", "a_abs_zscore"]:
        if column in df.columns:
            df[column] = pd.to_numeric(df[column], errors="coerce")
    if "chrom" in df.columns:
        df["chrom"] = df["chrom"].map(normalize_chrom)
    if "event_length_bp" not in df.columns and {"start", "end"}.issubset(df.columns):
        df["event_length_bp"] = (df["end"] - df["start"]).clip(lower=0)
    return df


def load_evidence_ledger(evidence_dir: str, sample_id: str) -> pd.DataFrame:
    if not evidence_dir:
        return pd.DataFrame()
    df = _read_table(sample_file(evidence_dir, sample_id, ".candidate_evidence.tsv"))
    if df.empty:
        return df
    for column in df.columns:
        if column.endswith("_fraction") or column.endswith("_burden") or column in {
            "hard_region_fraction",
            "region_risk_score_mean",
            "region_risk_score_max",
            "cnvseq_gap_centromere_bin_fraction",
            "a_abs_zscore",
        }:
            df[column] = pd.to_numeric(df[column], errors="coerce")
    if "chrom" in df.columns:
        df["chrom"] = df["chrom"].map(normalize_chrom)
    if "event_length_bp" not in df.columns and {"start", "end"}.issubset(df.columns):
        start = pd.to_numeric(df["start"], errors="coerce")
        end = pd.to_numeric(df["end"], errors="coerce")
        df["event_length_bp"] = (end - start).clip(lower=0)
    return df


def burden_max(df: pd.DataFrame, columns: list[str]) -> float:
    values = []
    for column in columns:
        if column in df.columns:
            values.append(pd.to_numeric(df[column], errors="coerce"))
    if not values:
        return 0.0
    merged = pd.concat(values, axis=1).max(axis=1).max()
    return _safe_float(merged, default=0.0)


def build_signal_keys(sample_candidates: dict[str, pd.DataFrame], broad_event_min_bp: int) -> dict[str, set[str]]:
    key_to_samples: dict[str, set[str]] = {}
    for sample_id, candidates in sample_candidates.items():
        if candidates.empty:
            continue
        for _, row in candidates.iterrows():
            length_bp = _safe_float(row.get("event_length_bp", _safe_float(row.get("end")) - _safe_float(row.get("start"))))
            if length_bp < broad_event_min_bp:
                continue
            chrom = normalize_chrom(row.get("chrom", ""))
            state = str(row.get("state", "")).lower()
            if not chrom or not state:
                continue
            key_to_samples.setdefault(f"{chrom}:{state}", set()).add(sample_id)
    return key_to_samples


def summarize_sample(
    sample_row: pd.Series,
    candidates: pd.DataFrame,
    ledger: pd.DataFrame,
    qc_status: str,
    sex_call: str,
    reference_members: set[str],
    shared_keys: dict[str, set[str]],
    genome_bp: float,
    args,
) -> dict[str, object]:
    sample_id = str(sample_row["sample_id"])
    evidence = ledger if not ledger.empty else candidates
    event_lengths = pd.to_numeric(candidates.get("event_length_bp", pd.Series(dtype=float)), errors="coerce").fillna(0)
    broad_count = int((event_lengths >= int(args.broad_event_min_bp)).sum())
    altered_fraction = float(event_lengths.sum() / genome_bp) if genome_bp else 0.0
    strong_count = int(candidates.get("a_support_level", pd.Series(dtype=str)).astype(str).str.lower().eq("strong").sum())
    sensitive_count = int(candidates.get("a_support_level", pd.Series(dtype=str)).astype(str).str.lower().eq("sensitive").sum())
    acrocentric_candidates = 0
    acrocentric_bp = 0.0
    if not candidates.empty and "chrom" in candidates.columns:
        acrocentric_mask = candidates["chrom"].map(normalize_chrom).isin(ACROCENTRIC_CHROMS)
        acrocentric_candidates = int(acrocentric_mask.sum())
        acrocentric_bp = float(event_lengths.loc[acrocentric_mask].sum())
    high_repeat_burden = burden_max(
        evidence,
        [
            "hard_region_fraction",
            "segmental_duplication_overlap_fraction",
            "low_mappability_overlap_fraction",
            "repeat_rich_overlap_fraction",
            "blacklist_overlap_fraction",
            "ambiguous_alignment_overlap_fraction",
            "region_risk_score_mean",
            "region_risk_score_max",
        ],
    )
    acrocentric_genome_fraction = float(acrocentric_bp / genome_bp) if genome_bp else 0.0
    acrocentric_or_high_repeat_burden = max(high_repeat_burden, acrocentric_genome_fraction)
    centromere_telomere_burden = burden_max(
        evidence,
        [
            "gap_centromere_telomere_overlap_fraction",
            "cnvseq_gap_centromere_bin_fraction",
        ],
    )
    shared_flag = "no"
    sample_specific_flag = "no"
    for _, row in candidates.iterrows():
        length_bp = _safe_float(row.get("event_length_bp", _safe_float(row.get("end")) - _safe_float(row.get("start"))))
        if length_bp < int(args.broad_event_min_bp):
            continue
        key = f"{normalize_chrom(row.get('chrom', ''))}:{str(row.get('state', '')).lower()}"
        sample_count = len(shared_keys.get(key, set()))
        if sample_count >= int(args.shared_signal_min_samples):
            shared_flag = "yes"
        elif _safe_float(row.get("a_abs_zscore"), default=0.0) >= float(args.strong_z):
            sample_specific_flag = "yes"
    if altered_fraction >= float(args.sample_specific_fraction_threshold) and shared_flag != "yes":
        sample_specific_flag = "yes"

    reference_membership = str(sample_row.get("reference_membership", "candidate_only") or "candidate_only")
    if sample_id in reference_members:
        reference_membership = "included_in_evaluated_reference"
    known_status = str(sample_row.get("known_status", "unknown") or "unknown").lower()
    manual_label = str(sample_row.get("manual_r_label", "") or "").upper()
    manual_reason = str(sample_row.get("manual_r_label_reason", "") or "")

    r_label, reason, action = assign_r_label(
        qc_status=qc_status,
        known_status=known_status,
        manual_label=manual_label,
        manual_reason=manual_reason,
        candidate_count=int(len(candidates)),
        broad_count=broad_count,
        altered_fraction=altered_fraction,
        high_risk_burden=acrocentric_or_high_repeat_burden,
        centromere_telomere_burden=centromere_telomere_burden,
        shared_flag=shared_flag,
        sample_specific_flag=sample_specific_flag,
        args=args,
    )

    return {
        "sample": sample_id,
        "batch": str(sample_row.get("batch", "unknown") or "unknown"),
        "qc_status": qc_status,
        "sex_call": sex_call,
        "reference_membership": reference_membership,
        "branch_a_candidate_count": int(len(candidates)),
        "branch_a_strong_count": strong_count,
        "branch_a_sensitive_count": sensitive_count,
        "broad_event_count": broad_count,
        "altered_genome_fraction": round(altered_fraction, 6),
        "acrocentric_or_high_repeat_burden": round(float(acrocentric_or_high_repeat_burden), 6),
        "centromere_telomere_proximal_burden": round(float(centromere_telomere_burden), 6),
        "shared_batch_signal_flag": shared_flag,
        "sample_specific_signal_flag": sample_specific_flag,
        "R_label": r_label,
        "R_label_reason": reason,
        "recommended_action": action,
    }


def assign_r_label(
    qc_status: str,
    known_status: str,
    manual_label: str,
    manual_reason: str,
    candidate_count: int,
    broad_count: int,
    altered_fraction: float,
    high_risk_burden: float,
    centromere_telomere_burden: float,
    shared_flag: str,
    sample_specific_flag: str,
    args,
) -> tuple[str, str, str]:
    if manual_label in {"R0", "R1", "R2"}:
        return manual_label, f"manual_label:{manual_reason or manual_label}", "use_manual_review_label"
    if known_status in {"positive", "truth_positive", "abnormal", "exclude"}:
        return "R2", f"known_status={known_status}", "exclude_from_reference"
    if qc_status not in {"PASS", "WARN"}:
        return "R2", f"qc_status={qc_status}", "exclude_until_qc_resolved"
    reasons = []
    if sample_specific_flag == "yes":
        reasons.append("sample_specific_branch_a_signal")
    if broad_count > 0 and shared_flag != "yes":
        reasons.append("broad_signal_not_shared_across_batch")
    if altered_fraction >= float(args.sample_specific_fraction_threshold):
        reasons.append(f"altered_genome_fraction>={args.sample_specific_fraction_threshold:g}")
    if high_risk_burden >= float(args.high_risk_burden_threshold):
        reasons.append(f"high_repeat_or_acrocentric_burden>={args.high_risk_burden_threshold:g}")
    if centromere_telomere_burden >= float(args.high_risk_burden_threshold):
        reasons.append(f"centromere_telomere_burden>={args.high_risk_burden_threshold:g}")
    if reasons:
        return "R1", ";".join(dict.fromkeys(reasons)), "shadow_or_ablation_only"
    if candidate_count == 0:
        return "R0", "qc_pass_warn_no_branch_a_candidates", "candidate_reference"
    if shared_flag == "yes":
        return "R0", "branch_a_signal_shared_across_batch_not_sample_specific", "candidate_reference_with_batch_context"
    return "R0", "qc_pass_warn_no_sample_specific_signal", "candidate_reference"


def build_reference_audit(args) -> tuple[pd.DataFrame, dict[str, object]]:
    metadata = load_sample_metadata(args.sample_metadata, args.sample_id)
    sample_ids = metadata["sample_id"].astype(str).tolist()
    reference_members = load_reference_membership(args.reference_samples_file)
    genome_bp = infer_genome_bp(args.bin_annotations)
    candidates_by_sample = {
        sample_id: load_candidates(args.a_candidates_dir, sample_id)
        for sample_id in sample_ids
    }
    shared_keys = build_signal_keys(candidates_by_sample, int(args.broad_event_min_bp))
    rows = []
    for _, sample_row in metadata.iterrows():
        sample_id = str(sample_row["sample_id"])
        qc_status, _qc_reason = read_qc_status(args.qc_dir, sample_id)
        sex_call = read_sex_call(args.gender_dir, sample_id)
        ledger = load_evidence_ledger(args.evidence_ledger_dir, sample_id)
        rows.append(
            summarize_sample(
                sample_row=sample_row,
                candidates=candidates_by_sample[sample_id],
                ledger=ledger,
                qc_status=qc_status,
                sex_call=sex_call,
                reference_members=reference_members,
                shared_keys=shared_keys,
                genome_bp=genome_bp,
                args=args,
            )
        )
    audit = pd.DataFrame(rows, columns=AUDIT_COLUMNS)
    label_counts = audit["R_label"].value_counts(dropna=False).to_dict() if not audit.empty else {}
    summary = {
        "reference_id": args.reference_id,
        "sample_count": int(len(audit)),
        "label_counts": {str(key): int(value) for key, value in label_counts.items()},
        "strong_z": float(args.strong_z),
        "broad_event_min_bp": int(args.broad_event_min_bp),
        "sample_specific_fraction_threshold": float(args.sample_specific_fraction_threshold),
        "high_risk_burden_threshold": float(args.high_risk_burden_threshold),
        "shared_signal_min_samples": int(args.shared_signal_min_samples),
        "formal_n0_used": False,
        "legacy_branch_b_kept_counts_used_for_r_label": False,
        "notes": [
            "R labels are audit labels for reference-rebuild eligibility, not N0/N1/N2 background labels.",
            "Branch B kept counts are not used for R-label assignment.",
        ],
    }
    return audit, summary


def main():
    args = parse_args()
    audit, summary = build_reference_audit(args)
    output_audit = Path(args.output_audit)
    output_summary = Path(args.output_summary)
    output_audit.parent.mkdir(parents=True, exist_ok=True)
    output_summary.parent.mkdir(parents=True, exist_ok=True)
    audit.to_csv(output_audit, sep="\t", index=False)
    output_summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
