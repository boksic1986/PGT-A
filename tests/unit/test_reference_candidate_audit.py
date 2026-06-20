from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from pgta.reference.audit import build_reference_audit


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def base_args(tmp_path: Path, sample_ids: list[str]) -> argparse.Namespace:
    a_dir = tmp_path / "a_branch"
    qc_dir = tmp_path / "qc"
    gender_dir = tmp_path / "gender"
    for sample_id in sample_ids:
        write_tsv(
            qc_dir / f"{sample_id}.qc.tsv",
            [{"sample_id": sample_id, "status": "PASS", "reason": ""}],
        )
        write_tsv(
            gender_dir / f"{sample_id}.gender.tsv",
            [{"sample_id": sample_id, "sex_call": "XY"}],
        )
    return argparse.Namespace(
        sample_id=sample_ids,
        sample_metadata="",
        a_candidates_dir=str(a_dir),
        qc_dir=str(qc_dir),
        gender_dir=str(gender_dir),
        evidence_ledger_dir="",
        reference_samples_file="",
        bin_annotations="",
        reference_id="test_ref",
        output_audit=str(tmp_path / "audit.tsv"),
        output_summary=str(tmp_path / "audit.summary.json"),
        strong_z=10.0,
        broad_event_min_bp=10_000_000,
        sample_specific_fraction_threshold=0.05,
        high_risk_burden_threshold=0.20,
        shared_signal_min_samples=3,
    )


def candidate(
    sample_id: str,
    chrom: str,
    start: int,
    end: int,
    state: str = "gain",
    zscore: float = 12.0,
) -> dict[str, object]:
    return {
        "candidate_id": f"{sample_id}.A0001",
        "sample_id": sample_id,
        "chrom": chrom,
        "start": start,
        "end": end,
        "state": state,
        "a_zscore": zscore,
        "a_abs_zscore": abs(zscore),
        "a_support_level": "strong" if abs(zscore) >= 10 else "sensitive",
    }


def test_reference_audit_does_not_use_legacy_branch_b_kept_counts(tmp_path: Path):
    args = base_args(tmp_path, ["H9"])
    write_tsv(Path(args.a_candidates_dir) / "H9.candidate_events.tsv", [])
    audit, summary = build_reference_audit(args)

    assert audit.loc[0, "R_label"] == "R0"
    assert summary["formal_n0_used"] is False
    assert summary["legacy_branch_b_kept_counts_used_for_r_label"] is False


def test_shared_batch_branch_a_signal_can_stay_r0_context(tmp_path: Path):
    args = base_args(tmp_path, ["H9", "H10", "H11"])
    for sample_id in args.sample_id:
        write_tsv(
            Path(args.a_candidates_dir) / f"{sample_id}.candidate_events.tsv",
            [candidate(sample_id, "chr19", 1_000_000, 20_000_000)],
        )

    audit, _summary = build_reference_audit(args)

    assert set(audit["shared_batch_signal_flag"]) == {"yes"}
    assert set(audit["sample_specific_signal_flag"]) == {"no"}
    assert set(audit["R_label"]) == {"R0"}


def test_sample_specific_broad_branch_a_signal_is_r1(tmp_path: Path):
    args = base_args(tmp_path, ["H9", "H10", "H11"])
    write_tsv(
        Path(args.a_candidates_dir) / "H9.candidate_events.tsv",
        [candidate("H9", "chr5", 1_000_000, 30_000_000)],
    )
    write_tsv(Path(args.a_candidates_dir) / "H10.candidate_events.tsv", [])
    write_tsv(Path(args.a_candidates_dir) / "H11.candidate_events.tsv", [])

    audit, _summary = build_reference_audit(args)
    labels = dict(zip(audit["sample"], audit["R_label"]))

    assert labels["H9"] == "R1"
    assert labels["H10"] == "R0"
    assert labels["H11"] == "R0"


def test_small_acrocentric_candidate_is_not_automatic_high_risk(tmp_path: Path):
    args = base_args(tmp_path, ["H9"])
    write_tsv(
        Path(args.a_candidates_dir) / "H9.candidate_events.tsv",
        [candidate("H9", "chr21", 20_000_000, 22_000_000)],
    )

    audit, _summary = build_reference_audit(args)

    assert audit.loc[0, "acrocentric_or_high_repeat_burden"] < 0.20
    assert audit.loc[0, "R_label"] == "R0"


def test_reference_audit_writes_expected_outputs(tmp_path: Path):
    args = base_args(tmp_path, ["H9"])
    write_tsv(Path(args.a_candidates_dir) / "H9.candidate_events.tsv", [])

    audit, summary = build_reference_audit(args)
    Path(args.output_audit).parent.mkdir(parents=True, exist_ok=True)
    audit.to_csv(args.output_audit, sep="\t", index=False)
    Path(args.output_summary).write_text(json.dumps(summary), encoding="utf-8")

    assert Path(args.output_audit).exists()
    assert Path(args.output_summary).exists()
