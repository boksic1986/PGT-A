from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from pgta.predict.branch_a_validation import build_branch_a_validation


def write_tsv(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def test_branch_a_validation_reports_truth_recall_and_h6_chr21_status(tmp_path: Path):
    y1_candidates = tmp_path / "a_branch" / "Y1.candidate_events.tsv"
    h6_candidates = tmp_path / "a_branch" / "H6.candidate_events.tsv"
    h7_candidates = tmp_path / "a_branch" / "H7.candidate_events.tsv"
    truth_tsv = tmp_path / "truth.tsv"
    out_samples = tmp_path / "branch_a_validation" / "sample_summary.tsv"
    out_truth = tmp_path / "branch_a_validation" / "truth_metrics.tsv"
    out_summary = tmp_path / "branch_a_validation" / "summary.json"

    write_tsv(
        y1_candidates,
        [
            {
                "candidate_id": "Y1.A0001",
                "sample_id": "Y1",
                "chrom": "chr21",
                "start": 100,
                "end": 500,
                "state": "gain",
                "a_abs_zscore": 12.0,
                "a_zscore": 12.0,
                "a_support_level": "strong",
            }
        ],
    )
    write_tsv(
        h6_candidates,
        [
            {
                "candidate_id": "H6.A0001",
                "sample_id": "H6",
                "chrom": "chr21",
                "start": 200,
                "end": 600,
                "state": "gain",
                "a_abs_zscore": 7.11,
                "a_zscore": 7.11,
                "a_support_level": "sensitive",
            }
        ],
    )
    write_tsv(
        h7_candidates,
        [
            {
                "candidate_id": "H7.A0001",
                "sample_id": "H7",
                "chrom": "chr1",
                "start": 1,
                "end": 100,
                "state": "loss",
                "a_abs_zscore": 5.5,
                "a_zscore": -5.5,
                "a_support_level": "sensitive",
            }
        ],
    )
    write_tsv(
        truth_tsv,
        [
            {"sample_id": "Y1", "chrom": "21", "start": 150, "end": 450, "expected_state": "gain"},
            {"sample_id": "H6", "chrom": "chr21", "start": 250, "end": 550, "expected_state": "gain"},
        ],
    )

    sample_summary, truth_metrics, summary = build_branch_a_validation(
        candidate_paths=[str(y1_candidates), str(h6_candidates), str(h7_candidates)],
        truth_tsv=str(truth_tsv),
        sample_ids=["Y1", "H6", "H7"],
        reference_id="h_r0_shadow_ref_20260619",
        binsize=100000,
        preprocess_mask_version="mask_only",
        wisecondorx_predict_command="WisecondorX predict ... --blacklist hard_mask.bed",
        blacklist_bed="hard_mask.bed",
        minrefbins=150,
        zscore=8.0,
        alpha=0.001,
        maskrepeats=5,
        output_sample_summary=str(out_samples),
        output_truth_metrics=str(out_truth),
        output_summary=str(out_summary),
    )

    assert out_samples.exists()
    assert out_truth.exists()
    assert out_summary.exists()
    assert json.loads(out_summary.read_text())["FN_count"] == 0

    y1 = sample_summary.set_index("sample").loc["Y1"]
    assert int(y1["truth_event_count"]) == 1
    assert int(y1["truth_detected_count"]) == 1
    assert int(y1["FN_count"]) == 0
    assert int(y1["branch_a_strong_count"]) == 1

    h6 = sample_summary.set_index("sample").loc["H6"]
    assert h6["H6_chr21_status"] == "detected"
    assert float(h6["top_branch_a_abs_z"]) == 7.11

    h7 = sample_summary.set_index("sample").loc["H7"]
    assert int(h7["truth_event_count"]) == 0
    assert int(h7["branch_a_candidate_count"]) == 1
    assert h7["reference_id"] == "h_r0_shadow_ref_20260619"

    assert set(truth_metrics["detected"].tolist()) == {1}
