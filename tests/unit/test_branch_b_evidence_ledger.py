from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from pgta.predict.branch_b.evidence_ledger import build_candidate_evidence_ledger, summarize_evidence_ledger


class BranchBEvidenceLedgerTest(unittest.TestCase):
    def test_every_a_candidate_gets_one_shadow_evidence_row(self):
        a_candidates = pd.DataFrame(
            [
                {
                    "candidate_id": "H1.A0001",
                    "sample_id": "H1",
                    "chrom": "chr21",
                    "start": 20_000_000,
                    "end": 22_000_000,
                    "state": "gain",
                    "a_zscore": 55.0,
                    "a_abs_zscore": 55.0,
                    "a_ratio": 0.11,
                    "a_support_level": "strong",
                },
                {
                    "candidate_id": "H1.A0002",
                    "sample_id": "H1",
                    "chrom": "chr7",
                    "start": 10_000_000,
                    "end": 12_000_000,
                    "state": "loss",
                    "a_zscore": -8.0,
                    "a_abs_zscore": 8.0,
                    "a_ratio": -0.04,
                    "a_support_level": "sensitive",
                },
            ]
        )
        branch_b_events = pd.DataFrame(
            [
                {
                    "event_id": "H1.H1.A0001.B_refined",
                    "a_candidate_id": "H1.A0001",
                    "sample_id": "H1",
                    "chrom": "chr21",
                    "state": "gain",
                    "keep_event": 1,
                    "report_class": "candidate_pass",
                    "artifact_status": "pass",
                    "segment_mean_robust_z": 12.0,
                    "calibrated_mean_z": 10.0,
                    "event_corr_adjusted_z": 9.0,
                    "clean_bin_fraction": 0.80,
                    "high_risk_bin_fraction": 0.10,
                    "wisecondorx_ref_bin_count": 180.0,
                    "low_refbin_fraction": 0.0,
                    "same_chrom_ref_bin_count": 0.0,
                }
            ]
        )
        bins = pd.DataFrame(
            {
                "chrom": ["chr21", "chr21", "chr7", "chr7"],
                "start": [20_000_000, 21_000_000, 10_000_000, 11_000_000],
                "end": [21_000_000, 22_000_000, 11_000_000, 12_000_000],
                "bin_index": [20, 21, 10, 11],
                "calibrated_z": [8.5, 9.5, -1.0, -1.2],
                "robust_z": [10.0, 11.0, -1.1, -1.3],
                "mask_label": ["pass", "pass", "pass", "pass"],
                "region_risk_class": ["clean", "clean", "clean", "clean"],
                "calibration_null_eligible": [1, 1, 1, 1],
                "wisecondorx_ref_bin_count": [180.0, 181.0, np.nan, np.nan],
                "low_refbin_fraction": [0.0, 0.0, np.nan, np.nan],
                "same_chrom_ref_bin_count": [0.0, 0.0, np.nan, np.nan],
            }
        )

        ledger = build_candidate_evidence_ledger(
            sample_id="H1",
            a_candidates=a_candidates,
            branch_b_events=branch_b_events,
            bins_df=bins,
            reference_id="ref_mask_only_100kb",
            negative_bank_version="negative_bank_v1",
        )

        self.assertEqual(len(ledger), 2)
        self.assertEqual(set(ledger["candidate_id"]), {"H1.A0001", "H1.A0002"})

        confirmed = ledger.set_index("candidate_id").loc["H1.A0001"]
        self.assertEqual(confirmed["final_disposition"], "CONFIRMED")
        self.assertEqual(confirmed["refmap_status"], "OK")
        self.assertEqual(confirmed["calibration_null_status"], "OK")
        self.assertGreater(float(confirmed["same_direction_fraction"]), 0.9)

        missing = ledger.set_index("candidate_id").loc["H1.A0002"]
        self.assertEqual(missing["final_disposition"], "REVIEW_REQUIRED")
        self.assertEqual(missing["refmap_status"], "UNKNOWN")
        self.assertEqual(missing["calibration_null_status"], "OK")
        self.assertIn("branch_b_event_missing", missing["evidence_missing_reason"])
        self.assertIn("refmap_missing", missing["evidence_missing_reason"])

    def test_unknown_refmap_and_calibration_null_are_not_clean_support(self):
        a_candidates = pd.DataFrame(
            [
                {
                    "candidate_id": "S1.A0001",
                    "sample_id": "S1",
                    "chrom": "chr1",
                    "start": 0,
                    "end": 2_000_000,
                    "state": "gain",
                    "a_zscore": 12.0,
                    "a_abs_zscore": 12.0,
                    "a_ratio": 0.05,
                    "a_support_level": "strong",
                }
            ]
        )
        branch_b_events = pd.DataFrame(
            [
                {
                    "event_id": "S1.A0001.B_refined",
                    "a_candidate_id": "S1.A0001",
                    "sample_id": "S1",
                    "chrom": "chr1",
                    "state": "gain",
                    "keep_event": 1,
                    "report_class": "candidate_review",
                    "artifact_status": "review",
                    "segment_mean_robust_z": 3.0,
                    "calibrated_mean_z": 2.5,
                }
            ]
        )

        ledger = build_candidate_evidence_ledger(
            sample_id="S1",
            a_candidates=a_candidates,
            branch_b_events=branch_b_events,
            bins_df=pd.DataFrame(),
            reference_id="ref1",
            negative_bank_version="nb1",
        )

        row = ledger.iloc[0]
        self.assertEqual(row["final_disposition"], "REVIEW_REQUIRED")
        self.assertEqual(row["refmap_status"], "UNKNOWN")
        self.assertEqual(row["calibration_null_status"], "UNKNOWN")
        self.assertNotEqual(row["refmap_status"], "OK")
        self.assertIn("calibration_null_missing", row["evidence_missing_reason"])

    def test_summary_reports_disposition_and_missingness_counts(self):
        ledger = pd.DataFrame(
            [
                {"final_disposition": "CONFIRMED", "evidence_missing_reason": ""},
                {"final_disposition": "REVIEW_REQUIRED", "evidence_missing_reason": "refmap_missing"},
            ]
        )

        summary = summarize_evidence_ledger("S1", ledger)

        self.assertEqual(summary["candidate_count"], 2)
        self.assertEqual(summary["disposition_counts"]["CONFIRMED"], 1)
        self.assertEqual(summary["disposition_counts"]["REVIEW_REQUIRED"], 1)
        self.assertEqual(summary["missing_evidence_candidate_count"], 1)


if __name__ == "__main__":
    unittest.main()
