from __future__ import annotations

import sys
import unittest
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from pgta.predict.branch_b.matched_negative import (
    build_matched_negative_percentiles,
    summarize_matched_negative,
)


class BranchBMatchedNegativeTest(unittest.TestCase):
    def test_same_region_n0_background_adds_shadow_percentile(self):
        query = pd.DataFrame(
            [
                {
                    "sample_id": "H1",
                    "candidate_id": "H1.A0001",
                    "chrom": "chr21",
                    "start": 20_000_000,
                    "end": 22_000_000,
                    "state": "gain",
                    "corrected_amplitude": 4.0,
                    "final_disposition": "REVIEW_REQUIRED",
                }
            ]
        )
        background = pd.DataFrame(
            [
                {
                    "sample_id": "N0_A",
                    "candidate_id": "N0_A.A0001",
                    "chrom": "chr21",
                    "start": 20_000_000,
                    "end": 22_000_000,
                    "state": "gain",
                    "corrected_amplitude": 1.0,
                },
                {
                    "sample_id": "N0_B",
                    "candidate_id": "N0_B.A0001",
                    "chrom": "chr21",
                    "start": 20_000_000,
                    "end": 22_000_000,
                    "state": "gain",
                    "corrected_amplitude": 3.0,
                },
                {
                    "sample_id": "N1_C",
                    "candidate_id": "N1_C.A0001",
                    "chrom": "chr21",
                    "start": 20_000_000,
                    "end": 22_000_000,
                    "state": "gain",
                    "corrected_amplitude": 100.0,
                },
            ]
        )
        labels = pd.DataFrame(
            [
                {"sample_id": "N0_A", "negative_bank_label": "N0", "matched_negative_eligible": 1},
                {"sample_id": "N0_B", "negative_bank_label": "N0", "matched_negative_eligible": 1},
                {"sample_id": "N1_C", "negative_bank_label": "N1", "matched_negative_eligible": 0},
            ]
        )

        result = build_matched_negative_percentiles(
            query,
            background,
            labels,
            min_background=2,
            feature_column="corrected_amplitude",
        )
        row = result.iloc[0]

        self.assertEqual(row["matched_negative_background_status"], "OK")
        self.assertEqual(row["matched_negative_scope"], "same_region")
        self.assertEqual(int(row["matched_negative_n"]), 2)
        self.assertAlmostEqual(float(row["matched_negative_abs_percentile"]), 1.0)
        self.assertEqual(row["matched_negative_action"], "BACKGROUND_SUPPORTED")

    def test_no_n0_background_is_unknown_review_not_artifact(self):
        query = pd.DataFrame(
            [
                {
                    "sample_id": "H1",
                    "candidate_id": "H1.A0002",
                    "chrom": "chr7",
                    "start": 10_000_000,
                    "end": 12_000_000,
                    "state": "loss",
                    "corrected_amplitude": -2.5,
                    "final_disposition": "REVIEW_REQUIRED",
                }
            ]
        )
        labels = pd.DataFrame(
            [
                {"sample_id": "H9", "negative_bank_label": "N1", "matched_negative_eligible": 0},
                {"sample_id": "H13", "negative_bank_label": "N2", "matched_negative_eligible": 0},
            ]
        )

        result = build_matched_negative_percentiles(
            query,
            pd.DataFrame(),
            labels,
            min_background=2,
            feature_column="corrected_amplitude",
        )
        row = result.iloc[0]

        self.assertEqual(row["matched_negative_background_status"], "UNKNOWN_BACKGROUND")
        self.assertEqual(row["matched_negative_scope"], "UNKNOWN_BACKGROUND")
        self.assertEqual(int(row["matched_negative_n"]), 0)
        self.assertTrue(pd.isna(row["matched_negative_abs_percentile"]))
        self.assertEqual(row["matched_negative_action"], "REVIEW_NO_CALL")
        self.assertEqual(row["final_disposition"], "REVIEW_REQUIRED")

    def test_summary_counts_unknown_background_without_promotion(self):
        annotated = pd.DataFrame(
            [
                {"matched_negative_background_status": "OK", "matched_negative_action": "BACKGROUND_SUPPORTED"},
                {"matched_negative_background_status": "UNKNOWN_BACKGROUND", "matched_negative_action": "REVIEW_NO_CALL"},
            ]
        )

        summary = summarize_matched_negative("H1", annotated, version="phase3_test")

        self.assertEqual(summary["sample_id"], "H1")
        self.assertEqual(summary["version"], "phase3_test")
        self.assertEqual(summary["candidate_count"], 2)
        self.assertEqual(summary["background_status_counts"]["UNKNOWN_BACKGROUND"], 1)
        self.assertEqual(summary["final_report_impact"], "none_shadow_only")


if __name__ == "__main__":
    unittest.main()
