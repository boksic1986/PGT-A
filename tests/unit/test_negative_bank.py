from __future__ import annotations

import sys
import unittest
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from pgta.predict.branch_b.negative_bank import assign_negative_bank_labels, summarize_negative_bank


class NegativeBankLabelTest(unittest.TestCase):
    def test_n0_is_the_only_matched_negative_eligible_label(self):
        samples = pd.DataFrame(
            [
                {
                    "sample_id": "N_LOCKED",
                    "qc_status": "PASS",
                    "sample_role": "locked_negative",
                    "branch_b_kept_count": 0,
                    "manual_review_status": "clean",
                },
                {
                    "sample_id": "H9",
                    "qc_status": "PASS",
                    "sample_role": "presumed_negative",
                    "branch_b_kept_count": 0,
                    "manual_negative_bank_label": "N1",
                    "manual_review_status": "ref_evaluation_only",
                },
                {
                    "sample_id": "H8",
                    "qc_status": "PASS",
                    "sample_role": "presumed_negative",
                    "branch_b_kept_count": 1,
                    "manual_review_status": "retained_review_event",
                },
            ]
        )

        labeled = assign_negative_bank_labels(samples, version="branch_ab_v2_2026-06-18")
        by_sample = labeled.set_index("sample_id")

        self.assertEqual(by_sample.loc["N_LOCKED", "negative_bank_label"], "N0")
        self.assertEqual(int(by_sample.loc["N_LOCKED", "matched_negative_eligible"]), 1)
        self.assertEqual(by_sample.loc["H9", "negative_bank_label"], "N1")
        self.assertEqual(int(by_sample.loc["H9", "matched_negative_eligible"]), 0)
        self.assertEqual(by_sample.loc["H8", "negative_bank_label"], "N2")
        self.assertIn("branch_b_retained_review_event", by_sample.loc["H8", "negative_bank_reason"])

    def test_h7_h16_initial_documented_labels_are_not_promoted_to_n0(self):
        samples = pd.DataFrame(
            [
                {"sample_id": "H7", "qc_status": "PASS", "manual_negative_bank_label": "N2", "manual_review_status": "recurring_review_artifact"},
                {"sample_id": "H9", "qc_status": "PASS", "manual_negative_bank_label": "N1", "manual_review_status": "ref_evaluation_only"},
                {"sample_id": "H10", "qc_status": "PASS", "manual_negative_bank_label": "N1", "manual_review_status": "ref_evaluation_only"},
                {"sample_id": "H11", "qc_status": "PASS", "manual_negative_bank_label": "N1", "manual_review_status": "ref_evaluation_only"},
                {"sample_id": "H12", "qc_status": "PASS", "manual_negative_bank_label": "N1", "manual_review_status": "ref_evaluation_only"},
                {"sample_id": "H13", "qc_status": "PASS", "manual_negative_bank_label": "N2", "manual_review_status": "retained_review_event"},
                {"sample_id": "H14", "qc_status": "PASS", "manual_negative_bank_label": "N2", "manual_review_status": "retained_review_event"},
                {"sample_id": "H15", "qc_status": "PASS", "manual_negative_bank_label": "N1", "manual_review_status": "ref_evaluation_only"},
                {"sample_id": "H16", "qc_status": "PASS", "manual_negative_bank_label": "N1", "manual_review_status": "ref_evaluation_only"},
            ]
        )

        labeled = assign_negative_bank_labels(samples, version="branch_ab_v2_2026-06-18")

        self.assertEqual(set(labeled["negative_bank_label"]), {"N1", "N2"})
        self.assertEqual(int(labeled["matched_negative_eligible"].sum()), 0)
        self.assertEqual(set(labeled[labeled["negative_bank_label"].eq("N1")]["sample_id"]), {"H9", "H10", "H11", "H12", "H15", "H16"})
        self.assertEqual(set(labeled[labeled["negative_bank_label"].eq("N2")]["sample_id"]), {"H7", "H13", "H14"})

    def test_summary_counts_labels_and_matched_negative_pool(self):
        labeled = pd.DataFrame(
            [
                {"negative_bank_label": "N0", "matched_negative_eligible": 1},
                {"negative_bank_label": "N1", "matched_negative_eligible": 0},
                {"negative_bank_label": "N2", "matched_negative_eligible": 0},
            ]
        )

        summary = summarize_negative_bank(labeled, version="v1")

        self.assertEqual(summary["version"], "v1")
        self.assertEqual(summary["label_counts"]["N0"], 1)
        self.assertEqual(summary["matched_negative_eligible_count"], 1)


if __name__ == "__main__":
    unittest.main()
