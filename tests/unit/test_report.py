from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

try:
    import pandas as pd
except ModuleNotFoundError as exc:  # pragma: no cover - environment-dependent skip
    pd = None
    IMPORT_ERROR = exc
else:
    IMPORT_ERROR = None
    from pgta.predict.report import format_biological_candidate_conclusion, format_technical_conclusion, load_a_branch


@unittest.skipIf(IMPORT_ERROR is not None, f"optional dependency missing: {IMPORT_ERROR}")
class ReportBranchRoleTest(unittest.TestCase):
    def test_branch_b_event_remains_final_candidate_when_present(self):
        row = pd.Series(
            {
                "branch_b_top_event": "chr21:9000000-48000000 gain [pass/high]",
                "a_branch_top_event": "chr21:9000001-48000000 gain",
            }
        )

        self.assertEqual(
            format_biological_candidate_conclusion(row),
            "chr21:9000000-48000000 gain [pass/high]",
        )

    def test_a_branch_only_signal_is_reported_as_sensitive_review_candidate(self):
        row = pd.Series(
            {
                "branch_b_top_event": "",
                "a_branch_top_event": "chr16:46000001-90000000 gain",
                "branch_b_suppressed_sex_review_events": 0,
            }
        )

        self.assertEqual(
            format_biological_candidate_conclusion(row),
            "A-branch sensitive signal only: chr16:46000001-90000000 gain; requires Branch B review",
        )

    def test_technical_conclusion_includes_a_branch_sensitive_signal_summary(self):
        row = pd.Series(
            {
                "branch_b_kept_events": 2,
                "branch_b_top_event": "chr7:63000000-159000000 gain [review/moderate]",
                "branch_b_top_fraction": "11.3%",
                "branch_b_top_flags": "clean_event_below_ultra_pass",
                "branch_b_top_downgrade_reason": "low_risk_event_preserved_below_ultra_pass",
                "branch_b_suppressed_sex_review_events": 0,
                "a_branch_event_count": 12,
                "a_branch_top_event": "chr7:63000001-159000000 gain",
                "a_branch_review_candidate_count": 3,
                "a_branch_review_shortlist": "chr7:63000001-159000000 gain z=12.40",
                "a_branch_strong_signal_count": 1,
                "qc_status": "PASS",
                "sex_call": "XY",
            }
        )

        conclusion = format_technical_conclusion(row)

        self.assertIn("Branch B kept 2 events", conclusion)
        self.assertIn("A-branch_sensitive_candidates=12", conclusion)
        self.assertIn("A-branch_review_shortlist_top3=3", conclusion)
        self.assertIn("A-branch_strong_signals_z10=1", conclusion)
        self.assertIn("A-branch_top_signal=chr7:63000001-159000000 gain", conclusion)

    def test_a_branch_loader_merges_and_ranks_review_shortlist(self):
        with tempfile.TemporaryDirectory() as temp_dir_value:
            temp_dir = Path(temp_dir_value)
            bed_path = temp_dir / "S1_aberrations.bed"
            pd.DataFrame(
                [
                    {"chr": "7", "start": 1, "end": 100, "type": "gain", "zscore": 6.0},
                    {"chr": "7", "start": 101, "end": 200, "type": "gain", "zscore": 12.0},
                    {"chr": "8", "start": 1, "end": 150, "type": "loss", "zscore": -8.0},
                    {"chr": "9", "start": 1, "end": 160, "type": "gain", "zscore": 11.0},
                    {"chr": "10", "start": 1, "end": 170, "type": "loss", "zscore": -7.0},
                ]
            ).to_csv(bed_path, sep="\t", index=False)
            a_branch_df = load_a_branch([str(bed_path)])

        row = a_branch_df.iloc[0]
        self.assertEqual(int(row["a_branch_event_count"]), 4)
        self.assertEqual(row["a_branch_top_event"], "chr7:1-200 gain")
        self.assertEqual(int(row["a_branch_review_candidate_count"]), 3)
        self.assertEqual(int(row["a_branch_strong_signal_count"]), 2)
        self.assertIn("chr7:1-200 gain z=12.00", row["a_branch_review_shortlist"])
        self.assertIn("chr9:1-160 gain z=11.00", row["a_branch_review_shortlist"])


if __name__ == "__main__":
    unittest.main()
