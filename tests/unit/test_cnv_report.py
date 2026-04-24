from __future__ import annotations

import sys
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
    from pgta.predict.report import (
        format_biological_candidate_conclusion,
        format_technical_conclusion,
        summarize_branch_b_events,
    )


@unittest.skipIf(IMPORT_ERROR is not None, f"optional dependency missing: {IMPORT_ERROR}")
class CnvReportRankingTest(unittest.TestCase):
    def test_sex_review_event_is_suppressed_when_autosomal_kept_event_exists(self):
        events_df = pd.DataFrame(
            [
                {
                    "sample_id": "Y6",
                    "event_id": "event_sex",
                    "keep_event": 1,
                    "artifact_status": "review",
                    "priority_score": 12.5,
                    "n_bins": 20,
                    "chrom": "chrX",
                    "start": 11_000_000,
                    "end": 37_000_000,
                    "state": "loss",
                    "technical_confidence": "moderate",
                    "artifact_flags": "sex_chromosome_event",
                    "downgrade_reason": "",
                    "filter_reason": "",
                    "retain_reason": "",
                },
                {
                    "sample_id": "Y6",
                    "event_id": "event_auto",
                    "keep_event": 1,
                    "artifact_status": "review",
                    "priority_score": 9.2,
                    "n_bins": 12,
                    "chrom": "chr7",
                    "start": 63_000_001,
                    "end": 159_000_000,
                    "state": "gain",
                    "technical_confidence": "moderate",
                    "artifact_flags": "low_level_signal",
                    "downgrade_reason": "",
                    "filter_reason": "",
                    "retain_reason": "",
                },
            ]
        )

        sample_df, top_branch_b = summarize_branch_b_events(events_df)
        merged = sample_df.merge(top_branch_b, on="sample_id", how="left")
        row = merged.iloc[0]

        self.assertEqual(int(row["branch_b_suppressed_sex_review_events"]), 1)
        self.assertEqual(row["branch_b_top_event"], "chr7:63000001-159000000 gain [review/moderate]")

    def test_only_sex_review_event_leaves_branch_b_top_event_blank(self):
        events_df = pd.DataFrame(
            [
                {
                    "sample_id": "Y3",
                    "event_id": "event_sex",
                    "keep_event": 1,
                    "artifact_status": "review",
                    "priority_score": 9.5,
                    "n_bins": 18,
                    "chrom": "chrX",
                    "start": 11_000_000,
                    "end": 37_000_000,
                    "state": "loss",
                    "technical_confidence": "moderate",
                    "artifact_flags": "sex_chromosome_event",
                    "downgrade_reason": "",
                    "filter_reason": "",
                    "retain_reason": "",
                }
            ]
        )

        sample_df, top_branch_b = summarize_branch_b_events(events_df)
        merged = sample_df.merge(top_branch_b, on="sample_id", how="left")
        row = merged.iloc[0]

        self.assertEqual(int(row["branch_b_suppressed_sex_review_events"]), 1)
        self.assertTrue(pd.isna(row["branch_b_top_event"]))

        conclusion_row = row.to_dict()
        conclusion_row["a_branch_top_event"] = ""
        technical = format_technical_conclusion(conclusion_row)
        biological = format_biological_candidate_conclusion(conclusion_row)

        self.assertIn("suppressed_sex_chromosome_review", technical)
        self.assertEqual(biological, "Branch B top event suppressed (sex-chromosome review only)")

    def test_pass_sex_chromosome_event_remains_displayable(self):
        events_df = pd.DataFrame(
            [
                {
                    "sample_id": "S1",
                    "event_id": "event_sex_pass",
                    "keep_event": 1,
                    "artifact_status": "pass",
                    "priority_score": 15.0,
                    "n_bins": 30,
                    "chrom": "chrX",
                    "start": 5_000_000,
                    "end": 40_000_000,
                    "state": "loss",
                    "technical_confidence": "high",
                    "artifact_flags": "sex_chromosome_event",
                    "downgrade_reason": "",
                    "filter_reason": "",
                    "retain_reason": "clean_support_and_statistical_support",
                }
            ]
        )

        sample_df, top_branch_b = summarize_branch_b_events(events_df)
        merged = sample_df.merge(top_branch_b, on="sample_id", how="left")
        row = merged.iloc[0]

        self.assertEqual(int(row["branch_b_suppressed_sex_review_events"]), 0)
        self.assertEqual(row["branch_b_top_event"], "chrX:5000000-40000000 loss [pass/high]")


if __name__ == "__main__":
    unittest.main()
