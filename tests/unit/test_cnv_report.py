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
        build_plot_lookup,
        format_biological_candidate_conclusion,
        format_technical_conclusion,
        summarize_branch_b_events,
    )


@unittest.skipIf(IMPORT_ERROR is not None, f"optional dependency missing: {IMPORT_ERROR}")
class CnvReportRankingTest(unittest.TestCase):
    def test_build_plot_lookup_indexes_existing_svg_by_sample_id(self):
        lookup = build_plot_lookup(
            [
                "/data/project/CNV/PGT-A/wisecondorx/cnv/plots/Y1.final_cnv.svg",
                "/data/project/CNV/PGT-A/wisecondorx/cnv/plots/Y2.branch_ab.svg",
            ]
        )

        self.assertEqual(
            lookup,
            {
                "Y1": "/data/project/CNV/PGT-A/wisecondorx/cnv/plots/Y1.final_cnv.svg",
                "Y2": "/data/project/CNV/PGT-A/wisecondorx/cnv/plots/Y2.branch_ab.svg",
            },
        )

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

    def test_only_sex_review_event_remains_displayable(self):
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

        self.assertEqual(int(row["branch_b_suppressed_sex_review_events"]), 0)
        self.assertEqual(row["branch_b_top_event"], "chrX:11000000-37000000 loss [review/moderate]")

        conclusion_row = row.to_dict()
        conclusion_row["a_branch_top_event"] = ""
        technical = format_technical_conclusion(conclusion_row)
        biological = format_biological_candidate_conclusion(conclusion_row)

        self.assertIn("chrX:11000000-37000000 loss", technical)
        self.assertEqual(biological, "chrX:11000000-37000000 loss [review/moderate]")

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

    def test_cnvseq_reportable_tier_is_counted_and_prioritized_for_top_event(self):
        events_df = pd.DataFrame(
            [
                {
                    "sample_id": "H1",
                    "event_id": "small_high_priority",
                    "keep_event": 1,
                    "artifact_status": "review",
                    "priority_score": 100.0,
                    "n_bins": 2,
                    "chrom": "chr8",
                    "start": 10_000_000,
                    "end": 11_500_000,
                    "state": "loss",
                    "technical_confidence": "moderate",
                    "artifact_flags": "review_1_2mb",
                    "downgrade_reason": "cnvseq_subreportable_size_review",
                    "filter_reason": "",
                    "retain_reason": "",
                    "cnvseq_report_tier": "review_1_2mb",
                    "cnvseq_reportable": 0,
                },
                {
                    "sample_id": "H1",
                    "event_id": "large_lower_priority",
                    "keep_event": 1,
                    "artifact_status": "review",
                    "priority_score": 10.0,
                    "n_bins": 4,
                    "chrom": "chr16",
                    "start": 46_000_001,
                    "end": 90_000_000,
                    "state": "gain",
                    "technical_confidence": "moderate",
                    "artifact_flags": "",
                    "downgrade_reason": "",
                    "filter_reason": "",
                    "retain_reason": "",
                    "cnvseq_report_tier": "reportable",
                    "cnvseq_reportable": 1,
                },
                {
                    "sample_id": "H1",
                    "event_id": "artifact",
                    "keep_event": 0,
                    "artifact_status": "artifact",
                    "priority_score": 999.0,
                    "n_bins": 20,
                    "chrom": "chr19",
                    "start": 1,
                    "end": 3_000_000,
                    "state": "gain",
                    "technical_confidence": "rejected",
                    "artifact_flags": "",
                    "downgrade_reason": "",
                    "filter_reason": "",
                    "retain_reason": "",
                    "cnvseq_report_tier": "reportable",
                    "cnvseq_reportable": 1,
                },
            ]
        )

        sample_df, top_branch_b = summarize_branch_b_events(events_df)
        merged = sample_df.merge(top_branch_b, on="sample_id", how="left")
        row = merged.iloc[0]

        self.assertEqual(int(row["branch_b_reportable_events"]), 1)
        self.assertEqual(int(row["branch_b_review_tier_events"]), 1)
        self.assertEqual(row["branch_b_top_event"], "chr16:46000001-90000000 gain [review/moderate; reportable]")


if __name__ == "__main__":
    unittest.main()
