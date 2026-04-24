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
    from pgta.predict.benchmark import overlap_fraction, select_matches
    from pgta.predict.evaluation import overlap_fraction as evaluation_overlap_fraction


@unittest.skipIf(IMPORT_ERROR is not None, f"optional dependency missing: {IMPORT_ERROR}")
class BenchmarkOverlapTest(unittest.TestCase):
    def test_overlap_fraction_treats_shared_boundary_as_positive_overlap(self):
        self.assertGreater(overlap_fraction(1_000_000, 58_000_000, 58_000_000, 61_000_000), 0.0)

    def test_select_matches_accepts_shared_boundary_but_not_one_base_gap(self):
        events_df = pd.DataFrame(
            [
                {
                    "sample_id": "Y6",
                    "chrom": "chr7",
                    "state": "loss",
                    "start": 58_000_000,
                    "end": 61_000_000,
                },
                {
                    "sample_id": "Y6",
                    "chrom": "chr7",
                    "state": "loss",
                    "start": 58_000_001,
                    "end": 61_000_000,
                },
            ]
        )

        shared_boundary = select_matches(events_df, "Y6", "chr7", "loss", 1_000_000, 58_000_000, "start", "end")
        one_base_gap = select_matches(
            events_df.iloc[[1]].copy(), "Y6", "chr7", "loss", 1_000_000, 58_000_000, "start", "end"
        )

        self.assertEqual(len(shared_boundary), 1)
        self.assertEqual(int(shared_boundary.iloc[0]["start"]), 58_000_000)
        self.assertTrue(one_base_gap.empty)

    def test_evaluation_overlap_fraction_uses_same_boundary_rule(self):
        self.assertGreater(evaluation_overlap_fraction(1_000_000, 58_000_000, 58_000_000, 61_000_000), 0.0)
        self.assertEqual(evaluation_overlap_fraction(1_000_000, 58_000_000, 58_000_001, 61_000_000), 0.0)


if __name__ == "__main__":
    unittest.main()
