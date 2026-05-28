from __future__ import annotations

import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

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
    from pgta.predict import benchmark, evaluation
    from pgta.predict.truth import event_detected, event_support_z, normalize_chrom, normalize_state


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

    def test_truth_helpers_are_shared_by_benchmark_and_evaluation(self):
        self.assertIs(benchmark.normalize_chrom, normalize_chrom)
        self.assertIs(evaluation.normalize_chrom, normalize_chrom)
        self.assertIs(benchmark.normalize_state, normalize_state)
        self.assertIs(evaluation.normalize_state, normalize_state)
        self.assertIs(benchmark.event_support_z, event_support_z)
        self.assertIs(evaluation.event_support_z, event_support_z)
        self.assertIs(benchmark.event_detected, event_detected)
        self.assertIs(evaluation.event_detected, event_detected)

    def test_evaluation_detects_low_level_broad_raw_gain_truth_match(self):
        events_df = pd.DataFrame(
            [
                {
                    "sample_id": "Y2",
                    "event_id": "Y2.chr14.raw_chr_dosage.gain",
                    "chrom": "chr14",
                    "state": "gain",
                    "caller": "raw_chromosome_dosage_detector",
                    "start": 0,
                    "end": 108_000_000,
                    "keep_event": 1,
                    "calibrated_mean_z": 1.86,
                    "calibrated_median_z": 1.86,
                    "event_corr_adjusted_z": 1.86,
                    "empirical_qvalue": 0.062,
                    "n_bins": 108,
                    "effective_bin_count": 86.0,
                    "clean_bin_fraction": 0.22,
                    "moderate_risk_bin_fraction": 0.77,
                    "high_risk_bin_fraction": 0.01,
                    "artifact_flags": "low_calibrated_signal,broad_chrom_fraction,edge_event",
                }
            ]
        )
        truth_df = pd.DataFrame(
            [
                {
                    "sample_id": "Y2",
                    "chrom": "14",
                    "expected_state": "gain",
                    "start": 0,
                    "end": 108_000_000,
                }
            ]
        )

        with TemporaryDirectory() as tmpdir:
            truth_tsv = Path(tmpdir) / "truth.tsv"
            truth_df.to_csv(truth_tsv, sep="\t", index=False)
            metrics_df, summary = evaluation.compute_truth_metrics(events_df, str(truth_tsv), branch_b_z_threshold=2.5)

        self.assertEqual(int(metrics_df.loc[0, "matched"]), 1)
        self.assertEqual(int(metrics_df.loc[0, "detected"]), 1)
        self.assertAlmostEqual(float(metrics_df.loc[0, "top_support_z"]), 1.86)
        self.assertEqual(summary["truth_detected_count"], 1)


if __name__ == "__main__":
    unittest.main()
