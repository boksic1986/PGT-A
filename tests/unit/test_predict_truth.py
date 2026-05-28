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
    from pgta.predict.truth import event_detected, event_detection_threshold

from pgta.predict.truth import normalize_chrom, normalize_state, overlap_fraction


class PredictTruthHelpersTest(unittest.TestCase):
    def test_normalizes_chromosome_and_state_aliases(self):
        self.assertEqual(normalize_chrom("7"), "chr7")
        self.assertEqual(normalize_chrom("chrX"), "chrX")
        self.assertEqual(normalize_state("duplication"), "gain")
        self.assertEqual(normalize_state("DEL"), "loss")

    def test_overlap_fraction_treats_shared_boundary_as_overlap(self):
        self.assertGreater(overlap_fraction(1_000_000, 58_000_000, 58_000_000, 61_000_000), 0.0)
        self.assertEqual(overlap_fraction(1_000_000, 58_000_000, 58_000_001, 61_000_000), 0.0)


@unittest.skipIf(IMPORT_ERROR is not None, f"optional dependency missing: {IMPORT_ERROR}")
class PredictTruthDetectionTest(unittest.TestCase):
    def test_low_level_broad_raw_autosome_gain_uses_rescue_threshold(self):
        events = pd.DataFrame(
            [
                {
                    "sample_id": "Y2",
                    "chrom": "chr14",
                    "state": "gain",
                    "caller": "raw_chromosome_dosage_detector",
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

        self.assertEqual(float(event_detection_threshold(events, standard_z_threshold=2.5).iloc[0]), 1.8)
        self.assertTrue(bool(event_detected(events, standard_z_threshold=2.5).iloc[0]))

    def test_low_level_focal_or_weak_qvalue_event_keeps_standard_threshold(self):
        events = pd.DataFrame(
            [
                {
                    "chrom": "chr14",
                    "state": "gain",
                    "caller": "segment_level_detector",
                    "keep_event": 1,
                    "calibrated_mean_z": 1.86,
                    "calibrated_median_z": 1.86,
                    "event_corr_adjusted_z": 1.86,
                    "empirical_qvalue": 0.062,
                    "effective_bin_count": 86.0,
                    "clean_bin_fraction": 0.22,
                    "moderate_risk_bin_fraction": 0.77,
                    "high_risk_bin_fraction": 0.01,
                    "artifact_flags": "edge_event",
                },
                {
                    "chrom": "chr21",
                    "state": "gain",
                    "caller": "raw_chromosome_dosage_detector",
                    "keep_event": 1,
                    "calibrated_mean_z": 1.86,
                    "calibrated_median_z": 1.86,
                    "event_corr_adjusted_z": 1.86,
                    "empirical_qvalue": 0.20,
                    "effective_bin_count": 33.0,
                    "clean_bin_fraction": 0.30,
                    "moderate_risk_bin_fraction": 0.66,
                    "high_risk_bin_fraction": 0.04,
                    "artifact_flags": "broad_chrom_fraction,edge_event",
                },
            ]
        )

        self.assertEqual(event_detection_threshold(events, standard_z_threshold=2.5).tolist(), [2.5, 2.5])
        self.assertEqual(event_detected(events, standard_z_threshold=2.5).tolist(), [False, False])

    def test_a_refined_candidate_uses_wisecondorx_a_support_for_detection(self):
        events = pd.DataFrame(
            [
                {
                    "chrom": "chr21",
                    "state": "gain",
                    "caller": "wisecondorx_a_branch",
                    "keep_event": 1,
                    "calibrated_mean_z": 1.1,
                    "calibrated_median_z": 1.0,
                    "event_corr_adjusted_z": 0.9,
                    "a_abs_zscore": 83.9,
                }
            ]
        )

        self.assertTrue(bool(event_detected(events, standard_z_threshold=2.5).iloc[0]))


if __name__ == "__main__":
    unittest.main()
