from __future__ import annotations

import sys
import types
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

try:
    from pgta.predict.branch_b.artifact_rules import classify_event
except ModuleNotFoundError as exc:  # pragma: no cover - environment-dependent skip
    classify_event = None
    IMPORT_ERROR = exc
else:
    IMPORT_ERROR = None


@unittest.skipIf(IMPORT_ERROR is not None, f"optional dependency missing: {IMPORT_ERROR}")
class BranchBArtifactRulesTest(unittest.TestCase):
    def build_args(self):
        return types.SimpleNamespace(
            min_event_bins=3,
            min_abs_calibrated_z=2.0,
            max_chrom_fraction=0.35,
            edge_bin_window=2,
            max_qvalue=0.25,
            keep_review=1,
            high_confidence_z=4.0,
            high_confidence_qvalue=0.05,
            broad_support_min_abs_z=4.0,
            broad_support_max_qvalue=0.25,
            broad_support_min_clean_fraction=0.30,
            broad_support_min_effective_bins=10.0,
            edge_review_min_priority=2.0,
            ultra_pass_z=15.0,
            ultra_pass_qvalue=0.001,
            ultra_pass_effective_bins=8.0,
            clean_review_min_support_fraction=0.50,
            clean_review_max_overlap_fraction=0.15,
            clean_review_max_region_risk=0.35,
            focal_review_min_support_z=6.0,
            focal_review_max_overlap_fraction=0.25,
            focal_review_max_region_risk=0.20,
        )

    def test_high_confidence_broad_event_is_preserved_for_review(self):
        row = types.SimpleNamespace(
            chrom="chr7",
            start=1,
            end=159000000,
            start_bin=0,
            end_bin=158,
            n_bins=159,
            calibrated_mean_z=6.2,
            calibrated_median_z=6.0,
            event_corr_adjusted_z=5.1,
            empirical_qvalue=0.02,
            clean_bin_fraction=0.92,
            high_risk_bin_fraction=0.0,
            effective_bin_count=140.0,
            region_risk_score_mean=0.0,
            region_risk_score_max=0.0,
            xtr_overlap_fraction=0.0,
            sex_homology_overlap_fraction=0.0,
            segmental_duplication_overlap_fraction=0.0,
            low_mappability_overlap_fraction=0.0,
            gap_centromere_telomere_overlap_fraction=0.0,
            repeat_rich_overlap_fraction=0.0,
            blacklist_overlap_fraction=0.0,
            ambiguous_alignment_overlap_fraction=0.0,
            high_risk_boundary_crossing=0,
            caller="weighted_cbs_hmm",
            segment_mean_robust_z=-6.3,
            segment_median_robust_z=-6.0,
            segment_abs_max_robust_z=7.1,
        )

        result = classify_event(
            row=row,
            chrom_bin_count=159,
            args=self.build_args(),
            sex_call="XY",
            par_regions={},
        )

        self.assertEqual(result["artifact_status"], "review")
        self.assertEqual(result["keep_event"], 1)
        self.assertIn("broad_event_preserved_by_internal_support", result["downgrade_reason"])

    def test_clean_non_ultra_event_with_internal_support_is_review_not_artifact(self):
        row = types.SimpleNamespace(
            chrom="chr15",
            start=23000000,
            end=29000000,
            start_bin=23,
            end_bin=29,
            n_bins=7,
            calibrated_mean_z=5.5,
            calibrated_median_z=5.2,
            event_corr_adjusted_z=4.8,
            empirical_qvalue=0.04,
            clean_bin_fraction=1.0,
            high_risk_bin_fraction=0.0,
            effective_bin_count=12.0,
            region_risk_score_mean=0.0,
            region_risk_score_max=0.0,
            xtr_overlap_fraction=0.0,
            sex_homology_overlap_fraction=0.0,
            segmental_duplication_overlap_fraction=0.0,
            low_mappability_overlap_fraction=0.0,
            gap_centromere_telomere_overlap_fraction=0.0,
            repeat_rich_overlap_fraction=0.0,
            blacklist_overlap_fraction=0.0,
            ambiguous_alignment_overlap_fraction=0.0,
            high_risk_boundary_crossing=0,
            caller="weighted_cbs_hmm",
            segment_mean_robust_z=-5.1,
            segment_median_robust_z=-5.0,
            segment_abs_max_robust_z=5.4,
        )

        result = classify_event(
            row=row,
            chrom_bin_count=120,
            args=self.build_args(),
            sex_call="XX",
            par_regions={},
        )

        self.assertEqual(result["artifact_status"], "review")
        self.assertEqual(result["keep_event"], 1)
        self.assertIn("low_risk_event_preserved_below_ultra_pass", result["downgrade_reason"])

    def test_moderate_only_non_high_risk_event_with_strong_signal_is_review(self):
        row = types.SimpleNamespace(
            chrom="chr21",
            start=20000000,
            end=27000000,
            start_bin=20,
            end_bin=26,
            n_bins=7,
            calibrated_mean_z=2.49,
            calibrated_median_z=2.45,
            event_corr_adjusted_z=6.60,
            empirical_qvalue=0.04,
            clean_bin_fraction=0.42,
            moderate_risk_bin_fraction=0.58,
            high_risk_bin_fraction=0.0,
            effective_bin_count=12.0,
            region_risk_score_mean=0.02,
            region_risk_score_max=0.08,
            xtr_overlap_fraction=0.0,
            sex_homology_overlap_fraction=0.0,
            segmental_duplication_overlap_fraction=0.01,
            low_mappability_overlap_fraction=0.0,
            gap_centromere_telomere_overlap_fraction=0.0,
            repeat_rich_overlap_fraction=0.0,
            blacklist_overlap_fraction=0.0,
            ambiguous_alignment_overlap_fraction=0.0,
            high_risk_boundary_crossing=0,
            caller="segment_level_detector",
            segment_mean_robust_z=2.44,
            segment_median_robust_z=2.40,
            segment_abs_max_robust_z=6.70,
        )

        result = classify_event(
            row=row,
            chrom_bin_count=48,
            args=self.build_args(),
            sex_call="XY",
            par_regions={},
        )

        self.assertEqual(result["artifact_status"], "review")
        self.assertEqual(result["keep_event"], 1)
        self.assertIn("clean_event_below_ultra_pass", result["artifact_flags"])
        self.assertIn("low_risk_event_preserved_below_ultra_pass", result["downgrade_reason"])

    def test_focal_low_risk_segmental_dup_event_with_strong_signal_is_review(self):
        row = types.SimpleNamespace(
            chrom="chr15",
            start=24000000,
            end=29000000,
            start_bin=24,
            end_bin=29,
            n_bins=5,
            calibrated_mean_z=3.21,
            calibrated_median_z=3.89,
            event_corr_adjusted_z=7.01,
            empirical_qvalue=0.02,
            clean_bin_fraction=0.24561,
            moderate_risk_bin_fraction=0.75439,
            high_risk_bin_fraction=0.0,
            effective_bin_count=5.8,
            region_risk_score_mean=0.048253,
            region_risk_score_max=0.170908,
            xtr_overlap_fraction=0.0,
            sex_homology_overlap_fraction=0.0,
            segmental_duplication_overlap_fraction=0.212207,
            low_mappability_overlap_fraction=0.0,
            gap_centromere_telomere_overlap_fraction=0.0,
            repeat_rich_overlap_fraction=0.0,
            blacklist_overlap_fraction=0.0,
            ambiguous_alignment_overlap_fraction=0.0,
            high_risk_boundary_crossing=0,
            caller="segment_level_detector",
            segment_mean_robust_z=-2.48,
            segment_median_robust_z=-3.05,
            segment_abs_max_robust_z=3.50,
        )

        result = classify_event(
            row=row,
            chrom_bin_count=103,
            args=self.build_args(),
            sex_call="XX",
            par_regions={},
        )

        self.assertEqual(result["artifact_status"], "review")
        self.assertEqual(result["keep_event"], 1)
        self.assertIn("focal_low_risk_event_preserved_below_ultra_pass", result["downgrade_reason"])

    def test_segmental_duplication_event_with_limited_clean_support_is_artifact(self):
        row = types.SimpleNamespace(
            chrom="chr5",
            start=68000000,
            end=73000000,
            start_bin=68,
            end_bin=72,
            n_bins=5,
            calibrated_mean_z=-3.4,
            calibrated_median_z=-3.8,
            event_corr_adjusted_z=-5.1,
            empirical_qvalue=0.04,
            clean_bin_fraction=0.33,
            moderate_risk_bin_fraction=0.67,
            high_risk_bin_fraction=0.0,
            effective_bin_count=6.0,
            region_risk_score_mean=0.05,
            region_risk_score_max=0.18,
            xtr_overlap_fraction=0.0,
            sex_homology_overlap_fraction=0.0,
            segmental_duplication_overlap_fraction=0.39,
            low_mappability_overlap_fraction=0.0,
            gap_centromere_telomere_overlap_fraction=0.0,
            repeat_rich_overlap_fraction=0.0,
            blacklist_overlap_fraction=0.0,
            ambiguous_alignment_overlap_fraction=0.0,
            high_risk_boundary_crossing=0,
            caller="segment_level_detector",
            segment_mean_robust_z=-3.2,
            segment_median_robust_z=-3.5,
            segment_abs_max_robust_z=5.2,
        )

        result = classify_event(
            row=row,
            chrom_bin_count=160,
            args=self.build_args(),
            sex_call="XX",
            par_regions={},
        )

        self.assertEqual(result["artifact_status"], "artifact")
        self.assertEqual(result["keep_event"], 0)
        self.assertIn("segmental_duplication_overlap", result["artifact_flags"])
        self.assertIn("segmental_duplication_overlap_with_limited_clean_support", result["filter_reason"])


if __name__ == "__main__":
    unittest.main()
