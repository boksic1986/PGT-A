from __future__ import annotations

import sys
import types
import unittest
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

try:
    import pgta.predict.branch_b.artifact_rules as artifact_rules
    from pgta.predict.branch_b.artifact_rules import classify_event
except ModuleNotFoundError as exc:  # pragma: no cover - environment-dependent skip
    artifact_rules = None
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
            broad_gain_rescue_min_abs_z=1.8,
            broad_gain_rescue_min_support_fraction=0.35,
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
            a_branch_review_min_abs_z=20.0,
            a_branch_discordant_protect_min_abs_z=50.0,
            branch_b_direction_min_abs_z=0.25,
            cnvseq_large_event_min_bp=10_000_000,
            cnvseq_boundary_max_abs_z=4.0,
            cnvseq_whole_chrom_available_fraction=0.90,
        )

    def test_cnvseq_context_reports_available_fraction_and_boundary_crossing(self):
        self.assertTrue(hasattr(artifact_rules, "summarize_cnvseq_event_context"))
        bins = pd.DataFrame(
            {
                "chrom": ["chr7"] * 8,
                "bin_index": list(range(8)),
                "mask_label": ["pass", "pass", "pass", "pass", "pass", "hard", "pass", "pass"],
                "blacklist_overlap_fraction": [0.0] * 8,
                "gap_centromere_telomere_overlap_fraction": [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                "ambiguous_alignment_overlap_fraction": [0.0] * 8,
                "region_risk_class": ["clean", "clean", "moderate", "high", "clean", "high", "clean", "clean"],
            }
        )
        row = types.SimpleNamespace(chrom="chr7", start_bin=1, end_bin=4)

        result = artifact_rules.summarize_cnvseq_event_context(row, bins)

        self.assertAlmostEqual(result["cnvseq_available_chrom_fraction"], 4 / 7)
        self.assertEqual(result["cnvseq_crosses_gap_or_centromere"], 1)
        self.assertGreater(result["cnvseq_gap_centromere_bin_fraction"], 0.0)

    def test_large_subchrom_boundary_event_with_weak_cnvseq_support_is_artifact(self):
        row = types.SimpleNamespace(
            chrom="chr7",
            start=58_000_000,
            end=104_000_000,
            start_bin=58,
            end_bin=103,
            n_bins=46,
            state="gain",
            calibrated_mean_z=3.1,
            calibrated_median_z=3.0,
            event_corr_adjusted_z=3.2,
            empirical_qvalue=0.04,
            clean_bin_fraction=0.88,
            moderate_risk_bin_fraction=0.10,
            high_risk_bin_fraction=0.02,
            effective_bin_count=42.0,
            region_risk_score_mean=0.04,
            region_risk_score_max=0.18,
            xtr_overlap_fraction=0.0,
            sex_homology_overlap_fraction=0.0,
            segmental_duplication_overlap_fraction=0.0,
            low_mappability_overlap_fraction=0.0,
            gap_centromere_telomere_overlap_fraction=0.02,
            repeat_rich_overlap_fraction=0.0,
            blacklist_overlap_fraction=0.0,
            ambiguous_alignment_overlap_fraction=0.0,
            high_risk_boundary_crossing=1,
            cnvseq_available_chrom_fraction=0.31,
            cnvseq_crosses_gap_or_centromere=1,
            cnvseq_gap_centromere_bin_fraction=0.05,
            caller="wisecondorx_a_branch",
            a_candidate_id="Y8.A0003",
            a_abs_zscore=8.5,
            segment_mean_robust_z=3.0,
            segment_median_robust_z=3.0,
            segment_abs_max_robust_z=3.6,
        )

        result = classify_event(
            row=row,
            chrom_bin_count=159,
            args=self.build_args(),
            sex_call="XY",
            par_regions={},
        )

        self.assertEqual(result["artifact_status"], "artifact")
        self.assertEqual(result["keep_event"], 0)
        self.assertIn("cnvseq_subchrom_boundary_weak_support", result["filter_reason"])

    def test_cnvseq_boundary_rule_uses_segment_level_support_not_single_bin_spike(self):
        row = types.SimpleNamespace(
            chrom="chr17",
            start=25_000_001,
            end=81_000_000,
            start_bin=25,
            end_bin=80,
            n_bins=56,
            state="loss",
            calibrated_mean_z=-0.7,
            calibrated_median_z=-0.6,
            event_corr_adjusted_z=-3.2,
            empirical_qvalue=0.05,
            clean_bin_fraction=0.80,
            moderate_risk_bin_fraction=0.15,
            high_risk_bin_fraction=0.05,
            effective_bin_count=49.0,
            region_risk_score_mean=0.04,
            region_risk_score_max=0.18,
            xtr_overlap_fraction=0.0,
            sex_homology_overlap_fraction=0.0,
            segmental_duplication_overlap_fraction=0.0,
            low_mappability_overlap_fraction=0.0,
            gap_centromere_telomere_overlap_fraction=0.0,
            repeat_rich_overlap_fraction=0.0,
            blacklist_overlap_fraction=0.0,
            ambiguous_alignment_overlap_fraction=0.0,
            high_risk_boundary_crossing=1,
            cnvseq_available_chrom_fraction=0.70,
            cnvseq_crosses_gap_or_centromere=0,
            cnvseq_gap_centromere_bin_fraction=0.0,
            cnvseq_region_class_transition=1,
            caller="wisecondorx_a_branch",
            a_candidate_id="Y2.A0004",
            a_abs_zscore=16.4,
            segment_mean_robust_z=-0.5,
            segment_median_robust_z=-0.4,
            segment_abs_max_robust_z=21.6,
        )

        result = classify_event(
            row=row,
            chrom_bin_count=78,
            args=self.build_args(),
            sex_call="XX",
            par_regions={},
        )

        self.assertEqual(result["artifact_status"], "artifact")
        self.assertEqual(result["keep_event"], 0)
        self.assertIn("cnvseq_subchrom_boundary_weak_support", result["filter_reason"])

    def test_large_boundary_event_with_strong_a_branch_support_remains_review(self):
        row = types.SimpleNamespace(
            chrom="chr21",
            start=14_000_000,
            end=42_000_000,
            start_bin=14,
            end_bin=41,
            n_bins=28,
            state="gain",
            calibrated_mean_z=1.5,
            calibrated_median_z=1.4,
            event_corr_adjusted_z=1.3,
            empirical_qvalue=0.60,
            clean_bin_fraction=0.40,
            moderate_risk_bin_fraction=0.55,
            high_risk_bin_fraction=0.05,
            effective_bin_count=22.0,
            region_risk_score_mean=0.05,
            region_risk_score_max=0.20,
            xtr_overlap_fraction=0.0,
            sex_homology_overlap_fraction=0.0,
            segmental_duplication_overlap_fraction=0.0,
            low_mappability_overlap_fraction=0.0,
            gap_centromere_telomere_overlap_fraction=0.04,
            repeat_rich_overlap_fraction=0.0,
            blacklist_overlap_fraction=0.0,
            ambiguous_alignment_overlap_fraction=0.0,
            high_risk_boundary_crossing=1,
            cnvseq_available_chrom_fraction=0.58,
            cnvseq_crosses_gap_or_centromere=1,
            cnvseq_gap_centromere_bin_fraction=0.08,
            caller="wisecondorx_a_branch",
            a_candidate_id="Y2.A0002",
            a_abs_zscore=83.9,
            segment_mean_robust_z=1.2,
            segment_median_robust_z=1.1,
            segment_abs_max_robust_z=1.6,
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
        self.assertIn("a_branch_strong_evidence_preserved_for_review", result["downgrade_reason"])

    def test_a_branch_very_strong_discordant_direction_remains_review_for_sensitivity(self):
        row = types.SimpleNamespace(
            chrom="chr21",
            start=14_000_000,
            end=42_000_000,
            start_bin=14,
            end_bin=41,
            n_bins=28,
            state="gain",
            calibrated_mean_z=-0.72,
            calibrated_median_z=-0.72,
            event_corr_adjusted_z=-2.85,
            empirical_qvalue=0.04,
            clean_bin_fraction=0.82,
            moderate_risk_bin_fraction=0.0,
            high_risk_bin_fraction=0.18,
            effective_bin_count=24.0,
            region_risk_score_mean=0.01,
            region_risk_score_max=0.06,
            xtr_overlap_fraction=0.0,
            sex_homology_overlap_fraction=0.0,
            segmental_duplication_overlap_fraction=0.0,
            low_mappability_overlap_fraction=0.0,
            gap_centromere_telomere_overlap_fraction=0.0,
            repeat_rich_overlap_fraction=0.0,
            blacklist_overlap_fraction=0.0,
            ambiguous_alignment_overlap_fraction=0.0,
            high_risk_boundary_crossing=1,
            cnvseq_available_chrom_fraction=0.79,
            cnvseq_crosses_gap_or_centromere=0,
            cnvseq_gap_centromere_bin_fraction=0.0,
            cnvseq_region_class_transition=1,
            caller="wisecondorx_a_branch",
            a_candidate_id="Y2.A0002",
            a_abs_zscore=83.9,
            segment_mean_robust_z=-0.72,
            segment_median_robust_z=-0.72,
            segment_abs_max_robust_z=1.64,
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
        self.assertNotIn("branch_b_direction_discordant_with_candidate_state", result["filter_reason"])

    def test_a_branch_strong_gain_with_branch_b_loss_direction_is_artifact(self):
        row = types.SimpleNamespace(
            chrom="chr21",
            start=15_000_001,
            end=42_000_000,
            start_bin=15,
            end_bin=41,
            n_bins=27,
            state="gain",
            calibrated_mean_z=-0.72,
            calibrated_median_z=-0.72,
            event_corr_adjusted_z=-2.85,
            empirical_qvalue=0.04,
            clean_bin_fraction=0.82,
            moderate_risk_bin_fraction=0.0,
            high_risk_bin_fraction=0.18,
            effective_bin_count=24.0,
            region_risk_score_mean=0.01,
            region_risk_score_max=0.06,
            xtr_overlap_fraction=0.0,
            sex_homology_overlap_fraction=0.0,
            segmental_duplication_overlap_fraction=0.0,
            low_mappability_overlap_fraction=0.0,
            gap_centromere_telomere_overlap_fraction=0.0,
            repeat_rich_overlap_fraction=0.0,
            blacklist_overlap_fraction=0.0,
            ambiguous_alignment_overlap_fraction=0.0,
            high_risk_boundary_crossing=1,
            cnvseq_available_chrom_fraction=0.79,
            cnvseq_crosses_gap_or_centromere=0,
            cnvseq_gap_centromere_bin_fraction=0.0,
            cnvseq_region_class_transition=1,
            caller="wisecondorx_a_branch",
            a_candidate_id="Y8.A0003",
            a_abs_zscore=24.17,
            segment_mean_robust_z=-0.72,
            segment_median_robust_z=-0.72,
            segment_abs_max_robust_z=1.64,
        )

        result = classify_event(
            row=row,
            chrom_bin_count=34,
            args=self.build_args(),
            sex_call="XY",
            par_regions={},
        )

        self.assertEqual(result["artifact_status"], "artifact")
        self.assertEqual(result["keep_event"], 0)
        self.assertIn("branch_b_direction_discordant_with_candidate_state", result["filter_reason"])

    def test_sex_chromosome_direction_discordance_is_artifact(self):
        row = types.SimpleNamespace(
            chrom="chrY",
            start=13_000_001,
            end=20_000_000,
            start_bin=13,
            end_bin=19,
            n_bins=7,
            state="gain",
            calibrated_mean_z=-2.16,
            calibrated_median_z=-2.24,
            event_corr_adjusted_z=-4.77,
            empirical_qvalue=0.001,
            clean_bin_fraction=0.89,
            moderate_risk_bin_fraction=0.0,
            high_risk_bin_fraction=0.11,
            effective_bin_count=6.9,
            region_risk_score_mean=0.01,
            region_risk_score_max=0.13,
            xtr_overlap_fraction=0.0,
            sex_homology_overlap_fraction=0.0,
            segmental_duplication_overlap_fraction=0.0,
            low_mappability_overlap_fraction=0.0,
            gap_centromere_telomere_overlap_fraction=0.0,
            repeat_rich_overlap_fraction=0.0,
            blacklist_overlap_fraction=0.0,
            ambiguous_alignment_overlap_fraction=0.0,
            high_risk_boundary_crossing=1,
            cnvseq_available_chrom_fraction=0.24,
            cnvseq_crosses_gap_or_centromere=0,
            cnvseq_gap_centromere_bin_fraction=0.0,
            cnvseq_region_class_transition=1,
            caller="wisecondorx_a_branch",
            a_candidate_id="Y1.A0005",
            a_abs_zscore=9.18,
            segment_mean_robust_z=-2.21,
            segment_median_robust_z=-2.30,
            segment_abs_max_robust_z=2.34,
        )

        result = classify_event(
            row=row,
            chrom_bin_count=25,
            args=self.build_args(),
            sex_call="XY",
            par_regions={},
        )

        self.assertEqual(result["artifact_status"], "artifact")
        self.assertEqual(result["keep_event"], 0)
        self.assertIn("branch_b_direction_discordant_with_candidate_state", result["filter_reason"])

    def test_true_chrX_loss_with_matching_branch_b_direction_remains_review(self):
        row = types.SimpleNamespace(
            chrom="chrX",
            start=61_000_001,
            end=155_000_000,
            start_bin=61,
            end_bin=154,
            n_bins=94,
            state="loss",
            calibrated_mean_z=-0.70,
            calibrated_median_z=-0.44,
            event_corr_adjusted_z=-5.66,
            empirical_qvalue=0.01,
            clean_bin_fraction=0.78,
            moderate_risk_bin_fraction=0.0,
            high_risk_bin_fraction=0.22,
            effective_bin_count=84.0,
            region_risk_score_mean=0.01,
            region_risk_score_max=0.05,
            xtr_overlap_fraction=0.0,
            sex_homology_overlap_fraction=0.0,
            segmental_duplication_overlap_fraction=0.0,
            low_mappability_overlap_fraction=0.0,
            gap_centromere_telomere_overlap_fraction=0.0,
            repeat_rich_overlap_fraction=0.0,
            blacklist_overlap_fraction=0.0,
            ambiguous_alignment_overlap_fraction=0.0,
            high_risk_boundary_crossing=1,
            cnvseq_available_chrom_fraction=0.62,
            cnvseq_crosses_gap_or_centromere=0,
            cnvseq_gap_centromere_bin_fraction=0.0,
            cnvseq_region_class_transition=1,
            caller="wisecondorx_a_branch",
            a_candidate_id="Y3.A0001",
            a_abs_zscore=73.74,
            segment_mean_robust_z=-0.61,
            segment_median_robust_z=-0.35,
            segment_abs_max_robust_z=26.66,
        )

        result = classify_event(
            row=row,
            chrom_bin_count=149,
            args=self.build_args(),
            sex_call="XX",
            par_regions={"chrX": [(60001, 2699520), (154931044, 155260560)]},
        )

        self.assertEqual(result["artifact_status"], "review")
        self.assertEqual(result["keep_event"], 1)
        self.assertNotIn("branch_b_direction_discordant_with_candidate_state", result["filter_reason"])

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

    def test_low_clean_broad_gain_detector_event_is_preserved_for_review(self):
        row = types.SimpleNamespace(
            chrom="chr21",
            start=1,
            end=48000000,
            start_bin=0,
            end_bin=48,
            n_bins=49,
            state="gain",
            calibrated_mean_z=2.2,
            calibrated_median_z=2.2,
            event_corr_adjusted_z=2.1,
            empirical_qvalue=0.04,
            clean_bin_fraction=0.18,
            moderate_risk_bin_fraction=0.45,
            high_risk_bin_fraction=0.37,
            effective_bin_count=33.0,
            region_risk_score_mean=0.04,
            region_risk_score_max=0.12,
            xtr_overlap_fraction=0.0,
            sex_homology_overlap_fraction=0.0,
            segmental_duplication_overlap_fraction=0.0,
            low_mappability_overlap_fraction=0.0,
            gap_centromere_telomere_overlap_fraction=0.05,
            repeat_rich_overlap_fraction=0.0,
            blacklist_overlap_fraction=0.0,
            ambiguous_alignment_overlap_fraction=0.0,
            high_risk_boundary_crossing=0,
            caller="raw_chromosome_dosage_detector",
            segment_mean_robust_z=2.17,
            segment_median_robust_z=2.17,
            segment_abs_max_robust_z=2.17,
        )

        result = classify_event(
            row=row,
            chrom_bin_count=49,
            args=self.build_args(),
            sex_call="XY",
            par_regions={},
        )

        self.assertEqual(result["artifact_status"], "review")
        self.assertEqual(result["keep_event"], 1)
        self.assertIn("broad_gain_preserved_by_raw_or_chromosome_support", result["downgrade_reason"])

    def test_low_level_raw_whole_chromosome_gain_gets_report_priority_floor(self):
        row = types.SimpleNamespace(
            chrom="chr14",
            start=0,
            end=108000000,
            start_bin=0,
            end_bin=107,
            n_bins=108,
            state="gain",
            calibrated_mean_z=1.863,
            calibrated_median_z=1.863,
            event_corr_adjusted_z=1.863,
            empirical_qvalue=0.062,
            clean_bin_fraction=0.220,
            moderate_risk_bin_fraction=0.774,
            high_risk_bin_fraction=0.006,
            effective_bin_count=86.3,
            region_risk_score_mean=0.04,
            region_risk_score_max=0.12,
            xtr_overlap_fraction=0.0,
            sex_homology_overlap_fraction=0.0,
            segmental_duplication_overlap_fraction=0.0,
            low_mappability_overlap_fraction=0.0,
            gap_centromere_telomere_overlap_fraction=0.05,
            repeat_rich_overlap_fraction=0.0,
            blacklist_overlap_fraction=0.0,
            ambiguous_alignment_overlap_fraction=0.0,
            high_risk_boundary_crossing=0,
            caller="raw_chromosome_dosage_detector",
            segment_mean_robust_z=1.863,
            segment_median_robust_z=1.863,
            segment_abs_max_robust_z=1.863,
        )

        result = classify_event(
            row=row,
            chrom_bin_count=108,
            args=self.build_args(),
            sex_call="XY",
            par_regions={},
        )

        self.assertEqual(result["artifact_status"], "review")
        self.assertEqual(result["keep_event"], 1)
        self.assertGreater(result["priority_score"], 7.0)

    def test_low_clean_broad_raw_loss_stays_artifact(self):
        row = types.SimpleNamespace(
            chrom="chr22",
            start=1,
            end=51000000,
            start_bin=0,
            end_bin=51,
            n_bins=52,
            state="loss",
            calibrated_mean_z=-2.6,
            calibrated_median_z=-2.6,
            event_corr_adjusted_z=-2.5,
            empirical_qvalue=0.02,
            clean_bin_fraction=0.12,
            moderate_risk_bin_fraction=0.48,
            high_risk_bin_fraction=0.40,
            effective_bin_count=31.0,
            region_risk_score_mean=0.05,
            region_risk_score_max=0.14,
            xtr_overlap_fraction=0.0,
            sex_homology_overlap_fraction=0.0,
            segmental_duplication_overlap_fraction=0.0,
            low_mappability_overlap_fraction=0.0,
            gap_centromere_telomere_overlap_fraction=0.05,
            repeat_rich_overlap_fraction=0.0,
            blacklist_overlap_fraction=0.0,
            ambiguous_alignment_overlap_fraction=0.0,
            high_risk_boundary_crossing=0,
            caller="raw_chromosome_dosage_detector",
            segment_mean_robust_z=-2.6,
            segment_median_robust_z=-2.6,
            segment_abs_max_robust_z=2.6,
        )

        result = classify_event(
            row=row,
            chrom_bin_count=52,
            args=self.build_args(),
            sex_call="XY",
            par_regions={},
        )

        self.assertEqual(result["artifact_status"], "artifact")
        self.assertEqual(result["keep_event"], 0)
        self.assertIn("chromosome_fraction_too_large", result["filter_reason"])

    def test_strong_a_branch_event_is_preserved_for_review_when_b_signal_is_weak(self):
        row = types.SimpleNamespace(
            chrom="chr21",
            start=15000001,
            end=42000000,
            start_bin=15,
            end_bin=42,
            n_bins=28,
            state="gain",
            calibrated_mean_z=1.2,
            calibrated_median_z=1.1,
            event_corr_adjusted_z=1.0,
            empirical_qvalue=0.70,
            clean_bin_fraction=0.22,
            moderate_risk_bin_fraction=0.74,
            high_risk_bin_fraction=0.04,
            effective_bin_count=24.0,
            region_risk_score_mean=0.04,
            region_risk_score_max=0.12,
            xtr_overlap_fraction=0.0,
            sex_homology_overlap_fraction=0.0,
            segmental_duplication_overlap_fraction=0.0,
            low_mappability_overlap_fraction=0.0,
            gap_centromere_telomere_overlap_fraction=0.0,
            repeat_rich_overlap_fraction=0.0,
            blacklist_overlap_fraction=0.0,
            ambiguous_alignment_overlap_fraction=0.0,
            high_risk_boundary_crossing=0,
            caller="wisecondorx_a_branch",
            a_candidate_id="Y2.A0002",
            a_abs_zscore=83.9,
            segment_mean_robust_z=1.1,
            segment_median_robust_z=1.0,
            segment_abs_max_robust_z=1.3,
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
        self.assertIn("a_branch_strong_evidence_preserved_for_review", result["downgrade_reason"])

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
