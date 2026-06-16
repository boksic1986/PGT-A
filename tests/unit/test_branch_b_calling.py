from __future__ import annotations

import sys
import tempfile
import types
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

try:
    import numpy as np
    import pandas as pd
except ModuleNotFoundError as exc:  # pragma: no cover - environment-dependent skip
    np = None
    pd = None
    IMPORT_ERROR = exc
else:
    IMPORT_ERROR = None
    from pgta.predict.branch_b.common import load_sample_bins
    from pgta.predict.branch_b.calling import (
        annotate_calling_z_scores,
        build_a_refined_candidate_events,
        build_candidate_events,
        build_chromosome_dosage_event,
        build_segment_level_events,
        iter_seed_blocks,
        merge_adjacent_segment_calls,
        reconcile_segment_state,
        segment_events_for_hmm_role,
    )


@unittest.skipIf(IMPORT_ERROR is not None, f"optional dependency missing: {IMPORT_ERROR}")
class BranchBCallingSignalTest(unittest.TestCase):
    def test_sidecar_hmm_does_not_emit_standalone_candidate(self):
        hmm_segments = [
            {
                "chrom": "chr1",
                "start": 1,
                "end": 1000000,
                "state": "gain",
                "caller": "segment_level_detector",
            }
        ]

        self.assertEqual(segment_events_for_hmm_role(hmm_segments, "sidecar"), [])
        self.assertEqual(segment_events_for_hmm_role(hmm_segments, "disabled"), [])
        self.assertEqual(segment_events_for_hmm_role(hmm_segments, "legacy_candidate"), hmm_segments)

    def test_a_candidates_are_refined_without_independent_branch_b_discovery(self):
        bins_df = pd.DataFrame(
            {
                "chrom": ["chr14"] * 5,
                "start": [0, 100, 200, 300, 400],
                "end": [100, 200, 300, 400, 500],
                "bin_index": [0, 1, 2, 3, 4],
                "signal_for_calling": [10.0, 10.1, 10.0, 10.1, 10.0],
                "normalized_signal": [10.0, 10.1, 10.0, 10.1, 10.0],
                "robust_z": [0.1, 0.2, 0.1, 0.2, 0.1],
                "calling_weight": [1.0] * 5,
                "bin_weight": [1.0] * 5,
                "region_risk_weight": [1.0] * 5,
            }
        )
        a_candidates = pd.DataFrame(
            [
                {
                    "candidate_id": "Y2.A0001",
                    "sample_id": "Y2",
                    "chrom": "chr14",
                    "start": 100,
                    "end": 400,
                    "state": "gain",
                    "svtype": "DUP",
                    "a_zscore": 121.5,
                    "a_abs_zscore": 121.5,
                    "a_ratio": 0.09,
                    "a_rank": 1,
                    "a_support_level": "strong",
                    "a_source_event_count": 1,
                }
            ]
        )
        args = types.SimpleNamespace(sample_id="Y2", branch="B", correction_model="2d_loess_gc_mappability")

        events = build_a_refined_candidate_events(args, bins_df, a_candidates)

        self.assertEqual(len(events), 1)
        event = events[0]
        self.assertEqual(event["caller"], "wisecondorx_a_branch")
        self.assertEqual(event["caller_stage"], "a_branch_refinement")
        self.assertEqual(event["a_candidate_id"], "Y2.A0001")
        self.assertEqual(event["state"], "gain")
        self.assertEqual((int(event["start_bin"]), int(event["end_bin"])), (1, 3))
        self.assertEqual(int(event["n_bins"]), 3)

    def test_empty_explicit_a_candidates_suppress_branch_b_discovery(self):
        bins_df = pd.DataFrame(
            {
                "chrom": ["chr21"] * 6,
                "start": [index * 100 for index in range(6)],
                "end": [(index + 1) * 100 for index in range(6)],
                "bin_index": list(range(6)),
                "signal_for_calling": [20.0] * 6,
                "normalized_signal": [20.0] * 6,
                "robust_z": [5.0] * 6,
                "calling_weight": [1.0] * 6,
                "bin_weight": [1.0] * 6,
                "region_risk_weight": [1.0] * 6,
                "variance_inflation": [1.0] * 6,
                "calling_seed_eligible": [1] * 6,
                "is_autosome": [1] * 6,
            }
        )
        args = types.SimpleNamespace(
            sample_id="H16",
            branch="B",
            correction_model="2d_loess_gc_mappability",
            hmm_role="sidecar",
            min_bins=3,
            max_segments_per_chrom=12,
            split_threshold=2.5,
            hmm_state_shift=2.5,
            hmm_stay_prob=0.995,
            min_event_bins=3,
            min_event_z=1.5,
            chromosome_z_threshold=1.8,
            chromosome_min_effective_bins=3.0,
            chromosome_min_clean_fraction=0.0,
            raw_chromosome_z_threshold=1.8,
            masked_gap_rescue_min_abs_local_z=1.8,
            masked_gap_rescue_min_median_z=4.0,
            masked_gap_rescue_min_chrom_shift_z=0.75,
            fusion_iou_threshold=0.8,
            a_candidates_provided=True,
        )
        logger = types.SimpleNamespace(info=lambda *args, **kwargs: None)

        _, events = build_candidate_events(args, bins_df, logger, a_candidates_df=pd.DataFrame())

        self.assertEqual(events, [])

    def test_chromosome_shift_survives_chrom_local_normalization(self):
        neutral_values = [10.0, 10.1, 9.9, 10.0, 10.2, 9.8]
        bins_df = pd.DataFrame(
            {
                "chrom": (["chr1"] * 6 + ["chr2"] * 6 + ["chr3"] * 6 + ["chr4"] * 6 + ["chr5"] * 6 + ["chr21"] * 6),
                "start": [index * 10 for index in range(36)],
                "end": [(index + 1) * 10 for index in range(36)],
                "bin_index": [0, 1, 2, 3, 4, 5] * 6,
                "signal_for_calling": neutral_values * 5 + [7.0, 7.1, 6.9, 7.0, 7.1, 6.8],
                "normalized_signal": neutral_values * 5 + [7.0, 7.1, 6.9, 7.0, 7.1, 6.8],
                "bin_weight": [1.0] * 36,
                "region_risk_weight": [1.0] * 36,
                "variance_inflation": [1.0] * 36,
                "calling_weight": [1.0] * 36,
                "calling_seed_eligible": [1] * 36,
                "is_autosome": [1] * 36,
            }
        )

        scored_df, background = annotate_calling_z_scores(bins_df, "signal_for_calling")

        chr21 = scored_df[scored_df["chrom"] == "chr21"].copy()
        chr1 = scored_df[scored_df["chrom"] == "chr1"].copy()

        self.assertLess(background["chrom_shift_z"]["chr21"], -2.0)
        self.assertLess(abs(float(chr21["raw_local_robust_z"].mean())), 0.5)
        self.assertLess(float(chr21["robust_z"].mean()), -2.0)
        self.assertLess(abs(float(chr1["robust_z"].mean())), 1.0)

    def test_balanced_chromosome_keeps_local_structure_without_large_shift(self):
        bins_df = pd.DataFrame(
            {
                "chrom": ["chr7"] * 8 + ["chr1"] * 8,
                "start": [index * 10 for index in range(16)],
                "end": [(index + 1) * 10 for index in range(16)],
                "bin_index": list(range(8)) + list(range(8)),
                "signal_for_calling": [8.5, 8.4, 8.6, 8.5, 11.5, 11.4, 11.6, 11.5] + [10.0, 10.1, 9.9, 10.0, 10.1, 9.8, 10.0, 10.2],
                "normalized_signal": [8.5, 8.4, 8.6, 8.5, 11.5, 11.4, 11.6, 11.5] + [10.0, 10.1, 9.9, 10.0, 10.1, 9.8, 10.0, 10.2],
                "bin_weight": [1.0] * 16,
                "region_risk_weight": [1.0] * 16,
                "variance_inflation": [1.0] * 16,
                "calling_weight": [1.0] * 16,
                "calling_seed_eligible": [1] * 16,
                "is_autosome": [1] * 16,
            }
        )

        scored_df, background = annotate_calling_z_scores(bins_df, "signal_for_calling")
        chr7 = scored_df[scored_df["chrom"] == "chr7"].copy()

        self.assertLess(abs(background["chrom_shift_z"]["chr7"]), 1.0)
        self.assertLess(float(chr7.iloc[:4]["robust_z"].mean()), -0.5)
        self.assertGreater(float(chr7.iloc[4:]["robust_z"].mean()), 0.5)

    def test_chromosome_dosage_detector_emits_broad_event_for_global_shift(self):
        neutral_values = [10.0, 10.1, 9.9, 10.0, 10.2, 9.8]
        bins_df = pd.DataFrame(
            {
                "chrom": (["chr1"] * 6 + ["chr2"] * 6 + ["chr3"] * 6 + ["chr4"] * 6 + ["chr5"] * 6 + ["chr21"] * 6),
                "start": [index * 10 for index in range(36)],
                "end": [(index + 1) * 10 for index in range(36)],
                "bin_index": [0, 1, 2, 3, 4, 5] * 6,
                "signal_for_calling": neutral_values * 5 + [7.0, 7.1, 6.9, 7.0, 7.1, 6.8],
                "normalized_signal": neutral_values * 5 + [7.0, 7.1, 6.9, 7.0, 7.1, 6.8],
                "bin_weight": [1.0] * 36,
                "region_risk_weight": [1.0] * 36,
                "variance_inflation": [1.0] * 36,
                "calling_weight": [1.0] * 36,
                "calling_seed_eligible": [1] * 36,
                "is_autosome": [1] * 36,
                "region_risk_class": ["clean"] * 36,
            }
        )

        scored_df, _ = annotate_calling_z_scores(bins_df, "signal_for_calling")
        args = types.SimpleNamespace(
            sample_id="TEST001",
            branch="B",
            correction_model="2d_loess_gc_mappability",
            chromosome_z_threshold=2.5,
            chromosome_min_effective_bins=5.0,
            chromosome_min_clean_fraction=0.30,
        )

        chr21 = scored_df[scored_df["chrom"] == "chr21"].copy()
        events = build_chromosome_dosage_event(args, chr21)

        self.assertEqual(len(events), 1)
        event = events[0]
        self.assertEqual(event["caller"], "chromosome_dosage_detector")
        self.assertEqual(event["caller_stage"], "chromosome_dosage")
        self.assertEqual(event["state"], "loss")
        self.assertEqual(event["chrom"], "chr21")
        self.assertEqual(event["start_bin"], 0)
        self.assertEqual(event["end_bin"], 5)

    def test_raw_chromosome_dosage_detector_recovers_shift_flattened_by_correction(self):
        neutral_values = [10.0, 10.1, 9.9, 10.0, 10.2, 9.8]
        corrected_neutral = [10.0, 10.0, 10.0, 10.0, 10.1, 9.9]
        bins_df = pd.DataFrame(
            {
                "chrom": (["chr1"] * 6 + ["chr2"] * 6 + ["chr3"] * 6 + ["chr4"] * 6 + ["chr5"] * 6 + ["chr21"] * 6),
                "start": [index * 10 for index in range(36)],
                "end": [(index + 1) * 10 for index in range(36)],
                "bin_index": [0, 1, 2, 3, 4, 5] * 6,
                "signal_for_calling": corrected_neutral * 6,
                "normalized_signal": neutral_values * 5 + [13.0, 13.1, 12.9, 13.0, 13.2, 12.8],
                "bin_weight": [1.0] * 36,
                "region_risk_weight": [1.0] * 36,
                "variance_inflation": [1.0] * 36,
                "calling_weight": [1.0] * 36,
                "calling_seed_eligible": [1] * 36,
                "is_autosome": [1] * 36,
                "region_risk_class": ["clean"] * 30 + ["moderate", "moderate", "moderate", "moderate", "moderate", "high"],
                "mask_label": ["pass"] * 36,
                "gap_centromere_telomere_overlap_fraction": [0.0] * 36,
                "segment_id": [""] * 36,
                "hmm_state": ["neutral"] * 36,
            }
        )

        scored_df, background = annotate_calling_z_scores(bins_df, "signal_for_calling")
        args = types.SimpleNamespace(
            sample_id="TEST_RAW_BROAD",
            branch="B",
            correction_model="2d_loess_gc_mappability",
            min_bins=2,
            split_threshold=1.0,
            max_segments_per_chrom=8,
            hmm_state_shift=0.3,
            hmm_stay_prob=0.90,
            min_event_bins=2,
            min_event_z=0.6,
            chromosome_z_threshold=2.5,
            chromosome_min_effective_bins=5.0,
            chromosome_min_clean_fraction=0.30,
            masked_gap_rescue_min_abs_local_z=1.0,
            masked_gap_rescue_min_median_z=2.0,
            masked_gap_rescue_min_chrom_shift_z=0.75,
            fusion_iou_threshold=0.8,
        )
        logger = types.SimpleNamespace(info=lambda *args, **kwargs: None)

        chr21 = scored_df[scored_df["chrom"] == "chr21"].copy()
        _, events = build_candidate_events(args, chr21, logger)

        self.assertLess(abs(background["chrom_shift_z"]["chr21"]), 1.0)
        raw_broad_events = [event for event in events if event["caller"] == "raw_chromosome_dosage_detector"]
        self.assertEqual(len(raw_broad_events), 1)
        event = raw_broad_events[0]
        self.assertEqual(event["state"], "gain")
        self.assertEqual(event["chrom"], "chr21")
        self.assertEqual(event["start_bin"], 0)
        self.assertEqual(event["end_bin"], 5)

    def test_raw_chromosome_dosage_detector_blocks_low_clean_raw_loss(self):
        neutral_values = [10.0, 10.1, 9.9, 10.0, 10.2, 9.8]
        corrected_neutral = [10.0, 10.0, 10.0, 10.0, 10.1, 9.9]
        bins_df = pd.DataFrame(
            {
                "chrom": (["chr1"] * 6 + ["chr2"] * 6 + ["chr3"] * 6 + ["chr4"] * 6 + ["chr5"] * 6 + ["chr21"] * 6),
                "start": [index * 10 for index in range(36)],
                "end": [(index + 1) * 10 for index in range(36)],
                "bin_index": [0, 1, 2, 3, 4, 5] * 6,
                "signal_for_calling": corrected_neutral * 6,
                "normalized_signal": neutral_values * 5 + [6.0, 6.1, 5.9, 6.0, 6.2, 5.8],
                "bin_weight": [1.0] * 36,
                "region_risk_weight": [1.0] * 36,
                "variance_inflation": [1.0] * 36,
                "calling_weight": [1.0] * 36,
                "calling_seed_eligible": [1] * 36,
                "is_autosome": [1] * 36,
                "region_risk_class": ["clean"] * 30 + ["moderate", "moderate", "moderate", "moderate", "moderate", "high"],
                "mask_label": ["pass"] * 36,
                "gap_centromere_telomere_overlap_fraction": [0.0] * 36,
                "segment_id": [""] * 36,
                "hmm_state": ["neutral"] * 36,
            }
        )

        scored_df, background = annotate_calling_z_scores(bins_df, "signal_for_calling")
        args = types.SimpleNamespace(
            sample_id="TEST_RAW_LOW_CLEAN_LOSS",
            branch="B",
            correction_model="2d_loess_gc_mappability",
            min_bins=2,
            split_threshold=1.0,
            max_segments_per_chrom=8,
            hmm_state_shift=0.3,
            hmm_stay_prob=0.90,
            min_event_bins=2,
            min_event_z=0.6,
            chromosome_z_threshold=2.5,
            chromosome_min_effective_bins=5.0,
            chromosome_min_clean_fraction=0.30,
            masked_gap_rescue_min_abs_local_z=1.0,
            masked_gap_rescue_min_median_z=2.0,
            masked_gap_rescue_min_chrom_shift_z=0.75,
            fusion_iou_threshold=0.8,
        )
        logger = types.SimpleNamespace(info=lambda *args, **kwargs: None)

        chr21 = scored_df[scored_df["chrom"] == "chr21"].copy()
        _, events = build_candidate_events(args, chr21, logger)

        self.assertLess(abs(background["chrom_shift_z"]["chr21"]), 1.0)
        self.assertLess(background["raw_chrom_shift_z"]["chr21"], -2.0)
        self.assertFalse([event for event in events if event["caller"] == "raw_chromosome_dosage_detector"])

    def test_segment_level_detector_emits_focal_event_without_large_chromosome_shift(self):
        focal_signal = [7.0, 7.1, 7.0, 7.2, 12.8, 13.0, 12.9, 13.1, 7.0, 7.1, 7.0, 7.2]
        background_signal = [10.0, 10.1, 9.9, 10.0, 10.1, 9.8, 10.0, 10.2, 9.9, 10.0, 10.1, 9.8]
        bins_df = pd.DataFrame(
            {
                "chrom": ["chr7"] * len(focal_signal) + ["chr1"] * len(background_signal) + ["chr2"] * len(background_signal) + ["chr3"] * len(background_signal),
                "start": [index * 10 for index in range(len(focal_signal) + len(background_signal) * 3)],
                "end": [(index + 1) * 10 for index in range(len(focal_signal) + len(background_signal) * 3)],
                "bin_index": list(range(len(focal_signal))) + list(range(len(background_signal))) * 3,
                "signal_for_calling": focal_signal + background_signal + background_signal + background_signal,
                "normalized_signal": focal_signal + background_signal + background_signal + background_signal,
                "bin_weight": [1.0] * (len(focal_signal) + len(background_signal) * 3),
                "region_risk_weight": [1.0] * (len(focal_signal) + len(background_signal) * 3),
                "variance_inflation": [1.0] * (len(focal_signal) + len(background_signal) * 3),
                "calling_weight": [1.0] * (len(focal_signal) + len(background_signal) * 3),
                "calling_seed_eligible": [1] * (len(focal_signal) + len(background_signal) * 3),
                "is_autosome": [1] * (len(focal_signal) + len(background_signal) * 3),
                "region_risk_class": ["clean"] * (len(focal_signal) + len(background_signal) * 3),
            }
        )

        scored_df, background = annotate_calling_z_scores(bins_df, "signal_for_calling")
        args = types.SimpleNamespace(
            sample_id="TEST002",
            branch="B",
            correction_model="2d_loess_gc_mappability",
            min_bins=2,
            split_threshold=1.0,
            max_segments_per_chrom=8,
            hmm_state_shift=0.3,
            hmm_stay_prob=0.90,
            min_event_bins=2,
            min_event_z=0.6,
        )
        logger = types.SimpleNamespace(info=lambda *args, **kwargs: None)

        chr7 = scored_df[scored_df["chrom"] == "chr7"].copy()
        _, segment_events, _ = build_segment_level_events(args, chr7, logger)

        self.assertLess(abs(background["chrom_shift_z"]["chr7"]), 1.0)
        self.assertTrue(segment_events)
        focal_event = segment_events[0]
        self.assertEqual(focal_event["caller"], "segment_level_detector")
        self.assertEqual(focal_event["caller_stage"], "segment_level")
        self.assertEqual(focal_event["chrom"], "chr7")
        self.assertGreaterEqual(focal_event["n_bins"], 2)

    def test_segment_state_reconcile_prefers_branch_b_signal_direction(self):
        self.assertEqual(reconcile_segment_state("loss", 3.5, 3.2, 0.8), "gain")
        self.assertEqual(reconcile_segment_state("neutral", -2.2, -1.8, 0.8), "loss")
        self.assertEqual(reconcile_segment_state("gain", 0.2, 0.1, 0.8), "neutral")

    def test_load_sample_bins_accepts_numeric_sex_chromosome_keys(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            npz_path = Path(tmpdir) / "sample.npz"
            sample_payload = {
                "1": np.array([10.0, 12.0], dtype=np.float64),
                "23": np.array([5.0, 4.0], dtype=np.float64),
                "24": np.array([1.0], dtype=np.float64),
            }
            np.savez(npz_path, sample=sample_payload, binsize=np.array(1000000), quality=np.array(0.8))

            bins_df, binsize, quality = load_sample_bins(npz_path)

        self.assertEqual(binsize, 1000000)
        self.assertAlmostEqual(float(quality), 0.8, places=6)
        self.assertEqual(sorted(bins_df["chrom"].unique().tolist()), ["chr1", "chrX", "chrY"])
        self.assertEqual(int((bins_df["chrom"] == "chrX").sum()), 2)
        self.assertEqual(int((bins_df["chrom"] == "chrY").sum()), 1)

    def test_merge_adjacent_segment_calls_merges_same_state_within_seed_block(self):
        merged = merge_adjacent_segment_calls(
            [
                {
                    "chrom": "chr21",
                    "block_index": 1,
                    "segment_id": "chr21_block1_seg1",
                    "left": 10,
                    "right": 14,
                    "state": "gain",
                    "mean_z": 2.2,
                    "median_z": 2.0,
                    "source_segment_ids": ["chr21_block1_seg1"],
                },
                {
                    "chrom": "chr21",
                    "block_index": 1,
                    "segment_id": "chr21_block1_seg2",
                    "left": 14,
                    "right": 18,
                    "state": "gain",
                    "mean_z": 2.8,
                    "median_z": 2.6,
                    "source_segment_ids": ["chr21_block1_seg2"],
                },
                {
                    "chrom": "chr21",
                    "block_index": 2,
                    "segment_id": "chr21_block2_seg1",
                    "left": 20,
                    "right": 24,
                    "state": "gain",
                    "mean_z": 3.1,
                    "median_z": 3.0,
                    "source_segment_ids": ["chr21_block2_seg1"],
                },
            ]
        )

        self.assertEqual(len(merged), 2)
        self.assertEqual((merged[0]["left"], merged[0]["right"]), (10, 18))
        self.assertEqual(merged[0]["state"], "gain")
        self.assertEqual(merged[0]["source_segment_ids"], ["chr21_block1_seg1", "chr21_block1_seg2"])
        self.assertEqual((merged[1]["left"], merged[1]["right"]), (20, 24))

    def test_segment_detector_skips_hard_masked_seed_block(self):
        chrom_df = pd.DataFrame(
            {
                "chrom": ["chr13"] * 12,
                "start": [index * 10 for index in range(12)],
                "end": [(index + 1) * 10 for index in range(12)],
                "bin_index": list(range(12)),
                "normalized_signal": [5.0] * 4 + [10.0] * 4 + [13.0] * 4,
                "calling_weight": [1.0] * 12,
                "calling_seed_eligible": [0] * 4 + [1] * 8,
                "robust_z": [-5.0] * 4 + [0.0] * 4 + [3.0] * 4,
                "segment_id": [""] * 12,
                "hmm_state": ["neutral"] * 12,
            }
        )
        args = types.SimpleNamespace(
            sample_id="TEST003",
            branch="B",
            correction_model="2d_loess_gc_mappability",
            min_bins=2,
            split_threshold=1.0,
            max_segments_per_chrom=8,
            hmm_state_shift=0.3,
            hmm_stay_prob=0.90,
            min_event_bins=2,
            min_event_z=0.6,
        )
        logger = types.SimpleNamespace(info=lambda *args, **kwargs: None)

        self.assertEqual(iter_seed_blocks(chrom_df, min_block_bins=2), [(4, 12)])
        _, segment_events, _ = build_segment_level_events(args, chrom_df, logger)

        self.assertTrue(segment_events)
        self.assertTrue(all(int(event["start_bin"]) >= 4 for event in segment_events))
        self.assertTrue(any(event["state"] == "gain" and int(event["start_bin"]) >= 8 for event in segment_events))

    def test_masked_gap_rescue_emits_loss_against_positive_chrom_shift(self):
        args = types.SimpleNamespace(
            sample_id="TEST004",
            branch="B",
            correction_model="2d_loess_gc_mappability",
            min_bins=3,
            split_threshold=2.5,
            max_segments_per_chrom=12,
            hmm_state_shift=2.5,
            hmm_stay_prob=0.995,
            min_event_bins=3,
            min_event_z=1.5,
            chromosome_z_threshold=2.5,
            chromosome_min_effective_bins=10.0,
            chromosome_min_clean_fraction=0.30,
            masked_gap_rescue_min_abs_local_z=1.0,
            masked_gap_rescue_min_median_z=2.0,
            masked_gap_rescue_min_chrom_shift_z=0.75,
            fusion_iou_threshold=0.8,
        )
        logger = types.SimpleNamespace(info=lambda *args, **kwargs: None)
        chrom_df = pd.DataFrame(
            {
                "chrom": ["chr7"] * 12,
                "start": [index * 1000000 for index in range(12)],
                "end": [(index + 1) * 1000000 for index in range(12)],
                "bin_index": list(range(12)),
                "normalized_signal": [6.0, 6.1, 6.0, 6.2, 4.4, 0.0, 0.0, 7.5, 6.2, 6.1, 6.0, 6.1],
                "signal_for_calling": [6.0, 6.1, 6.0, 6.2, 8.8, 4.5, 4.5, 8.9, 6.2, 6.1, 6.0, 6.1],
                "bin_weight": [10.0] * 12,
                "region_risk_weight": [1.0, 1.0, 1.0, 1.0, 0.1, 0.1, 0.1, 0.2, 1.0, 1.0, 1.0, 1.0],
                "variance_inflation": [1.0, 1.0, 1.0, 1.0, 4.0, 4.0, 4.0, 3.25, 1.0, 1.0, 1.0, 1.0],
                "calling_weight": [10.0, 10.0, 10.0, 10.0, 1.0, 1.0, 1.0, 2.0, 10.0, 10.0, 10.0, 10.0],
                "calling_seed_eligible": [1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1],
                "is_autosome": [1] * 12,
                "region_risk_class": ["clean", "clean", "clean", "clean", "high", "high", "high", "high", "clean", "clean", "clean", "clean"],
                "mask_label": ["pass", "pass", "pass", "pass", "hard", "hard", "hard", "hard", "pass", "pass", "pass", "pass"],
                "gap_centromere_telomere_overlap_fraction": [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.7, 0.0, 0.0, 0.0, 0.0],
                "chrom_shift_z": [1.2] * 12,
                "robust_z": [0.2, 0.3, 0.1, 0.2, 2.4, -0.8, -0.8, 3.0, 0.4, 0.3, 0.2, 0.1],
                "segment_id": [""] * 12,
                "hmm_state": ["neutral"] * 12,
            }
        )

        _, events = build_candidate_events(args, chrom_df, logger)
        rescue_events = [event for event in events if event["caller"] == "masked_gap_rescue_detector"]

        self.assertEqual(len(rescue_events), 1)
        rescue = rescue_events[0]
        self.assertEqual(rescue["state"], "loss")
        self.assertEqual((int(rescue["start_bin"]), int(rescue["end_bin"])), (4, 6))
        self.assertGreater(abs(float(rescue["segment_median_robust_z"])), 4.0)


if __name__ == "__main__":
    unittest.main()
