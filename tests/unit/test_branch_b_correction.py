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
    from pgta.predict.branch_b.correction import (
        aggregate_mask_to_sample_bins,
        aggregate_reference_to_sample_bins,
        escalate_hotspot_masks,
    )


@unittest.skipIf(IMPORT_ERROR is not None, f"optional dependency missing: {IMPORT_ERROR}")
class BranchBCorrectionAnnotationProjectionTest(unittest.TestCase):
    def test_projection_preserves_gap_overlap_and_hard_mask_across_bin_size_mismatch(self):
        sample_bins = pd.DataFrame(
            {
                "chrom": ["chr13"],
                "start": [0],
                "end": [1_000_000],
            }
        )
        annotations = pd.DataFrame(
            {
                "chrom": ["chr13"] * 10,
                "start": [index * 100_000 for index in range(10)],
                "end": [(index + 1) * 100_000 for index in range(10)],
                "gc_fraction": [0.4] * 10,
                "atgc_fraction": [0.8] * 10,
                "n_fraction": [1.0] * 10,
                "effective_size": [100_000] * 10,
                "mappability_score": [0.8] * 10,
                "gap_centromere_telomere_overlap_fraction": [1.0] * 10,
                "segmental_duplication_overlap_fraction": [0.0] * 10,
                "par_overlap_fraction": [0.0] * 10,
                "xtr_overlap_fraction": [0.0] * 10,
                "sex_homology_overlap_fraction": [0.0] * 10,
                "low_mappability_overlap_fraction": [0.0] * 10,
                "repeat_rich_overlap_fraction": [0.0] * 10,
                "blacklist_overlap_fraction": [0.0] * 10,
                "ambiguous_alignment_overlap_fraction": [0.0] * 10,
                "is_autosome": [1] * 10,
                "is_gap_centromere_telomere": [1] * 10,
            }
        )
        masks = pd.DataFrame(
            {
                "chrom": ["chr13"] * 10,
                "start": [index * 100_000 for index in range(10)],
                "end": [(index + 1) * 100_000 for index in range(10)],
                "mask_label": ["hard"] * 10,
                "mask_reason": ["n_fraction>=0.2"] * 10,
            }
        )

        projected_annotations = aggregate_reference_to_sample_bins(
            sample_bins_df=sample_bins,
            source_df=annotations,
            value_columns=[
                "gc_fraction",
                "atgc_fraction",
                "n_fraction",
                "effective_size",
                "mappability_score",
                "gap_centromere_telomere_overlap_fraction",
                "segmental_duplication_overlap_fraction",
                "par_overlap_fraction",
                "xtr_overlap_fraction",
                "sex_homology_overlap_fraction",
                "low_mappability_overlap_fraction",
                "repeat_rich_overlap_fraction",
                "blacklist_overlap_fraction",
                "ambiguous_alignment_overlap_fraction",
            ],
            bool_columns=["is_autosome", "is_gap_centromere_telomere"],
        )
        projected_masks = aggregate_mask_to_sample_bins(sample_bins_df=sample_bins, mask_df=masks)
        projected = projected_annotations.merge(projected_masks, on=["chrom", "start", "end"], how="left")

        self.assertEqual(len(projected), 1)
        row = projected.iloc[0]
        self.assertEqual(row["mask_label"], "hard")
        self.assertAlmostEqual(float(row["gap_centromere_telomere_overlap_fraction"]), 1.0, places=6)
        self.assertEqual(int(row["is_gap_centromere_telomere"]), 1)

    def test_hotspot_mask_escalates_xy_homology(self):
        reference_df = pd.DataFrame(
            {
                "chrom": ["chrX", "chr16"],
                "start": [0, 0],
                "end": [1_000_000, 1_000_000],
                "mask_label": ["pass", "pass"],
                "mask_reason": ["", ""],
                "par_overlap_fraction": [0.0, 0.0],
                "xtr_overlap_fraction": [0.0, 0.0],
                "sex_homology_overlap_fraction": [0.65, 0.0],
                "ambiguous_alignment_overlap_fraction": [0.0, 0.0],
                "gap_centromere_telomere_overlap_fraction": [0.0, 0.0],
                "segmental_duplication_overlap_fraction": [0.0, 0.45],
            }
        )

        escalated = escalate_hotspot_masks(reference_df)

        chrx = escalated[escalated["chrom"] == "chrX"].iloc[0]
        chr16 = escalated[escalated["chrom"] == "chr16"].iloc[0]
        self.assertEqual(chrx["mask_label"], "hard")
        self.assertIn("branch_b_hotspot:xy_homology", chrx["mask_reason"])
        self.assertEqual(chr16["mask_label"], "dynamic")
        self.assertIn("branch_b_hotspot:repeat_hotspot", chr16["mask_reason"])


if __name__ == "__main__":
    unittest.main()
