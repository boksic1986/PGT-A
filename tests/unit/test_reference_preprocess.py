from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import numpy as np
import pandas as pd

from pgta.reference.assets import classify_reference_mask, write_hard_mask_bed
from pgta.reference.preprocess import build_gc_rc_fit_mask, mask_npz_sample


def test_chrM_is_hard_excluded_from_reference_mask():
    label, reason = classify_reference_mask({"chrom": "chrM", "start": 1, "end": 1000})

    assert label == "hard"
    assert "chrM" in reason


def test_par_is_annotation_not_global_hard_mask():
    label, _reason = classify_reference_mask(
        {
            "chrom": "chrX",
            "start": 60001,
            "end": 2699520,
            "par_overlap_fraction": 1.0,
            "blacklist_overlap_fraction": 0.0,
            "gap_centromere_telomere_overlap_fraction": 0.0,
        }
    )

    assert label != "hard"


def test_gap_blacklist_is_hard_mask():
    label, reason = classify_reference_mask(
        {
            "chrom": "chr1",
            "start": 1,
            "end": 1000,
            "gap_centromere_telomere_overlap_fraction": 1.0,
        }
    )

    assert label == "hard"
    assert "gap" in reason


def test_low_wisecondorx_ref_bins_are_hard_masked():
    label, reason = classify_reference_mask(
        {
            "chrom": "chr16",
            "start": 1_000_000,
            "end": 2_000_000,
            "ref_size_after_cutoff": 12,
            "segmental_duplication_overlap_fraction": 0.0,
            "repeat_rich_overlap_fraction": 0.0,
            "low_mappability_overlap_fraction": 0.0,
        }
    )

    assert label == "hard"
    assert "low_ref_bins" in reason


def test_intermediate_wisecondorx_ref_bins_are_dynamic_masked():
    label, reason = classify_reference_mask(
        {
            "chrom": "chr9",
            "start": 1_000_000,
            "end": 2_000_000,
            "ref_size_after_cutoff": 120,
            "segmental_duplication_overlap_fraction": 0.0,
            "repeat_rich_overlap_fraction": 0.0,
            "low_mappability_overlap_fraction": 0.0,
        }
    )

    assert label == "dynamic"
    assert "low_ref_bins" in reason


def test_low_ref_bins_in_repeat_region_are_hard_masked():
    label, reason = classify_reference_mask(
        {
            "chrom": "chr16",
            "start": 1_000_000,
            "end": 2_000_000,
            "ref_size_after_cutoff": 90,
            "segmental_duplication_overlap_fraction": 0.40,
        }
    )

    assert label == "hard"
    assert "repeat_or_lowmap" in reason


def test_proximal_low_ref_bins_are_hard_masked():
    label, reason = classify_reference_mask(
        {
            "chrom": "chr21",
            "start": 15_000_000,
            "end": 16_000_000,
            "is_near_centromere": 1,
            "ref_size_after_cutoff": 120,
            "segmental_duplication_overlap_fraction": 0.0,
            "repeat_rich_overlap_fraction": 0.0,
            "low_mappability_overlap_fraction": 0.0,
        }
    )

    assert label == "hard"
    assert "proximal_low_ref_bins" in reason


def test_proximal_high_dynamic_noise_is_hard_masked():
    label, reason = classify_reference_mask(
        {
            "chrom": "chr16",
            "start": 46_000_000,
            "end": 47_000_000,
            "nearest_telomere_distance_bp": 2_000_000,
            "dynamic_median_abs_z": 2.0,
            "dynamic_z_frac": 0.0,
        }
    )

    assert label == "hard"
    assert "proximal_high_dynamic_noise" in reason


def test_proximal_bin_without_ref_or_noise_risk_is_not_hard_masked():
    label, _reason = classify_reference_mask(
        {
            "chrom": "chr12",
            "start": 123_000_000,
            "end": 124_000_000,
            "is_near_telomere": 1,
            "ref_size_after_cutoff": 220,
            "dynamic_median_abs_z": 0.2,
            "dynamic_z_frac": 0.0,
        }
    )

    assert label != "hard"


def test_gc_rc_fit_uses_only_clean_autosomal_bins():
    bins = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chrX", "chr2"],
            "normalized_signal": [10.0, 11.0, 30.0, 12.0],
            "gc_fraction": [0.40, 0.45, 0.50, 0.60],
            "mappability_score": [1.0, 1.0, 1.0, 1.0],
            "mask_label": ["pass", "soft", "pass", "pass"],
            "is_autosome": [1, 1, 0, 1],
        }
    )

    fit_mask = build_gc_rc_fit_mask(bins, include_mask_labels={"pass"})

    assert fit_mask.tolist() == [True, False, False, True]


def test_mask_only_npz_zeroes_hard_bins_and_keeps_count_like_vectors():
    sample = {"1": np.asarray([10.0, 20.0, 30.0]), "2": np.asarray([5.0, 6.0, 7.0])}
    combined_mask = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr2"],
            "start": [100, 200, 100],
            "end": [200, 300, 200],
            "bin_size": [100, 100, 100],
            "mask_label": ["hard", "soft", "pass"],
        }
    )

    masked, metadata = mask_npz_sample(sample, 100, combined_mask)

    assert masked["1"].tolist() == [10.0, 0.0, 30.0]
    assert masked["2"].tolist() == [5.0, 6.0, 7.0]
    assert metadata["preprocess_strategy"] == "mask_only"
    assert metadata["mask_policy"] == "hard_mask_to_zero_soft_keep"
    assert metadata["count_like_contract_validated"] is True


def test_hard_mask_bed_exports_atomic_bed_without_header(tmp_path):
    hard_mask = pd.DataFrame(
        {
            "bin_level": ["atomic", "analysis"],
            "chrom": ["chr1", "chr1"],
            "start": [100, 0],
            "end": [200, 1000],
            "mask_reason": ["blacklist", "coarse"],
        }
    )
    output = tmp_path / "hard_mask.bed"

    write_hard_mask_bed(output, hard_mask)

    assert output.read_text(encoding="utf-8").splitlines() == ["chr1\t100\t200"]
