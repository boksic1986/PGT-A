from __future__ import annotations

import numpy as np
import pandas as pd

from pgta.predict.branch_b.ref_stability import (
    compute_ref_bin_stability,
    summarize_event_ref_stability,
)


def _write_npz(path, chrom_counts):
    sample = {chrom: np.asarray(counts, dtype=float) for chrom, counts in chrom_counts.items()}
    np.savez(path, sample=sample, binsize=np.array(1_000_000), quality=np.array(0.9))


def test_ref_bin_stability_computes_per_bin_mad(tmp_path):
    r1 = tmp_path / "R1.npz"
    r2 = tmp_path / "R2.npz"
    r3 = tmp_path / "R3.npz"
    _write_npz(r1, {"1": [100, 100, 100], "2": [100, 100, 100]})
    _write_npz(r2, {"1": [102, 500, 101], "2": [101, 100, 101]})
    _write_npz(r3, {"1": [98, 20, 99], "2": [99, 100, 99]})

    bins, summary = compute_ref_bin_stability(
        [r1, r2, r3],
        sample_ids=["R1", "R2", "R3"],
        high_mad_z=1.0,
    )

    stable = bins.loc[(bins["chrom"] == "chr1") & (bins["bin_index"] == 0)].iloc[0]
    variable = bins.loc[(bins["chrom"] == "chr1") & (bins["bin_index"] == 1)].iloc[0]

    assert int(variable["ref_sample_count"]) == 3
    assert variable["ref_mad"] > stable["ref_mad"]
    assert variable["ref_stability_label"] == "high_ref_mad"
    assert summary["ref_sample_count"] == 3
    assert summary["binsize"] == 1_000_000


def test_event_ref_stability_summarizes_high_mad_fraction():
    ref_bins = pd.DataFrame(
        [
            {
                "chrom": "chr1",
                "start": 0,
                "end": 1_000_000,
                "bin_index": 0,
                "ref_mad": 0.02,
                "ref_mad_z": 0.5,
                "ref_stability_label": "stable",
            },
            {
                "chrom": "chr1",
                "start": 1_000_000,
                "end": 2_000_000,
                "bin_index": 1,
                "ref_mad": 0.50,
                "ref_mad_z": 5.0,
                "ref_stability_label": "high_ref_mad",
            },
        ]
    )
    events = pd.DataFrame(
        [
            {
                "sample_id": "S1",
                "candidate_id": "S1_chr1_gain",
                "chrom": "chr1",
                "start": 0,
                "end": 2_000_000,
                "state": "gain",
            }
        ]
    )

    annotated = summarize_event_ref_stability(events, ref_bins, high_mad_z=4.0)
    row = annotated.iloc[0]

    assert row["event_ref_mad_median"] == 0.26
    assert row["high_ref_mad_bin_fraction"] == 0.5
    assert row["ref_stability_context"] == "REF_STABILITY_HIGH_MAD_REVIEW"
