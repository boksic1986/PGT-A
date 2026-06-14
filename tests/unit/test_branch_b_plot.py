import re
from pathlib import Path

import pandas as pd

from pgta.predict.branch_b.plot import build_cnv_plot_svg


def test_build_cnv_plot_svg_renders_final_gain_loss_only(tmp_path):
    bins = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1", "chr2", "chr2"],
            "bin_index": [0, 1, 2, 0, 1],
            "start": [0, 1_000_000, 2_000_000, 0, 1_000_000],
            "end": [1_000_000, 2_000_000, 3_000_000, 1_000_000, 2_000_000],
            "calibrated_z": [0.1, 4.2, 3.8, -3.4, -3.1],
            "mask_label": ["pass", "pass", "soft", "pass", "hard"],
            "region_risk_class": ["clean", "clean", "moderate", "clean", "high"],
        }
    )
    events = pd.DataFrame(
        {
            "event_id": ["evt1", "evt2", "evt3"],
            "sample_id": ["Y1", "Y1", "Y1"],
            "chrom": ["chr1", "chr2", "chr2"],
            "start": [1_000_000, 0, 0],
            "end": [3_000_000, 2_000_000, 1_000_000],
            "state": ["gain", "loss", "gain"],
            "artifact_status": ["review", "pass", "artifact"],
            "keep_event": [1, 1, 0],
            "priority_score": [12.0, 11.0, 1.5],
        }
    )
    a_branch = pd.DataFrame(
        {
            "chr": ["1"],
            "start": [1_000_000],
            "end": [3_000_000],
            "type": ["gain"],
            "zscore": [22.0],
        }
    )

    output_svg = tmp_path / "Y1.final_cnv.svg"
    build_cnv_plot_svg(
        sample_id="Y1",
        bins_df=bins,
        branch_b_events_df=events,
        a_branch_df=a_branch,
        output_svg=output_svg,
    )

    svg = output_svg.read_text(encoding="utf-8")
    assert "<svg" in svg
    assert "CNV final profile - Y1" in svg
    assert "Final gain" in svg
    assert "Final loss" in svg
    assert "Smoothed signal trend" in svg
    assert "Branch A" not in svg
    assert "Branch B" not in svg
    assert "rejected" not in svg
    assert "chr1" in svg
    assert '<circle cx="' in svg
    assert '<circle cx="' in svg and 'fill="#475569"' in svg
    assert re.search(r'<circle[^>]+fill="#dc2626"', svg) is None
    assert re.search(r'<circle[^>]+fill="#2563eb"', svg) is None
    assert "soft mask" not in svg
    assert "hard/high-risk mask" not in svg
