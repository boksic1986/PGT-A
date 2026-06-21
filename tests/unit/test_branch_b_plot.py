import re
import re
from pathlib import Path

import pandas as pd
import pytest

from pgta.predict.branch_b.plot import build_cnv_plot_svg


def test_build_cnv_plot_svg_uses_calibrated_z_and_writes_plot_bins(tmp_path):
    bins = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1", "chr2", "chr2", "chr2"],
            "bin_index": [0, 1, 2, 0, 1, 2],
            "start": [0, 1_000_000, 2_000_000, 0, 1_000_000, 2_000_000],
            "end": [1_000_000, 2_000_000, 3_000_000, 1_000_000, 2_000_000, 3_000_000],
            "calibrated_z": [0.1, 4.2, 3.8, -3.4, -3.1, pd.NA],
            "robust_z": [9.9, 9.9, 9.9, 9.9, 9.9, 9.9],
            "mask_label": ["pass", "pass", "soft", "pass", "hard", "pass"],
            "region_risk_class": ["clean", "clean", "moderate", "clean", "high", "clean"],
        }
    )
    events = pd.DataFrame(
        {
            "event_id": ["evt1", "evt2", "evt3"],
            "sample_id": ["Y1", "Y1", "Y2"],
            "chrom": ["chr1", "chr2", "chr1"],
            "start": [1_000_000, 0, 0],
            "end": [3_000_000, 2_000_000, 1_000_000],
            "state": ["gain", "loss", "gain"],
            "v2_report_layer_class": ["report_event", "report_event", "report_event"],
            "v2_report_visibility": ["report_strong_event", "report_weak_event", "report_strong_event"],
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
    output_bins = tmp_path / "Y1.plot_bins.tsv"
    build_cnv_plot_svg(
        sample_id="Y1",
        bins_df=bins,
        branch_b_events_df=events,
        a_branch_df=a_branch,
        output_svg=output_svg,
        output_bins_tsv=output_bins,
    )

    svg = output_svg.read_text(encoding="utf-8")
    plot_bins = pd.read_csv(output_bins, sep="\t")

    assert "<svg" in svg
    assert "CNV final profile - Y1" in svg
    assert set(["chrom", "start", "end", "genome_pos", "z", "report_state"]).issubset(plot_bins.columns)
    assert set(plot_bins["report_state"]) == {"dup", "del", "neutral"}
    assert plot_bins.loc[plot_bins["report_state"].eq("dup"), "z"].tolist() == [4.2, 3.8]
    assert plot_bins.loc[plot_bins["report_state"].eq("del"), "z"].tolist() == [-3.4, -3.1]
    assert 2_000_000 not in plot_bins.loc[plot_bins["chrom"].eq("chr2"), "start"].tolist()
    assert 9.9 not in plot_bins["z"].tolist()
    assert "smooth z trend" in svg
    assert "Branch A" not in svg
    assert "Branch B" not in svg
    assert "report_strong" not in svg
    assert "report_weak" not in svg
    assert "internal" not in svg
    assert "filtered" not in svg
    assert "rejected" not in svg
    assert "mask" not in svg.lower()
    assert "chr1" in svg
    assert '<circle cx="' in svg
    assert re.search(r'<circle[^>]+fill="#64748b"', svg)
    assert re.search(r'<circle[^>]+fill="#dc2626"', svg)
    assert re.search(r'<circle[^>]+fill="#2563eb"', svg)
    assert "dup" in svg
    assert "del" in svg
    assert "neutral bin" in svg


def test_build_cnv_plot_svg_requires_calibrated_z_without_fallback(tmp_path):
    bins = pd.DataFrame(
        {
            "chrom": ["chr1"],
            "bin_index": [0],
            "start": [0],
            "end": [1_000_000],
            "robust_z": [8.0],
            "raw_robust_z": [8.0],
            "normalized_signal": [8.0],
        }
    )

    with pytest.raises(ValueError, match="calibrated_z"):
        build_cnv_plot_svg(
            sample_id="Y1",
            bins_df=bins,
            branch_b_events_df=pd.DataFrame(),
            a_branch_df=pd.DataFrame(),
            output_svg=tmp_path / "Y1.final_cnv.svg",
            output_bins_tsv=tmp_path / "Y1.plot_bins.tsv",
        )
