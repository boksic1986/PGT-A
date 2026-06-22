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
            "copy_number_estimate": [2.64, 1.42, 2.12],
            "a_ratio": [0.32, -0.29, 0.06],
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
    output_cn_svg = tmp_path / "Y1.final_cnv_cn.svg"
    output_cn_bins = tmp_path / "Y1.plot_bins_cn.tsv"
    build_cnv_plot_svg(
        sample_id="Y1",
        bins_df=bins,
        branch_b_events_df=events,
        a_branch_df=a_branch,
        output_svg=output_svg,
        output_bins_tsv=output_bins,
        output_copy_number_svg=output_cn_svg,
        output_copy_number_bins_tsv=output_cn_bins,
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
    assert "report z trend" in svg
    assert "smooth z trend" not in svg
    assert "<polyline" not in svg
    assert svg.count('class="report-z-trend"') == 2
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
    assert re.search(r'<circle[^>]+fill="#facc15"', svg)
    assert re.search(r'<circle[^>]+fill="#2563eb"', svg)
    assert re.search(r'<line[^>]+class="report-z-trend"[^>]+stroke="#dc2626"', svg)
    assert "dup" in svg
    assert "del" in svg
    assert "neutral bin" in svg

    cn_svg = output_cn_svg.read_text(encoding="utf-8")
    cn_bins = pd.read_csv(output_cn_bins, sep="\t")
    assert set(["chrom", "start", "end", "genome_pos", "copy_number", "report_state", "copy_number_source"]).issubset(
        cn_bins.columns
    )
    assert "Copy number" in cn_svg
    assert "calibrated z" not in cn_svg.lower()
    assert "report CN trend" in cn_svg
    assert "smooth" not in cn_svg.lower()
    assert "<polyline" not in cn_svg
    assert cn_svg.count('class="report-cn-trend"') == 2
    assert re.search(r'<line[^>]+class="report-cn-trend"[^>]+stroke="#dc2626"', cn_svg)
    assert re.search(r'<circle[^>]+fill="#facc15"', cn_svg)
    assert re.search(r'<circle[^>]+fill="#2563eb"', cn_svg)
    assert re.search(r'<circle[^>]+fill="#64748b"', cn_svg)
    assert "Branch A" not in cn_svg
    assert "Branch B" not in cn_svg
    assert "internal" not in cn_svg
    assert "filtered" not in cn_svg
    assert "mask" not in cn_svg.lower()
    assert cn_bins.loc[cn_bins["report_state"].eq("neutral"), "copy_number"].tolist() == [2.0]
    assert cn_bins.loc[cn_bins["report_state"].eq("dup"), "copy_number"].tolist() == [2.64, 2.64]
    assert cn_bins.loc[cn_bins["report_state"].eq("del"), "copy_number"].tolist() == [1.42, 1.42]
    assert set(cn_bins.loc[cn_bins["report_state"].eq("dup"), "copy_number_source"]) == {"copy_number_estimate"}
    assert set(cn_bins.loc[cn_bins["report_state"].eq("del"), "copy_number_source"]) == {"copy_number_estimate"}


def test_copy_number_plot_falls_back_to_a_ratio_and_refuses_missing_cn(tmp_path):
    bins = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "bin_index": [0, 1],
            "start": [0, 1_000_000],
            "end": [1_000_000, 2_000_000],
            "calibrated_z": [3.1, 3.2],
            "normalized_signal": [99.0, 100.0],
        }
    )
    events = pd.DataFrame(
        {
            "sample_id": ["Y1"],
            "chrom": ["chr1"],
            "start": [0],
            "end": [2_000_000],
            "state": ["gain"],
            "v2_report_layer_class": ["report_event"],
            "v2_report_visibility": ["report_weak_event"],
            "a_ratio": [0.25],
        }
    )

    output_cn_svg = tmp_path / "Y1.final_cnv_cn.svg"
    output_cn_bins = tmp_path / "Y1.plot_bins_cn.tsv"
    build_cnv_plot_svg(
        sample_id="Y1",
        bins_df=bins,
        branch_b_events_df=events,
        a_branch_df=pd.DataFrame(),
        output_svg=tmp_path / "Y1.final_cnv.svg",
        output_bins_tsv=tmp_path / "Y1.plot_bins.tsv",
        output_copy_number_svg=output_cn_svg,
        output_copy_number_bins_tsv=output_cn_bins,
    )
    cn_bins = pd.read_csv(output_cn_bins, sep="\t")
    assert cn_bins["copy_number"].tolist() == [2.5, 2.5]
    assert set(cn_bins["copy_number_source"]) == {"a_ratio"}

    missing_cn_events = events.drop(columns=["a_ratio"])
    with pytest.raises(ValueError, match="copy_number_estimate|a_ratio"):
        build_cnv_plot_svg(
            sample_id="Y1",
            bins_df=bins,
            branch_b_events_df=missing_cn_events,
            a_branch_df=pd.DataFrame(),
            output_svg=tmp_path / "Y1.no_cn.final_cnv.svg",
            output_bins_tsv=tmp_path / "Y1.no_cn.plot_bins.tsv",
            output_copy_number_svg=tmp_path / "Y1.no_cn.final_cnv_cn.svg",
            output_copy_number_bins_tsv=tmp_path / "Y1.no_cn.plot_bins_cn.tsv",
        )


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
