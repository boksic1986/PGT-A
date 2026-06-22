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
            "calibrated_z": [0.1, 7.2, 6.8, -7.1, -6.6, pd.NA],
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
    assert plot_bins.loc[plot_bins["report_state"].eq("dup"), "z"].tolist() == [7.2, 6.8]
    assert plot_bins.loc[plot_bins["report_state"].eq("del"), "z"].tolist() == [-7.1, -6.6]
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
    assert set(["chrom", "start", "end", "genome_pos", "z", "copy_number", "report_state", "event_report_state", "copy_number_source"]).issubset(
        cn_bins.columns
    )
    assert "Copy number" in cn_svg
    assert "calibrated z" not in cn_svg.lower()
    assert 'width="2560"' in cn_svg
    assert "report CN trend" in cn_svg
    assert "smooth" not in cn_svg.lower()
    assert "<polyline" not in cn_svg
    assert cn_svg.count('class="report-cn-trend"') == 2
    assert 'fill="#1f2937"' not in cn_svg
    assert 'fill="#263244"' not in cn_svg
    assert 'fill="#202b3a"' not in cn_svg
    assert 'class="chrom-background"' not in cn_svg
    assert 'fill="#f8fafc"' not in cn_svg
    assert 'fill="#eef2f7"' not in cn_svg
    assert re.search(r'<line[^>]+class="report-cn-trend"[^>]+stroke="#1d4ed8"', cn_svg)
    assert re.search(r'<line[^>]+class="report-cn-trend"[^>]+stroke="#ef4444"', cn_svg)
    assert not re.search(r'<line[^>]+class="report-cn-trend"[^>]+stroke="#dc2626"', cn_svg)
    assert re.search(r'<circle[^>]+class="cn-bin-scatter"[^>]+fill="#1d4ed8"', cn_svg)
    assert re.search(r'<circle[^>]+class="cn-bin-scatter"[^>]+fill="#ef4444"', cn_svg)
    assert re.search(r'<circle[^>]+class="cn-bin-scatter"[^>]+fill="#64748b"', cn_svg)
    assert "Branch A" not in cn_svg
    assert "Branch B" not in cn_svg
    assert "internal" not in cn_svg
    assert "filtered" not in cn_svg
    assert "mask" not in cn_svg.lower()
    assert "cytoband" not in cn_svg.lower()
    assert "neutral bin" not in cn_svg
    assert "report CN trend" not in re.sub(r"<title>.*?</title>", "", cn_svg)
    assert cn_bins.loc[cn_bins["report_state"].eq("neutral"), "copy_number"].tolist() == pytest.approx([2.005], abs=0.001)
    assert cn_bins.loc[cn_bins["report_state"].eq("dup"), "copy_number"].tolist() == pytest.approx([2.36, 2.34], abs=0.001)
    assert cn_bins.loc[cn_bins["report_state"].eq("del"), "copy_number"].tolist() == pytest.approx([1.645, 1.67], abs=0.001)
    assert set(cn_bins.loc[cn_bins["report_state"].eq("dup"), "copy_number_source"]) == {
        "calibrated_z_mosaic30_cn_proxy"
    }
    assert set(cn_bins.loc[cn_bins["report_state"].eq("del"), "copy_number_source"]) == {
        "calibrated_z_mosaic30_cn_proxy"
    }
    assert cn_bins["copy_number"].max() <= 4.0
    assert cn_bins["copy_number"].min() >= 0.0


def test_copy_number_plot_v2_uses_structural_gap_blanks_and_50mb_ticks(tmp_path):
    starts = [idx * 1_000_000 for idx in range(60)]
    bins = pd.DataFrame(
        {
            "chrom": ["chr1"] * 60,
            "bin_index": list(range(60)),
            "start": starts,
            "end": [start + 1_000_000 for start in starts],
            "calibrated_z": [3.0 + (idx % 6) * 0.4 for idx in range(60)],
            "is_gap_centromere_telomere": [idx == 10 for idx in range(60)],
            "gap_centromere_telomere_overlap_fraction": [1.0 if idx == 10 else 0.0 for idx in range(60)],
            "is_near_centromere": [idx == 30 for idx in range(60)],
        }
    )
    events = pd.DataFrame(
        {
            "sample_id": ["Y1"],
            "chrom": ["chr1"],
            "start": [0],
            "end": [60_000_000],
            "state": ["gain"],
            "v2_report_layer_class": ["report_event"],
            "v2_report_visibility": ["report_strong_event"],
            "copy_number_estimate": [2.5],
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

    cn_svg = output_cn_svg.read_text(encoding="utf-8")
    cn_bins = pd.read_csv(output_cn_bins, sep="\t")
    gap_row = cn_bins.loc[cn_bins["start"].eq(30_000_000)].iloc[0]

    assert 'class="structure-gap-blank"' in cn_svg
    assert re.search(r'<rect[^>]+class="structure-gap-blank"[^>]+fill="#cbd5e1"', cn_svg)
    assert 'class="structure-gap-blank" x="' in cn_svg
    assert cn_svg.count('class="structure-gap-blank"') == 1
    assert 'fill="#0f172a" opacity="0.96"' not in cn_svg
    assert "50Mb" in cn_svg
    assert cn_svg.count('class="cn-bin-scatter"') == 59
    assert re.search(r'<circle[^>]+class="cn-bin-scatter"[^>]+r="1\.35"', cn_svg)
    assert cn_svg.count('class="report-cn-trend"') == 2
    assert re.search(r'<line[^>]+class="report-cn-trend"[^>]+stroke="#1d4ed8"', cn_svg)
    non_centromere_gap_row = cn_bins.loc[cn_bins["start"].eq(10_000_000)].iloc[0]
    assert non_centromere_gap_row["copy_number_source"] != "structure_gap_blank"
    assert pd.isna(gap_row["copy_number"])
    assert gap_row["copy_number_source"] == "structure_gap_blank"
    assert set(cn_bins.loc[~cn_bins["copy_number_source"].eq("structure_gap_blank"), "copy_number_source"]) == {
        "calibrated_z_mosaic30_cn_proxy"
    }
    assert cn_bins.loc[~cn_bins["copy_number_source"].eq("structure_gap_blank"), "copy_number"].nunique() > 1
    assert "#1d4ed8" in cn_svg
    assert 'class="chrom-background"' not in cn_svg
    assert 'fill="#f8fafc"' not in cn_svg
    assert 'fill="#eef2f7"' not in cn_svg
    assert 'class="chrom-separator"' in cn_svg
    assert not re.search(r'<circle[^>]+class="cn-bin-scatter"[^>]+fill="#ef4444"', cn_svg)


def test_copy_number_plot_uses_hg19_centromere_fallback_without_telomere_gap_shading(tmp_path):
    bins = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1"],
            "bin_index": [0, 1, 2],
            "start": [0, 121_500_000, 200_000_000],
            "end": [1_000_000, 122_500_000, 201_000_000],
            "calibrated_z": [0.1, 3.0, 0.2],
            "is_gap_centromere_telomere": [1, 1, 0],
            "gap_centromere_telomere_overlap_fraction": [1.0, 1.0, 0.0],
        }
    )
    events = pd.DataFrame(
        {
            "sample_id": ["Y1"],
            "chrom": ["chr1"],
            "start": [121_000_000],
            "end": [123_000_000],
            "state": ["gain"],
            "v2_report_layer_class": ["report_event"],
            "v2_report_visibility": ["report_strong_event"],
            "copy_number_estimate": [2.5],
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

    cn_svg = output_cn_svg.read_text(encoding="utf-8")
    cn_bins = pd.read_csv(output_cn_bins, sep="\t")
    telomere_like_row = cn_bins.loc[cn_bins["start"].eq(0)].iloc[0]
    centromere_row = cn_bins.loc[cn_bins["start"].eq(121_500_000)].iloc[0]

    assert cn_svg.count('class="structure-gap-blank"') == 1
    assert telomere_like_row["copy_number_source"] != "structure_gap_blank"
    assert centromere_row["copy_number_source"] == "structure_gap_blank"
    assert cn_svg.count('class="cn-bin-scatter"') == 2


def test_copy_number_scatter_uses_bin_cn_thresholds_independent_of_event_direction(tmp_path):
    bins = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1"],
            "bin_index": [0, 1, 2],
            "start": [0, 1_000_000, 2_000_000],
            "end": [1_000_000, 2_000_000, 3_000_000],
            "calibrated_z": [0.0, 7.0, 8.0],
        }
    )
    events = pd.DataFrame(
        {
            "sample_id": ["Y1"],
            "chrom": ["chr1"],
            "start": [1_000_000],
            "end": [3_000_000],
            "state": ["loss"],
            "v2_report_layer_class": ["report_event"],
            "v2_report_visibility": ["report_weak_event"],
            "copy_number_estimate": [1.5],
        }
    )

    output_cn_bins = tmp_path / "Y1.plot_bins_cn.tsv"
    build_cnv_plot_svg(
        sample_id="Y1",
        bins_df=bins,
        branch_b_events_df=events,
        a_branch_df=pd.DataFrame(),
        output_svg=tmp_path / "Y1.final_cnv.svg",
        output_bins_tsv=tmp_path / "Y1.plot_bins.tsv",
        output_copy_number_svg=tmp_path / "Y1.final_cnv_cn.svg",
        output_copy_number_bins_tsv=output_cn_bins,
    )

    cn_bins = pd.read_csv(output_cn_bins, sep="\t")
    event_bins = cn_bins.loc[cn_bins["event_report_state"].eq("del")].copy()

    assert event_bins["copy_number"].nunique() == 2
    assert event_bins["copy_number"].tolist() == pytest.approx([2.35, 2.4], abs=0.001)
    assert event_bins["report_state"].tolist() == ["dup", "dup"]
    assert set(event_bins["copy_number_source"]) == {"calibrated_z_mosaic30_cn_proxy"}
    assert event_bins["z"].tolist() == [7.0, 8.0]


def test_copy_number_scatter_covers_all_bins_and_clips_extreme_z_without_events(tmp_path):
    bins = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1", "chr1"],
            "bin_index": [0, 1, 2, 3],
            "start": [0, 1_000_000, 2_000_000, 3_000_000],
            "end": [1_000_000, 2_000_000, 3_000_000, 4_000_000],
            "calibrated_z": [-8.0, 0.0, 8.0, 50.0],
        }
    )
    output_cn_svg = tmp_path / "Y1.final_cnv_cn.svg"
    output_cn_bins = tmp_path / "Y1.plot_bins_cn.tsv"

    build_cnv_plot_svg(
        sample_id="Y1",
        bins_df=bins,
        branch_b_events_df=pd.DataFrame(),
        a_branch_df=pd.DataFrame(),
        output_svg=tmp_path / "Y1.final_cnv.svg",
        output_bins_tsv=tmp_path / "Y1.plot_bins.tsv",
        output_copy_number_svg=output_cn_svg,
        output_copy_number_bins_tsv=output_cn_bins,
    )

    cn_svg = output_cn_svg.read_text(encoding="utf-8")
    cn_bins = pd.read_csv(output_cn_bins, sep="\t")

    assert cn_svg.count('class="cn-bin-scatter"') == 4
    assert re.search(r'<circle[^>]+class="cn-bin-scatter"[^>]+fill="#ef4444"', cn_svg)
    assert re.search(r'<circle[^>]+class="cn-bin-scatter"[^>]+fill="#64748b"', cn_svg)
    assert re.search(r'<circle[^>]+class="cn-bin-scatter"[^>]+fill="#1d4ed8"', cn_svg)
    assert cn_bins["copy_number"].tolist() == pytest.approx([1.6, 2.0, 2.4, 4.0], abs=0.001)
    assert cn_bins["report_state"].tolist() == ["del", "neutral", "dup", "dup"]
    assert set(cn_bins["event_report_state"]) == {"neutral"}
    assert set(cn_bins["copy_number_source"]) == {"calibrated_z_mosaic30_cn_proxy"}
    assert 'class="report-cn-trend"' not in cn_svg


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
    assert cn_bins["copy_number"].tolist() == pytest.approx([2.155, 2.16], abs=0.001)
    assert set(cn_bins["copy_number_source"]) == {"calibrated_z_mosaic30_cn_proxy"}
    assert output_cn_svg.read_text(encoding="utf-8").count('class="report-cn-trend"') == 1

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
