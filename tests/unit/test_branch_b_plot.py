import re
from pathlib import Path

import pandas as pd
import pytest

from pgta.predict.branch_b.plot import build_cnv_plot_svg


def _norm_signal(cpm):
    import numpy as np

    return np.log1p(cpm)


def _ref_bins(rows):
    return pd.DataFrame(
        {
            "chrom": [row[0] for row in rows],
            "bin_index": [row[1] for row in rows],
            "start": [row[2] for row in rows],
            "end": [row[3] for row in rows],
            "ref_median": [_norm_signal(row[4]) for row in rows],
        }
    )


def test_build_cnv_plot_svg_uses_calibrated_z_and_writes_plot_bins(tmp_path):
    bins = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1", "chr2", "chr2", "chr2"],
            "bin_index": [0, 1, 2, 0, 1, 2],
            "start": [0, 1_000_000, 2_000_000, 0, 1_000_000, 2_000_000],
            "end": [1_000_000, 2_000_000, 3_000_000, 1_000_000, 2_000_000, 3_000_000],
            "calibrated_z": [0.1, 7.2, 6.8, -7.1, -6.6, pd.NA],
            "normalized_signal": [
                _norm_signal(100),
                _norm_signal(118),
                _norm_signal(119),
                _norm_signal(84),
                _norm_signal(83),
                _norm_signal(100),
            ],
            "robust_z": [9.9, 9.9, 9.9, 9.9, 9.9, 9.9],
            "mask_label": ["pass", "pass", "soft", "pass", "hard", "pass"],
            "region_risk_class": ["clean", "clean", "moderate", "clean", "high", "clean"],
        }
    )
    ref_bins = _ref_bins(
        [
            ("chr1", 0, 0, 1_000_000, 100),
            ("chr1", 1, 1_000_000, 2_000_000, 100),
            ("chr1", 2, 2_000_000, 3_000_000, 100),
            ("chr2", 0, 0, 1_000_000, 100),
            ("chr2", 1, 1_000_000, 2_000_000, 100),
            ("chr2", 2, 2_000_000, 3_000_000, 100),
        ]
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
        ref_bins_df=ref_bins,
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
    assert set(
        [
            "chrom",
            "start",
            "end",
            "genome_pos",
            "z",
            "raw_log2r",
            "log2r",
            "copy_number",
            "copy_number_centering_log2_shift",
            "copy_number_source",
            "cn_scatter_state",
            "report_state",
            "event_report_state",
            "is_structure_gap_blank",
        ]
    ).issubset(cn_bins.columns)
    assert "Copy number" in cn_svg
    assert "calibrated z" not in cn_svg.lower()
    assert "ratio-derived" in cn_svg
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
    assert re.search(r'<rect[^>]+class="report-cn-region"[^>]+fill="#1d4ed8"', cn_svg)
    assert re.search(r'<rect[^>]+class="report-cn-region"[^>]+fill="#ef4444"', cn_svg)
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
    assert cn_bins.loc[cn_bins["cn_scatter_state"].eq("neutral"), "copy_number"].tolist() == pytest.approx([2.0], abs=0.01)
    assert cn_bins.loc[cn_bins["cn_scatter_state"].eq("dup"), "copy_number"].tolist() == pytest.approx([2.36, 2.38], abs=0.01)
    assert cn_bins.loc[cn_bins["cn_scatter_state"].eq("del"), "copy_number"].tolist() == pytest.approx([1.68, 1.66], abs=0.01)
    assert set(cn_bins.loc[cn_bins["cn_scatter_state"].eq("dup"), "copy_number_source"]) == {
        "normalized_signal_ref_median_log2r_autosome_centered"
    }
    assert set(cn_bins.loc[cn_bins["cn_scatter_state"].eq("del"), "copy_number_source"]) == {
        "normalized_signal_ref_median_log2r_autosome_centered"
    }
    assert cn_bins["copy_number_centering_log2_shift"].dropna().abs().max() == pytest.approx(0.0)
    assert "calibrated_z_mosaic30_cn_proxy" not in set(cn_bins["copy_number_source"])


def test_copy_number_plot_v2_uses_structural_gap_blanks_and_50mb_ticks(tmp_path):
    starts = [idx * 1_000_000 for idx in range(60)]
    bins = pd.DataFrame(
        {
            "chrom": ["chr1"] * 60,
            "bin_index": list(range(60)),
            "start": starts,
            "end": [start + 1_000_000 for start in starts],
            "calibrated_z": [3.0 + (idx % 6) * 0.4 for idx in range(60)],
            "normalized_signal": [_norm_signal(130 + (idx % 6)) for idx in range(60)],
            "is_gap_centromere_telomere": [idx == 10 for idx in range(60)],
            "gap_centromere_telomere_overlap_fraction": [1.0 if idx == 10 else 0.0 for idx in range(60)],
            "is_near_centromere": [idx == 30 for idx in range(60)],
        }
    )
    ref_bins = _ref_bins([("chr1", idx, start, start + 1_000_000, 100) for idx, start in enumerate(starts)])
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
        ref_bins_df=ref_bins,
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
        "normalized_signal_ref_median_log2r_autosome_centered"
    }
    assert cn_bins.loc[~cn_bins["copy_number_source"].eq("structure_gap_blank"), "copy_number"].nunique() > 1
    assert "#1d4ed8" in cn_svg
    assert 'class="chrom-background"' not in cn_svg
    assert 'fill="#f8fafc"' not in cn_svg


def test_copy_number_plot_recenters_ratio_cn_to_autosomal_median(tmp_path):
    bins = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1", "chr2", "chr2"],
            "bin_index": [0, 1, 2, 0, 1],
            "start": [0, 1_000_000, 2_000_000, 0, 1_000_000],
            "end": [1_000_000, 2_000_000, 3_000_000, 1_000_000, 2_000_000],
            "calibrated_z": [0.0, 6.0, -6.0, 0.0, 0.0],
            "normalized_signal": [
                _norm_signal(85),
                _norm_signal(102),
                _norm_signal(68),
                _norm_signal(85),
                _norm_signal(85),
            ],
        }
    )
    ref_bins = _ref_bins(
        [
            ("chr1", 0, 0, 1_000_000, 100),
            ("chr1", 1, 1_000_000, 2_000_000, 100),
            ("chr1", 2, 2_000_000, 3_000_000, 100),
            ("chr2", 0, 0, 1_000_000, 100),
            ("chr2", 1, 1_000_000, 2_000_000, 100),
        ]
    )

    output_cn_svg = tmp_path / "S1.final_cnv_cn.svg"
    output_cn_bins = tmp_path / "S1.plot_bins_cn.tsv"
    build_cnv_plot_svg(
        sample_id="S1",
        bins_df=bins,
        branch_b_events_df=pd.DataFrame(),
        a_branch_df=pd.DataFrame(),
        ref_bins_df=ref_bins,
        output_svg=tmp_path / "S1.final_cnv.svg",
        output_bins_tsv=tmp_path / "S1.plot_bins.tsv",
        output_copy_number_svg=output_cn_svg,
        output_copy_number_bins_tsv=output_cn_bins,
    )

    cn_bins = pd.read_csv(output_cn_bins, sep="\t")
    assert cn_bins["copy_number"].median() == pytest.approx(2.0, abs=0.01)
    assert cn_bins["copy_number_centering_log2_shift"].dropna().iloc[0] < 0.0
    assert set(cn_bins["copy_number_source"]) == {"normalized_signal_ref_median_log2r_autosome_centered"}


def test_copy_number_plot_uses_hg19_centromere_fallback_without_telomere_gap_shading(tmp_path):
    bins = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1"],
            "bin_index": [0, 1, 2],
            "start": [0, 121_500_000, 200_000_000],
            "end": [1_000_000, 122_500_000, 201_000_000],
            "calibrated_z": [0.1, 3.0, 0.2],
            "normalized_signal": [_norm_signal(100), _norm_signal(130), _norm_signal(101)],
            "is_gap_centromere_telomere": [1, 1, 0],
            "gap_centromere_telomere_overlap_fraction": [1.0, 1.0, 0.0],
        }
    )
    ref_bins = _ref_bins(
        [
            ("chr1", 0, 0, 1_000_000, 100),
            ("chr1", 1, 121_500_000, 122_500_000, 100),
            ("chr1", 2, 200_000_000, 201_000_000, 100),
        ]
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
        ref_bins_df=ref_bins,
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


def test_copy_number_scatter_uses_ratio_derived_bin_cn_independent_of_event_direction(tmp_path):
    bins = pd.DataFrame(
        {
            "chrom": ["chr1"] * 7,
            "bin_index": list(range(7)),
            "start": [idx * 1_000_000 for idx in range(7)],
            "end": [(idx + 1) * 1_000_000 for idx in range(7)],
            "calibrated_z": [0.0, 0.0, 0.0, 0.0, 0.0, 7.0, 8.0],
            "normalized_signal": [_norm_signal(value) for value in [100, 100, 100, 100, 100, 118, 120]],
        }
    )
    ref_bins = _ref_bins(
        [
            ("chr1", idx, idx * 1_000_000, (idx + 1) * 1_000_000, 100)
            for idx in range(7)
        ]
    )
    events = pd.DataFrame(
        {
            "sample_id": ["Y1"],
            "chrom": ["chr1"],
            "start": [5_000_000],
            "end": [7_000_000],
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
        ref_bins_df=ref_bins,
        output_svg=tmp_path / "Y1.final_cnv.svg",
        output_bins_tsv=tmp_path / "Y1.plot_bins.tsv",
        output_copy_number_svg=tmp_path / "Y1.final_cnv_cn.svg",
        output_copy_number_bins_tsv=output_cn_bins,
    )

    cn_bins = pd.read_csv(output_cn_bins, sep="\t")
    event_bins = cn_bins.loc[cn_bins["event_report_state"].eq("del")].copy()

    assert event_bins["copy_number"].nunique() == 2
    assert event_bins["copy_number"].tolist() == pytest.approx([2.36, 2.40], abs=0.01)
    assert event_bins["cn_scatter_state"].tolist() == ["dup", "dup"]
    assert event_bins["event_report_state"].tolist() == ["del", "del"]
    assert set(event_bins["copy_number_source"]) == {"normalized_signal_ref_median_log2r_autosome_centered"}
    assert event_bins["z"].tolist() == [7.0, 8.0]


def test_copy_number_scatter_covers_all_bins_and_does_not_clip_tsv_values(tmp_path):
    bins = pd.DataFrame(
        {
            "chrom": ["chr1"] * 5,
            "bin_index": [0, 1, 2, 3, 4],
            "start": [0, 1_000_000, 2_000_000, 3_000_000, 4_000_000],
            "end": [1_000_000, 2_000_000, 3_000_000, 4_000_000, 5_000_000],
            "calibrated_z": [-8.0, 0.0, 0.0, 8.0, 50.0],
            "normalized_signal": [_norm_signal(value) for value in [80, 100, 100, 120, 300]],
        }
    )
    ref_bins = _ref_bins(
        [
            ("chr1", idx, idx * 1_000_000, (idx + 1) * 1_000_000, 100)
            for idx in range(5)
        ]
    )
    output_cn_svg = tmp_path / "Y1.final_cnv_cn.svg"
    output_cn_bins = tmp_path / "Y1.plot_bins_cn.tsv"

    build_cnv_plot_svg(
        sample_id="Y1",
        bins_df=bins,
        branch_b_events_df=pd.DataFrame(),
        a_branch_df=pd.DataFrame(),
        ref_bins_df=ref_bins,
        output_svg=tmp_path / "Y1.final_cnv.svg",
        output_bins_tsv=tmp_path / "Y1.plot_bins.tsv",
        output_copy_number_svg=output_cn_svg,
        output_copy_number_bins_tsv=output_cn_bins,
    )

    cn_svg = output_cn_svg.read_text(encoding="utf-8")
    cn_bins = pd.read_csv(output_cn_bins, sep="\t")

    assert cn_svg.count('class="cn-bin-scatter"') == 5
    assert re.search(r'<circle[^>]+class="cn-bin-scatter"[^>]+fill="#ef4444"', cn_svg)
    assert re.search(r'<circle[^>]+class="cn-bin-scatter"[^>]+fill="#64748b"', cn_svg)
    assert re.search(r'<circle[^>]+class="cn-bin-scatter"[^>]+fill="#1d4ed8"', cn_svg)
    assert cn_bins["copy_number"].tolist() == pytest.approx([1.6, 2.0, 2.0, 2.4, 6.0], abs=0.01)
    assert cn_bins["cn_scatter_state"].tolist() == ["del", "neutral", "neutral", "dup", "dup"]
    assert set(cn_bins["event_report_state"]) == {"neutral"}
    assert set(cn_bins["copy_number_source"]) == {"normalized_signal_ref_median_log2r_autosome_centered"}
    assert 'data-copy-number-out-of-range="true"' in cn_svg
    assert 'class="report-cn-trend"' not in cn_svg


def test_copy_number_plot_uses_ref_median_and_refuses_z_only_cn(tmp_path):
    bins = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "bin_index": [0, 1],
            "start": [0, 1_000_000],
            "end": [1_000_000, 2_000_000],
            "calibrated_z": [3.1, 3.2],
            "normalized_signal": [_norm_signal(112), _norm_signal(114)],
        }
    )
    ref_bins = _ref_bins(
        [
            ("chr1", 0, 0, 1_000_000, 100),
            ("chr1", 1, 1_000_000, 2_000_000, 100),
        ]
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
        ref_bins_df=ref_bins,
        output_svg=tmp_path / "Y1.final_cnv.svg",
        output_bins_tsv=tmp_path / "Y1.plot_bins.tsv",
        output_copy_number_svg=output_cn_svg,
        output_copy_number_bins_tsv=output_cn_bins,
    )
    cn_bins = pd.read_csv(output_cn_bins, sep="\t")
    assert cn_bins["copy_number"].tolist() == pytest.approx([1.98, 2.02], abs=0.01)
    assert set(cn_bins["copy_number_source"]) == {"normalized_signal_ref_median_log2r_autosome_centered"}
    assert output_cn_svg.read_text(encoding="utf-8").count('class="report-cn-trend"') == 1

    with pytest.raises(ValueError, match="reference bin medians|log2r|copy_number"):
        build_cnv_plot_svg(
            sample_id="Y1",
            bins_df=bins,
            branch_b_events_df=events,
            a_branch_df=pd.DataFrame(),
            output_svg=tmp_path / "Y1.no_cn.final_cnv.svg",
            output_bins_tsv=tmp_path / "Y1.no_cn.plot_bins.tsv",
            output_copy_number_svg=tmp_path / "Y1.no_cn.final_cnv_cn.svg",
            output_copy_number_bins_tsv=tmp_path / "Y1.no_cn.plot_bins_cn.tsv",
        )


def test_copy_number_event_support_counts_use_ratio_cn_not_z(tmp_path):
    bins = pd.DataFrame(
        {
            "chrom": ["chr1"] * 9,
            "bin_index": list(range(9)),
            "start": [idx * 1_000_000 for idx in range(9)],
            "end": [(idx + 1) * 1_000_000 for idx in range(9)],
            "calibrated_z": [-8.0, -7.0, -6.0, -5.0, -4.0, 12.0, 8.0, 7.0, 6.0],
            "normalized_signal": [_norm_signal(value) for value in [100, 100, 100, 100, 100, 100, 180, 200, 190]],
            "is_near_centromere": [False, False, False, False, False, True, False, False, False],
        }
    )
    ref_bins = _ref_bins(
        [
            ("chr1", idx, idx * 1_000_000, (idx + 1) * 1_000_000, 100)
            for idx in range(9)
        ]
    )
    events = pd.DataFrame(
        {
            "sample_id": ["Y1"],
            "chrom": ["chr1"],
            "start": [5_000_000],
            "end": [9_000_000],
            "state": ["gain"],
            "v2_report_layer_class": ["report_event"],
            "v2_report_visibility": ["report_strong_event"],
            "copy_number_estimate": [2.36],
        }
    )
    support_tsv = tmp_path / "Y1.plot_event_support.tsv"

    build_cnv_plot_svg(
        sample_id="Y1",
        bins_df=bins,
        branch_b_events_df=events,
        a_branch_df=pd.DataFrame(),
        ref_bins_df=ref_bins,
        output_svg=tmp_path / "Y1.final_cnv.svg",
        output_bins_tsv=tmp_path / "Y1.plot_bins.tsv",
        output_copy_number_svg=tmp_path / "Y1.final_cnv_cn.svg",
        output_copy_number_bins_tsv=tmp_path / "Y1.plot_bins_cn.tsv",
        output_copy_number_event_support_tsv=support_tsv,
    )

    support = pd.read_csv(support_tsv, sep="\t")
    row = support.iloc[0]
    assert row["valid_bin_count"] == 3
    assert row["centromere_gap_bin_count"] == 1
    assert row["cn_support_bin_count"] == 3
    assert row["cn_same_direction_fraction"] == pytest.approx(1.0, abs=0.001)
    assert row["cn_direction_consistency_status"] == "CN_DIRECTION_SUPPORTED"
    assert row["z_support_bin_count"] == 3


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


def test_copy_number_scatter_is_sex_aware_for_xy_chrx_and_chry(tmp_path):
    bins = pd.DataFrame(
        {
            "chrom": ["chr1", "chrX", "chrY"],
            "bin_index": [0, 0, 0],
            "start": [0, 0, 0],
            "end": [1_000_000, 1_000_000, 1_000_000],
            "calibrated_z": [0.0, -0.1, 0.1],
            "normalized_signal": [_norm_signal(100), _norm_signal(50), _norm_signal(50)],
        }
    )
    ref_bins = _ref_bins(
        [
            ("chr1", 0, 0, 1_000_000, 100),
            ("chrX", 0, 0, 1_000_000, 100),
            ("chrY", 0, 0, 1_000_000, 50),
        ]
    )
    gender = pd.DataFrame(
        {
            "sample_id": ["XY1"],
            "sex_call": ["XY"],
            "predict_gender": ["M"],
            "sex_call_source": ["unit_test"],
        }
    )
    output_cn_svg = tmp_path / "XY1.final_cnv_cn.svg"
    output_cn_bins = tmp_path / "XY1.plot_bins_cn.tsv"

    build_cnv_plot_svg(
        sample_id="XY1",
        bins_df=bins,
        branch_b_events_df=pd.DataFrame(),
        a_branch_df=pd.DataFrame(),
        gender_df=gender,
        ref_bins_df=ref_bins,
        output_svg=tmp_path / "XY1.final_cnv.svg",
        output_bins_tsv=tmp_path / "XY1.plot_bins.tsv",
        output_copy_number_svg=output_cn_svg,
        output_copy_number_bins_tsv=output_cn_bins,
    )

    cn_bins = pd.read_csv(output_cn_bins, sep="\t")
    assert set(
        [
            "sex_call",
            "expected_copy_number",
            "copy_number_delta",
            "cn_scatter_state_sex_aware",
            "sex_chrom_region_class",
            "copy_number_interpretation_status",
            "event_layer",
        ]
    ).issubset(cn_bins.columns)
    chr_x = cn_bins.loc[cn_bins["chrom"].eq("chrX")].iloc[0]
    chr_y = cn_bins.loc[cn_bins["chrom"].eq("chrY")].iloc[0]
    assert chr_x["expected_copy_number"] == pytest.approx(1.0)
    assert chr_y["expected_copy_number"] == pytest.approx(1.0)
    assert chr_x["copy_number"] == pytest.approx(1.0, abs=0.02)
    assert chr_y["copy_number"] == pytest.approx(1.0, abs=0.02)
    assert chr_x["cn_scatter_state_sex_aware"] == "neutral"
    assert chr_y["cn_scatter_state_sex_aware"] == "neutral"
    assert chr_x["report_state"] == "neutral"
    assert chr_y["report_state"] == "neutral"
    assert chr_x["copy_number_interpretation_status"] == "sex_aware_interpretable"
    assert chr_y["copy_number_interpretation_status"] == "sex_aware_interpretable"


def test_copy_number_scatter_does_not_turn_xx_absent_chry_into_deletion(tmp_path):
    bins = pd.DataFrame(
        {
            "chrom": ["chr1", "chrY"],
            "bin_index": [0, 0],
            "start": [0, 0],
            "end": [1_000_000, 1_000_000],
            "calibrated_z": [0.0, -0.1],
            "normalized_signal": [_norm_signal(100), _norm_signal(0.01)],
        }
    )
    ref_bins = _ref_bins(
        [
            ("chr1", 0, 0, 1_000_000, 100),
            ("chrY", 0, 0, 1_000_000, 50),
        ]
    )
    gender = pd.DataFrame({"sample_id": ["XX1"], "sex_call": ["XX"], "predict_gender": ["F"]})
    output_cn_bins = tmp_path / "XX1.plot_bins_cn.tsv"

    build_cnv_plot_svg(
        sample_id="XX1",
        bins_df=bins,
        branch_b_events_df=pd.DataFrame(),
        a_branch_df=pd.DataFrame(),
        gender_df=gender,
        ref_bins_df=ref_bins,
        output_svg=tmp_path / "XX1.final_cnv.svg",
        output_bins_tsv=tmp_path / "XX1.plot_bins.tsv",
        output_copy_number_svg=tmp_path / "XX1.final_cnv_cn.svg",
        output_copy_number_bins_tsv=output_cn_bins,
    )

    cn_bins = pd.read_csv(output_cn_bins, sep="\t")
    chr_y = cn_bins.loc[cn_bins["chrom"].eq("chrY")].iloc[0]
    assert chr_y["expected_copy_number"] == pytest.approx(0.0)
    assert chr_y["cn_scatter_state_sex_aware"] == "neutral"
    assert chr_y["report_state"] == "neutral"
    assert chr_y["copy_number_interpretation_status"] == "sex_aware_absent_expected"


def test_copy_number_scatter_marks_chry_zero_ref_as_not_interpretable(tmp_path):
    bins = pd.DataFrame(
        {
            "chrom": ["chr1", "chrY"],
            "bin_index": [0, 0],
            "start": [0, 0],
            "end": [1_000_000, 1_000_000],
            "calibrated_z": [0.0, 0.1],
            "normalized_signal": [_norm_signal(100), _norm_signal(50)],
        }
    )
    ref_bins = _ref_bins(
        [
            ("chr1", 0, 0, 1_000_000, 100),
            ("chrY", 0, 0, 1_000_000, 0),
        ]
    )
    gender = pd.DataFrame({"sample_id": ["XY1"], "sex_call": ["XY"], "predict_gender": ["M"]})
    output_cn_bins = tmp_path / "XY1.plot_bins_cn.tsv"

    build_cnv_plot_svg(
        sample_id="XY1",
        bins_df=bins,
        branch_b_events_df=pd.DataFrame(),
        a_branch_df=pd.DataFrame(),
        gender_df=gender,
        ref_bins_df=ref_bins,
        output_svg=tmp_path / "XY1.final_cnv.svg",
        output_bins_tsv=tmp_path / "XY1.plot_bins.tsv",
        output_copy_number_svg=tmp_path / "XY1.final_cnv_cn.svg",
        output_copy_number_bins_tsv=output_cn_bins,
    )

    cn_bins = pd.read_csv(output_cn_bins, sep="\t")
    chr_y = cn_bins.loc[cn_bins["chrom"].eq("chrY")].iloc[0]
    assert pd.isna(chr_y["copy_number"])
    assert chr_y["cn_scatter_state_sex_aware"] == "neutral"
    assert chr_y["report_state"] == "neutral"
    assert chr_y["copy_number_interpretation_status"] == "sex_chrom_ref_ratio_not_interpretable"


def test_copy_number_scatter_marks_chry_median_zero_ref_as_not_interpretable(tmp_path):
    bins = pd.DataFrame(
        {
            "chrom": ["chr1", "chrY", "chrY", "chrY"],
            "bin_index": [0, 0, 1, 2],
            "start": [0, 0, 1_000_000, 2_000_000],
            "end": [1_000_000, 1_000_000, 2_000_000, 3_000_000],
            "calibrated_z": [0.0, 0.1, 0.2, 0.1],
            "normalized_signal": [_norm_signal(100), _norm_signal(50), _norm_signal(50), _norm_signal(50)],
        }
    )
    ref_bins = _ref_bins(
        [
            ("chr1", 0, 0, 1_000_000, 100),
            ("chrY", 0, 0, 1_000_000, 0),
            ("chrY", 1, 1_000_000, 2_000_000, 0),
            ("chrY", 2, 2_000_000, 3_000_000, 5),
        ]
    )
    gender = pd.DataFrame({"sample_id": ["XY1"], "sex_call": ["XY"], "predict_gender": ["M"]})
    output_cn_bins = tmp_path / "XY1.plot_bins_cn.tsv"

    build_cnv_plot_svg(
        sample_id="XY1",
        bins_df=bins,
        branch_b_events_df=pd.DataFrame(),
        a_branch_df=pd.DataFrame(),
        gender_df=gender,
        ref_bins_df=ref_bins,
        output_svg=tmp_path / "XY1.final_cnv.svg",
        output_bins_tsv=tmp_path / "XY1.plot_bins.tsv",
        output_copy_number_svg=tmp_path / "XY1.final_cnv_cn.svg",
        output_copy_number_bins_tsv=output_cn_bins,
    )

    cn_bins = pd.read_csv(output_cn_bins, sep="\t")
    chr_y = cn_bins.loc[cn_bins["chrom"].eq("chrY")]
    assert chr_y["copy_number"].isna().all()
    assert set(chr_y["cn_scatter_state_sex_aware"]) == {"neutral"}
    assert set(chr_y["copy_number_interpretation_status"]) == {"sex_chrom_ref_ratio_not_interpretable"}


def test_branch_s_report_review_event_is_overlaid_in_combined_genome_plot(tmp_path):
    bins = pd.DataFrame(
        {
            "chrom": ["chrX", "chrX", "chrX"],
            "bin_index": [0, 1, 2],
            "start": [0, 1_000_000, 2_000_000],
            "end": [1_000_000, 2_000_000, 3_000_000],
            "calibrated_z": [-0.1, -0.2, -0.1],
            "normalized_signal": [_norm_signal(85), _norm_signal(84), _norm_signal(86)],
        }
    )
    ref_bins = _ref_bins(
        [
            ("chrX", 0, 0, 1_000_000, 100),
            ("chrX", 1, 1_000_000, 2_000_000, 100),
            ("chrX", 2, 2_000_000, 3_000_000, 100),
        ]
    )
    gender = pd.DataFrame({"sample_id": ["G3"], "sex_call": ["XX"], "predict_gender": ["F"]})
    branch_s_summary = pd.DataFrame(
        {
            "sample_id": ["G3"],
            "sex_call": ["XX"],
            "sca_candidate_state": ["X_LOSS"],
            "sca_report_layer_class": ["sca_report_review_event"],
            "sca_report_layer_reason": ["sca_review_weak_report_visible_with_branch_a_support"],
        }
    )
    branch_s_scores = pd.DataFrame(
        {
            "sample_id": ["G3"],
            "sca_state": ["X_LOSS"],
            "state_score": [24.16],
            "state_score_reason": ["branch_a_candidate_zscore_sex_call_compatible_uncorroborated_review"],
        }
    )
    branch_s_evidence = pd.DataFrame(
        {
            "sample_id": ["G3"],
            "region_class": ["X_NONPAR"],
            "chrom": ["chrX"],
            "start": [0],
            "end": [3_000_000],
            "bin_count": [3],
            "median_calibrated_z": [-0.1],
        }
    )
    output_svg = tmp_path / "G3.final_cnv.svg"
    output_cn_svg = tmp_path / "G3.final_cnv_cn.svg"
    output_cn_bins = tmp_path / "G3.plot_bins_cn.tsv"

    build_cnv_plot_svg(
        sample_id="G3",
        bins_df=bins,
        branch_b_events_df=pd.DataFrame(),
        a_branch_df=pd.DataFrame(),
        gender_df=gender,
        branch_s_summary_df=branch_s_summary,
        branch_s_scores_df=branch_s_scores,
        branch_s_evidence_df=branch_s_evidence,
        ref_bins_df=ref_bins,
        output_svg=output_svg,
        output_bins_tsv=tmp_path / "G3.plot_bins.tsv",
        output_copy_number_svg=output_cn_svg,
        output_copy_number_bins_tsv=output_cn_bins,
    )

    svg = output_svg.read_text(encoding="utf-8")
    cn_svg = output_cn_svg.read_text(encoding="utf-8")
    cn_bins = pd.read_csv(output_cn_bins, sep="\t")
    assert 'class="branch-s-review-region"' in svg
    assert 'class="branch-s-cn-region"' in cn_svg
    assert "event-level score=24.160" in svg
    assert set(cn_bins["event_layer"]) == {"branch_s_review"}
