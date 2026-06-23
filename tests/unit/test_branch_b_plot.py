import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pgta.predict.branch_b.plot import build_cnv_plot_svg


def _norm_signal(cpm):
    import numpy as np

    return np.log1p(cpm)


def _ref_bins(rows):
    payload = {
        "chrom": [row[0] for row in rows],
        "bin_index": [row[1] for row in rows],
        "start": [row[2] for row in rows],
        "end": [row[3] for row in rows],
        "ref_median": [_norm_signal(row[4]) for row in rows],
        "ref_mad": [row[5] if len(row) > 5 else 0.1 for row in rows],
    }
    if any(len(row) > 6 for row in rows):
        payload["ref_group"] = [row[6] if len(row) > 6 else "mixed" for row in rows]
    return pd.DataFrame(payload)


def _ref_z(sample_cpm, ref_cpm=100, ref_mad=0.1):
    return (_norm_signal(sample_cpm) - _norm_signal(ref_cpm)) / (1.4826 * ref_mad)


def test_build_cnv_plot_svg_uses_branch_a_ref_z_and_writes_plot_bins(tmp_path):
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
    assert set(
        [
            "chrom",
            "start",
            "end",
            "genome_pos",
            "z",
            "report_state",
            "branch_a_ref_z",
            "display_ref_z",
            "residual_calibrated_z",
            "z_source",
            "display_z_source",
            "autosomal_neutral_ref_z_median",
            "ref_z_scale",
            "ref_z_scale_source",
        ]
    ).issubset(plot_bins.columns)
    assert set(plot_bins["report_state"]) == {"dup", "del", "neutral"}
    chr1_z = plot_bins.loc[plot_bins["chrom"].eq("chr1"), "z"].tolist()
    assert chr1_z == pytest.approx([0.0, _ref_z(118), _ref_z(119)], abs=1e-6)
    assert set(plot_bins.loc[plot_bins["chrom"].eq("chr1"), "z_source"]) == {
        "branch_a_ref_median_mad_z",
    }
    del_bins = plot_bins.loc[plot_bins["report_state"].eq("del")].sort_values("start")
    assert del_bins["z"].iloc[0] == pytest.approx(_ref_z(84), abs=1e-6)
    assert pd.isna(del_bins["z"].iloc[1])
    assert del_bins["z_source"].iloc[1] == "ref_z_unavailable_masked_or_structure_gap"
    assert 9.9 not in plot_bins["z"].tolist()
    assert "same-direction median display_ref_z" in svg
    assert "Q75" not in svg
    assert "Q25" not in svg
    assert "report z trend" not in svg
    assert "median z" not in svg
    assert "smooth z trend" not in svg
    assert "event ref-z trend" not in svg
    assert "<polyline" not in svg
    assert svg.count('class="report-ref-z-trend"') == 2
    assert "Branch B" not in svg
    assert "report_strong" not in svg
    assert "report_weak" not in svg
    assert "internal" not in svg
    assert "filtered" not in svg
    assert "rejected" not in svg
    assert "mask" not in svg.lower()
    assert "chr1" in svg
    assert 'class="z-bin-scatter"' in svg
    assert re.search(r'<circle[^>]+fill="#64748b"', svg)
    assert re.search(r'<circle[^>]+fill="#facc15"', svg)
    assert re.search(r'<circle[^>]+fill="#2563eb"', svg)
    assert re.search(r'<line[^>]+class="report-ref-z-trend"[^>]+stroke="#dc2626"', svg)
    assert "dup" in svg
    assert "del" in svg


def test_ref_z_event_trend_uses_same_direction_median_without_filling_merge_gap(tmp_path):
    bins = pd.DataFrame(
        {
            "chrom": ["chr1"] * 6,
            "bin_index": list(range(6)),
            "start": [idx * 1_000_000 for idx in range(6)],
            "end": [(idx + 1) * 1_000_000 for idx in range(6)],
            "calibrated_z": [0.0] * 6,
            "normalized_signal": [_norm_signal(value) for value in [130, 130, 100, 100, 100, 100]],
        }
    )
    ref_bins = _ref_bins(
        [("chr1", idx, idx * 1_000_000, (idx + 1) * 1_000_000, 100, 0.05) for idx in range(6)]
    )
    events = pd.DataFrame(
        {
            "event_id": ["merged_evt"],
            "sample_id": ["Y1"],
            "chrom": ["chr1"],
            "start": [0],
            "end": [6_000_000],
            "state": ["gain"],
            "v2_report_layer_class": ["report_event"],
            "v2_report_visibility": ["report_strong_event"],
            "a_zscore": [31.0],
        }
    )
    output_bins = tmp_path / "Y1.plot_bins.tsv"
    output_svg = tmp_path / "Y1.final_cnv.svg"

    build_cnv_plot_svg(
        sample_id="Y1",
        bins_df=bins,
        branch_b_events_df=events,
        a_branch_df=pd.DataFrame(),
        ref_bins_df=ref_bins,
        output_svg=output_svg,
        output_bins_tsv=output_bins,
    )

    plot_bins = pd.read_csv(output_bins, sep="\t")
    assert plot_bins["z"].tolist() == pytest.approx([_ref_z(130, ref_mad=0.05), _ref_z(130, ref_mad=0.05), 0, 0, 0, 0], abs=1e-6)
    assert plot_bins.loc[plot_bins["start"].isin([2_000_000, 3_000_000]), "z"].tolist() == pytest.approx([0.0, 0.0])
    svg = output_svg.read_text(encoding="utf-8")
    assert svg.count('class="report-ref-z-trend"') == 1
    assert "same-direction median display_ref_z" in svg
    assert "same-direction bins=2; valid bins=6" in svg
    assert "Q75 branch_a_ref_z" not in svg
    assert "a_zscore=31.000" in svg
    assert "report z trend" not in svg
    assert "neutral bin" in svg


def test_display_ref_z_recenters_neutral_background_and_support_uses_same_direction_median(tmp_path):
    bins = pd.DataFrame(
        {
            "chrom": ["chr1"] * 8,
            "bin_index": list(range(8)),
            "start": [idx * 1_000_000 for idx in range(8)],
            "end": [(idx + 1) * 1_000_000 for idx in range(8)],
            "calibrated_z": [0.0] * 8,
            "normalized_signal": [_norm_signal(value) for value in [90, 90, 140, 140, 90, 90, 90, 90]],
        }
    )
    ref_bins = _ref_bins(
        [("chr1", idx, idx * 1_000_000, (idx + 1) * 1_000_000, 100, 0.05) for idx in range(8)]
    )
    events = pd.DataFrame(
        {
            "sample_id": ["Y1"],
            "chrom": ["chr1"],
            "start": [2_000_000],
            "end": [6_000_000],
            "state": ["gain"],
            "v2_report_layer_class": ["report_event"],
            "v2_report_visibility": ["report_strong_event"],
            "a_zscore": [28.0],
            "copy_number_estimate": [2.4],
        }
    )
    output_bins = tmp_path / "Y1.plot_bins.tsv"
    output_svg = tmp_path / "Y1.final_cnv.svg"
    support_tsv = tmp_path / "Y1.plot_event_support.tsv"

    build_cnv_plot_svg(
        sample_id="Y1",
        bins_df=bins,
        branch_b_events_df=events,
        a_branch_df=pd.DataFrame(),
        ref_bins_df=ref_bins,
        output_svg=output_svg,
        output_bins_tsv=output_bins,
        output_copy_number_svg=tmp_path / "Y1.final_cnv_cn.svg",
        output_copy_number_bins_tsv=tmp_path / "Y1.plot_bins_cn.tsv",
        output_copy_number_event_support_tsv=support_tsv,
    )

    plot_bins = pd.read_csv(output_bins, sep="\t")
    support = pd.read_csv(support_tsv, sep="\t")
    neutral_bins = plot_bins.loc[plot_bins["report_state"].eq("neutral")]
    event_bins = plot_bins.loc[plot_bins["report_state"].eq("dup")].sort_values("start")
    row = support.iloc[0]

    assert set(["display_ref_z", "autosomal_neutral_ref_z_median", "display_z_source"]).issubset(plot_bins.columns)
    assert neutral_bins["display_ref_z"].median() == pytest.approx(0.0, abs=1e-6)
    assert event_bins["branch_a_ref_z"].iloc[2] == pytest.approx(_ref_z(90, ref_mad=0.05), abs=1e-6)
    assert event_bins["display_ref_z"].iloc[2] == pytest.approx(0.0, abs=1e-6)
    assert event_bins["z"].tolist() == pytest.approx(event_bins["display_ref_z"].tolist(), abs=1e-6)
    assert row["event_layer"] == "autosomal_report"
    assert row["valid_bin_count"] == 4
    assert row["same_direction_ref_z_bin_count"] == 2
    assert row["near_zero_ref_z_bin_count"] == 2
    assert row["same_direction_median_display_ref_z"] == pytest.approx(event_bins["display_ref_z"].iloc[:2].median(), abs=1e-6)
    assert row["median_raw_branch_a_ref_z"] != pytest.approx(row["same_direction_median_display_ref_z"], abs=1e-3)
    svg = output_svg.read_text(encoding="utf-8")
    assert "same-direction median display_ref_z" in svg
    assert "Q75 branch_a_ref_z" not in svg


def test_copy_number_plot_v2_leaves_structure_gaps_blank_and_uses_50mb_ticks(tmp_path):
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

    assert 'class="structure-gap-blank"' not in cn_svg
    assert 'fill="#0f172a" opacity="0.96"' not in cn_svg
    assert "50Mb" in cn_svg
    assert cn_svg.count('class="cn-bin-scatter"') == 59
    assert re.search(r'<circle[^>]+class="cn-bin-scatter"[^>]+r="1\.35"', cn_svg)
    assert cn_svg.count('class="report-cn-trend"') == 2
    assert re.search(r'<line[^>]+class="report-cn-trend"[^>]+stroke="#1d4ed8"', cn_svg)
    assert 'class="chrom-separator"' not in cn_svg
    non_centromere_gap_row = cn_bins.loc[cn_bins["start"].eq(10_000_000)].iloc[0]
    assert non_centromere_gap_row["copy_number_source"] != "structure_gap_blank"
    assert pd.isna(gap_row["copy_number"])
    assert gap_row["copy_number_source"] == "structure_gap_blank"
    assert set(cn_bins.loc[~cn_bins["copy_number_source"].eq("structure_gap_blank"), "copy_number_source"]) == {
        "normalized_signal_ref_median_log2r_autosome_centered"
    }
    assert cn_bins.loc[~cn_bins["copy_number_source"].eq("structure_gap_blank"), "copy_number"].nunique() > 1
    assert "#1d4ed8" in cn_svg
    assert 'class="chrom-background"' in cn_svg


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

    assert 'class="structure-gap-blank"' not in cn_svg
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
    assert event_bins["z"].tolist() == pytest.approx([_ref_z(118), _ref_z(120)], abs=1e-6)


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

    with pytest.raises(ValueError, match="reference bin stability|reference bin medians|log2r|copy_number"):
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
    assert row["same_direction_cn_bin_count"] == 3
    assert "median_residual_calibrated_z" in support.columns
    assert "median_raw_branch_a_ref_z" in support.columns
    assert "same_direction_median_display_ref_z" in support.columns


def test_build_cnv_plot_svg_requires_ref_stability_for_ref_z_without_fallback(tmp_path):
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

    with pytest.raises(ValueError, match="reference bin stability|ref_median|ref_mad"):
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
            ("chrY", 0, 0, 1_000_000, 50, 0.1, "XY"),
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
    assert chr_y["chrY_display_mode"] == "xx_absent_expected_hidden"
    assert pd.isna(chr_y["copy_number"])


def test_xx_chry_is_hidden_from_standard_z_and_cn_scatter(tmp_path):
    bins = pd.DataFrame(
        {
            "chrom": ["chr1", "chrY"],
            "bin_index": [0, 0],
            "start": [0, 0],
            "end": [1_000_000, 1_000_000],
            "calibrated_z": [0.0, 6.0],
            "normalized_signal": [_norm_signal(100), _norm_signal(40)],
        }
    )
    ref_bins = _ref_bins(
        [
            ("chr1", 0, 0, 1_000_000, 100),
            ("chrY", 0, 0, 1_000_000, 5),
        ]
    )
    gender = pd.DataFrame(
        {
            "sample_id": ["XX1"],
            "sex_call": ["XX"],
            "predict_gender": ["F"],
            "bam_y_relative_depth": [0.01],
            "bam_y_to_x_ratio": [0.02],
            "bam_inferred_sex": ["XX"],
        }
    )
    output_bins = tmp_path / "XX1.plot_bins.tsv"
    output_cn_bins = tmp_path / "XX1.plot_bins_cn.tsv"

    build_cnv_plot_svg(
        sample_id="XX1",
        bins_df=bins,
        branch_b_events_df=pd.DataFrame(),
        a_branch_df=pd.DataFrame(),
        gender_df=gender,
        ref_bins_df=ref_bins,
        output_svg=tmp_path / "XX1.final_cnv.svg",
        output_bins_tsv=output_bins,
        output_copy_number_svg=tmp_path / "XX1.final_cnv_cn.svg",
        output_copy_number_bins_tsv=output_cn_bins,
    )

    plot_bins = pd.read_csv(output_bins, sep="\t")
    cn_bins = pd.read_csv(output_cn_bins, sep="\t")
    z_chr_y = plot_bins.loc[plot_bins["chrom"].eq("chrY")].iloc[0]
    cn_chr_y = cn_bins.loc[cn_bins["chrom"].eq("chrY")].iloc[0]

    assert z_chr_y["z_display_mode"] == "xx_absent_expected_hidden"
    assert z_chr_y["chrY_display_mode"] == "xx_absent_expected_hidden"
    assert pd.isna(z_chr_y["z"])
    assert pd.isna(z_chr_y["plot_z"])
    assert np.isfinite(z_chr_y["branch_a_ref_z"])
    assert cn_chr_y["chrY_display_mode"] == "xx_absent_expected_hidden"
    assert cn_chr_y["copy_number_interpretation_status"] == "sex_aware_absent_expected"
    assert pd.isna(cn_chr_y["copy_number"])


def test_xy_chry_zero_mixed_ref_is_unavailable_not_fixed_cn_one(tmp_path):
    bins = pd.DataFrame(
        {
            "chrom": ["chr1", "chrY"],
            "bin_index": [0, 0],
            "start": [0, 0],
            "end": [1_000_000, 1_000_000],
            "calibrated_z": [0.0, 12.0],
            "normalized_signal": [_norm_signal(100), _norm_signal(50)],
        }
    )
    ref_bins = _ref_bins(
        [
            ("chr1", 0, 0, 1_000_000, 100),
            ("chrY", 0, 0, 1_000_000, 0),
        ]
    )
    gender = pd.DataFrame(
        {
            "sample_id": ["XY1"],
            "sex_call": ["XY"],
            "predict_gender": ["M"],
            "bam_y_relative_depth": [0.32],
            "bam_y_to_x_ratio": [0.53],
            "bam_inferred_sex": ["XY"],
        }
    )
    output_svg = tmp_path / "XY1.final_cnv.svg"
    output_bins = tmp_path / "XY1.plot_bins.tsv"
    output_cn_svg = tmp_path / "XY1.final_cnv_cn.svg"
    output_cn_bins = tmp_path / "XY1.plot_bins_cn.tsv"

    build_cnv_plot_svg(
        sample_id="XY1",
        bins_df=bins,
        branch_b_events_df=pd.DataFrame(),
        a_branch_df=pd.DataFrame(),
        gender_df=gender,
        ref_bins_df=ref_bins,
        output_svg=output_svg,
        output_bins_tsv=output_bins,
        output_copy_number_svg=output_cn_svg,
        output_copy_number_bins_tsv=output_cn_bins,
    )

    plot_bins = pd.read_csv(output_bins, sep="\t")
    cn_bins = pd.read_csv(output_cn_bins, sep="\t")
    z_chr_y = plot_bins.loc[plot_bins["chrom"].eq("chrY")].iloc[0]
    cn_chr_y = cn_bins.loc[cn_bins["chrom"].eq("chrY")].iloc[0]
    z_svg = output_svg.read_text(encoding="utf-8")
    cn_svg = output_cn_svg.read_text(encoding="utf-8")

    assert z_chr_y["z_display_mode"] == "sex_chrom_ref_z_unavailable"
    assert z_chr_y["chrY_display_mode"] == "sex_chrom_ref_z_unavailable"
    assert pd.isna(z_chr_y["z"])
    assert pd.isna(z_chr_y["plot_z"])
    assert z_chr_y["z_source"] == "sex_chrom_ref_z_unavailable"
    assert pd.isna(cn_chr_y["copy_number"])
    assert cn_chr_y["copy_number_interpretation_status"] == "sex_chrom_cn_unavailable"
    assert cn_chr_y["copy_number_source"] == "sex_chrom_ref_ratio_not_interpretable"
    assert cn_chr_y["cn_scatter_state_sex_aware"] == "neutral"
    assert 'chry-presence-guide' not in z_svg
    assert 'chry-presence-guide' not in cn_svg
    assert not re.search(r'data-chrom="chrY"[^>]+class="z-bin-scatter"', z_svg)
    assert not re.search(r'data-chrom="chrY"[^>]+class="cn-bin-scatter"', cn_svg)


def test_xy_chry_uses_sex_specific_ref_when_available(tmp_path):
    bins = pd.DataFrame(
        {
            "chrom": ["chr1", "chrY"],
            "bin_index": [0, 0],
            "start": [0, 0],
            "end": [1_000_000, 1_000_000],
            "calibrated_z": [0.0, 0.5],
            "normalized_signal": [_norm_signal(100), _norm_signal(50)],
        }
    )
    ref_bins = _ref_bins(
        [
            ("chr1", 0, 0, 1_000_000, 100, 0.1, "mixed"),
            ("chrY", 0, 0, 1_000_000, 0, 0.1, "mixed"),
            ("chrY", 0, 0, 1_000_000, 50, 0.1, "XY"),
        ]
    )
    gender = pd.DataFrame({"sample_id": ["XY1"], "sex_call": ["XY"], "predict_gender": ["M"]})
    output_bins = tmp_path / "XY1.plot_bins.tsv"
    output_cn_bins = tmp_path / "XY1.plot_bins_cn.tsv"

    build_cnv_plot_svg(
        sample_id="XY1",
        bins_df=bins,
        branch_b_events_df=pd.DataFrame(),
        a_branch_df=pd.DataFrame(),
        gender_df=gender,
        ref_bins_df=ref_bins,
        output_svg=tmp_path / "XY1.final_cnv.svg",
        output_bins_tsv=output_bins,
        output_copy_number_svg=tmp_path / "XY1.final_cnv_cn.svg",
        output_copy_number_bins_tsv=output_cn_bins,
    )

    plot_bins = pd.read_csv(output_bins, sep="\t")
    cn_bins = pd.read_csv(output_cn_bins, sep="\t")
    z_chr_y = plot_bins.loc[plot_bins["chrom"].eq("chrY")].iloc[0]
    cn_chr_y = cn_bins.loc[cn_bins["chrom"].eq("chrY")].iloc[0]

    assert z_chr_y["z_source"] == "branch_a_ref_median_mad_z_sex_specific_XY"
    assert z_chr_y["z_display_mode"] == "branch_a_ref_z_display"
    assert z_chr_y["z"] == pytest.approx(0.0, abs=1e-6)
    assert cn_chr_y["copy_number"] == pytest.approx(1.0, abs=0.02)
    assert cn_chr_y["copy_number_source"] == "normalized_signal_ref_median_log2r_y_haploid_sex_specific_XY"
    assert cn_chr_y["copy_number_interpretation_status"] == "sex_aware_interpretable"


def test_ref_z_plot_hides_out_of_range_points_but_keeps_raw_ref_z(tmp_path):
    bins = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "bin_index": [0, 1],
            "start": [0, 1_000_000],
            "end": [1_000_000, 2_000_000],
            "calibrated_z": [0.0, 0.0],
            "normalized_signal": [_norm_signal(100), _norm_signal(300)],
        }
    )
    ref_bins = _ref_bins(
        [
            ("chr1", 0, 0, 1_000_000, 100, 0.01),
            ("chr1", 1, 1_000_000, 2_000_000, 100, 0.01),
        ]
    )
    output_svg = tmp_path / "Y1.final_cnv.svg"
    output_bins = tmp_path / "Y1.plot_bins.tsv"

    build_cnv_plot_svg(
        sample_id="Y1",
        bins_df=bins,
        branch_b_events_df=pd.DataFrame(),
        a_branch_df=pd.DataFrame(),
        ref_bins_df=ref_bins,
        output_svg=output_svg,
        output_bins_tsv=output_bins,
    )

    plot_bins = pd.read_csv(output_bins, sep="\t")
    svg = output_svg.read_text(encoding="utf-8")
    extreme = plot_bins.loc[plot_bins["start"].eq(1_000_000)].iloc[0]
    assert extreme["branch_a_ref_z"] > 8.0
    assert pd.isna(extreme["plot_z"])
    assert extreme["z_plot_status"] == "out_of_range_hidden"
    assert bool(extreme["z_plot_out_of_range"]) is True
    assert 'data-z-plot-clipped="true"' not in svg
    assert 'data-z-plot-status="out_of_range_hidden"' not in svg


def test_z_plot_is_wide_and_uses_chromosome_background_without_separator_lines(tmp_path):
    bins = pd.DataFrame(
        {
            "chrom": ["chr1", "chr2"],
            "bin_index": [0, 0],
            "start": [0, 0],
            "end": [60_000_000, 60_000_000],
            "calibrated_z": [0.0, 0.0],
            "normalized_signal": [_norm_signal(100), _norm_signal(100)],
        }
    )
    ref_bins = _ref_bins(
        [
            ("chr1", 0, 0, 60_000_000, 100),
            ("chr2", 0, 0, 60_000_000, 100),
        ]
    )
    output_svg = tmp_path / "Y1.final_cnv.svg"

    build_cnv_plot_svg(
        sample_id="Y1",
        bins_df=bins,
        branch_b_events_df=pd.DataFrame(),
        a_branch_df=pd.DataFrame(),
        ref_bins_df=ref_bins,
        output_svg=output_svg,
        output_bins_tsv=tmp_path / "Y1.plot_bins.tsv",
    )

    svg = output_svg.read_text(encoding="utf-8")
    assert 'width="3200"' in svg
    assert 'class="chrom-background"' in svg
    assert 'class="chrom-separator"' not in svg
    assert "50Mb" in svg


def test_z_and_cn_axis_labels_are_bold_and_position_ticks_are_above_chrom_labels(tmp_path):
    bins = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "bin_index": [0, 1],
            "start": [0, 50_000_000],
            "end": [50_000_000, 100_000_000],
            "calibrated_z": [0.0, 0.0],
            "normalized_signal": [_norm_signal(100), _norm_signal(100)],
        }
    )
    ref_bins = _ref_bins(
        [
            ("chr1", 0, 0, 50_000_000, 100),
            ("chr1", 1, 50_000_000, 100_000_000, 100),
        ]
    )
    output_svg = tmp_path / "Y1.final_cnv.svg"
    output_cn_svg = tmp_path / "Y1.final_cnv_cn.svg"

    build_cnv_plot_svg(
        sample_id="Y1",
        bins_df=bins,
        branch_b_events_df=pd.DataFrame(),
        a_branch_df=pd.DataFrame(),
        ref_bins_df=ref_bins,
        output_svg=output_svg,
        output_bins_tsv=tmp_path / "Y1.plot_bins.tsv",
        output_copy_number_svg=output_cn_svg,
        output_copy_number_bins_tsv=tmp_path / "Y1.plot_bins_cn.tsv",
    )

    for svg_text_content in [output_svg.read_text(encoding="utf-8"), output_cn_svg.read_text(encoding="utf-8")]:
        chrom = re.search(
            r'<text[^>]+class="chrom-axis-label"[^>]+y="(?P<y>[0-9.]+)"[^>]+font-size="14"[^>]+font-weight="bold"[^>]*>chr1</text>',
            svg_text_content,
        )
        position = re.search(
            r'<text[^>]+class="chrom-position-label"[^>]+y="(?P<y>[0-9.]+)"[^>]+font-size="11"[^>]+font-weight="bold"[^>]*>50Mb</text>',
            svg_text_content,
        )
        assert chrom is not None
        assert position is not None
        assert float(position.group("y")) < float(chrom.group("y"))
    cn_svg = output_cn_svg.read_text(encoding="utf-8")
    assert re.search(r'<text[^>]+font-size="16"[^>]+font-weight="bold"[^>]*>Copy number</text>', cn_svg)
    assert re.search(r'<text[^>]+font-size="14"[^>]+font-weight="bold"[^>]*>CN=2</text>', cn_svg)


def test_internal_review_event_is_drawn_in_review_plot_without_promoting_report_table(tmp_path):
    bins = pd.DataFrame(
        {
            "chrom": ["chr15", "chr15", "chr15"],
            "bin_index": [0, 1, 2],
            "start": [25_000_000, 26_000_000, 27_000_000],
            "end": [26_000_000, 27_000_000, 28_000_000],
            "calibrated_z": [0.0, 0.0, 0.0],
            "normalized_signal": [_norm_signal(140), _norm_signal(141), _norm_signal(142)],
        }
    )
    ref_bins = _ref_bins(
        [
            ("chr15", 0, 25_000_000, 26_000_000, 100),
            ("chr15", 1, 26_000_000, 27_000_000, 100),
            ("chr15", 2, 27_000_000, 28_000_000, 100),
        ]
    )
    report_events = pd.DataFrame(columns=["sample_id", "chrom", "start", "end", "state"])
    review_events = pd.DataFrame(
        {
            "sample_id": ["H4"],
            "chrom": ["chr15"],
            "start": [25_000_000],
            "end": [28_000_000],
            "state": ["gain"],
            "v2_report_visibility": ["internal_review_event"],
            "v2_report_layer_class": ["internal_review_event"],
            "copy_number_estimate": [2.5],
            "a_zscore": [111.5],
        }
    )
    output_svg = tmp_path / "H4.final_cnv.svg"
    support_tsv = tmp_path / "H4.plot_event_support.tsv"

    build_cnv_plot_svg(
        sample_id="H4",
        bins_df=bins,
        branch_b_events_df=report_events,
        review_events_df=review_events,
        a_branch_df=pd.DataFrame(),
        ref_bins_df=ref_bins,
        output_svg=output_svg,
        output_bins_tsv=tmp_path / "H4.plot_bins.tsv",
        output_copy_number_svg=tmp_path / "H4.final_cnv_cn.svg",
        output_copy_number_bins_tsv=tmp_path / "H4.plot_bins_cn.tsv",
        output_copy_number_event_support_tsv=support_tsv,
    )

    svg = output_svg.read_text(encoding="utf-8")
    support = pd.read_csv(support_tsv, sep="\t")
    assert 'class="internal-review-region"' not in svg
    assert 'plot_event_layer=internal_review_event' not in svg
    assert "internal_review_event" in set(support["event_layer"])
    assert "autosomal_report" not in set(support["event_layer"])
    assert set(support["plot_visibility"]) == {"review_plot_only"}


def test_plot_event_manifest_drives_support_and_hides_z_supported_cn_not_supported_from_final_plot(tmp_path):
    bins = pd.DataFrame(
        {
            "chrom": ["chr1"] * 8,
            "bin_index": list(range(8)),
            "start": [idx * 1_000_000 for idx in range(8)],
            "end": [(idx + 1) * 1_000_000 for idx in range(8)],
            "calibrated_z": [0.0] * 8,
            "normalized_signal": [_norm_signal(value) for value in [100, 100, 110, 111, 110, 100, 100, 100]],
        }
    )
    ref_bins = _ref_bins(
        [("chr1", idx, idx * 1_000_000, (idx + 1) * 1_000_000, 100, 0.01) for idx in range(8)]
    )
    events = pd.DataFrame(
        {
            "event_id": ["G1_chr4_like"],
            "sample_id": ["G1"],
            "chrom": ["chr1"],
            "start": [2_000_000],
            "end": [5_000_000],
            "state": ["gain"],
            "v2_report_layer_class": ["report_event"],
            "v2_report_visibility": ["report_strong_event"],
            "copy_number_estimate": [2.05],
            "a_zscore": [26.0],
        }
    )
    output_svg = tmp_path / "G1.final_cnv.svg"
    output_cn_svg = tmp_path / "G1.final_cnv_cn.svg"
    support_tsv = tmp_path / "G1.plot_event_support.tsv"
    manifest_tsv = tmp_path / "G1.plot_event_manifest.tsv"

    build_cnv_plot_svg(
        sample_id="G1",
        bins_df=bins,
        branch_b_events_df=events,
        a_branch_df=pd.DataFrame(),
        ref_bins_df=ref_bins,
        output_svg=output_svg,
        output_bins_tsv=tmp_path / "G1.plot_bins.tsv",
        output_copy_number_svg=output_cn_svg,
        output_copy_number_bins_tsv=tmp_path / "G1.plot_bins_cn.tsv",
        output_copy_number_event_support_tsv=support_tsv,
        output_event_manifest_tsv=manifest_tsv,
    )

    manifest = pd.read_csv(manifest_tsv, sep="\t")
    support = pd.read_csv(support_tsv, sep="\t")
    svg = output_svg.read_text(encoding="utf-8")
    cn_svg = output_cn_svg.read_text(encoding="utf-8")

    assert len(manifest) == 1
    row = manifest.iloc[0]
    assert row["event_id"] == "G1_chr4_like"
    assert row["candidate_id"] == "G1_chr4_like"
    assert row["event_layer"] == "autosomal_report"
    assert row["plot_support_class"] == "Z_SUPPORTED_CN_NOT_SUPPORTED"
    assert row["plot_layer_class"] == "internal_review_event_candidate"
    assert row["plot_visibility"] == "review_plot_only"
    assert row["plot_visibility_reason"] == "z_supported_cn_not_supported_plot_ablation"
    support_row = support.iloc[0]
    assert support_row["event_id"] == row["event_id"]
    assert support_row["plot_support_class"] == "Z_SUPPORTED_CN_NOT_SUPPORTED"
    assert support_row["plot_visibility"] == "review_plot_only"
    assert support_row["support_interpretation_status"] == "Z_DIRECTION_SUPPORTED"
    assert support_row["cn_direction_consistency_status"] == "CN_DIRECTION_NOT_SUPPORTED"
    assert 'class="report-region"' not in svg
    assert 'class="report-ref-z-trend"' not in svg
    assert 'class="report-cn-region"' not in cn_svg
    assert 'class="report-cn-trend"' not in cn_svg


def test_cn_weak_or_mixed_truth_sensitive_event_stays_visible_in_final_plot(tmp_path):
    bins = pd.DataFrame(
        {
            "chrom": ["chr21"] * 6,
            "bin_index": list(range(6)),
            "start": [idx * 1_000_000 for idx in range(6)],
            "end": [(idx + 1) * 1_000_000 for idx in range(6)],
            "calibrated_z": [0.0] * 6,
            "normalized_signal": [_norm_signal(value) for value in [100, 118, 108, 107, 100, 100]],
        }
    )
    ref_bins = _ref_bins(
        [("chr21", idx, idx * 1_000_000, (idx + 1) * 1_000_000, 100, 0.01) for idx in range(6)]
    )
    events = pd.DataFrame(
        {
            "event_id": ["H6_chr21_gain"],
            "sample_id": ["H6"],
            "chrom": ["chr21"],
            "start": [1_000_000],
            "end": [4_000_000],
            "state": ["gain"],
            "v2_report_layer_class": ["report_event"],
            "v2_report_visibility": ["report_weak_event"],
            "copy_number_estimate": [2.16],
            "a_zscore": [7.11],
        }
    )
    support_tsv = tmp_path / "H6.plot_event_support.tsv"
    manifest_tsv = tmp_path / "H6.plot_event_manifest.tsv"
    output_svg = tmp_path / "H6.final_cnv.svg"

    build_cnv_plot_svg(
        sample_id="H6",
        bins_df=bins,
        branch_b_events_df=events,
        a_branch_df=pd.DataFrame(),
        ref_bins_df=ref_bins,
        output_svg=output_svg,
        output_bins_tsv=tmp_path / "H6.plot_bins.tsv",
        output_copy_number_svg=tmp_path / "H6.final_cnv_cn.svg",
        output_copy_number_bins_tsv=tmp_path / "H6.plot_bins_cn.tsv",
        output_copy_number_event_support_tsv=support_tsv,
        output_event_manifest_tsv=manifest_tsv,
    )

    manifest = pd.read_csv(manifest_tsv, sep="\t")
    support = pd.read_csv(support_tsv, sep="\t")
    svg = output_svg.read_text(encoding="utf-8")

    assert manifest.iloc[0]["plot_support_class"] == "Z_SUPPORTED_CN_WEAK"
    assert manifest.iloc[0]["plot_visibility"] == "final_report_plot"
    assert support.iloc[0]["cn_direction_consistency_status"] == "CN_DIRECTION_WEAK_OR_MIXED"
    assert 'class="report-region"' in svg
    assert 'class="report-ref-z-trend"' in svg


def test_manifest_merges_same_direction_events_across_structure_gap_but_plot_lines_are_split(tmp_path):
    bins = pd.DataFrame(
        {
            "chrom": ["chr1"] * 5,
            "bin_index": list(range(5)),
            "start": [idx * 1_000_000 for idx in range(5)],
            "end": [(idx + 1) * 1_000_000 for idx in range(5)],
            "calibrated_z": [0.0] * 5,
            "normalized_signal": [_norm_signal(value) for value in [130, 131, 100, 132, 133]],
            "is_near_centromere": [0, 0, 1, 0, 0],
        }
    )
    ref_bins = _ref_bins(
        [("chr1", idx, idx * 1_000_000, (idx + 1) * 1_000_000, 100, 0.02) for idx in range(5)]
    )
    events = pd.DataFrame(
        {
            "event_id": ["left_arm", "right_arm"],
            "sample_id": ["S1", "S1"],
            "chrom": ["chr1", "chr1"],
            "start": [0, 3_000_000],
            "end": [2_000_000, 5_000_000],
            "state": ["gain", "gain"],
            "v2_report_layer_class": ["report_event", "report_event"],
            "v2_report_visibility": ["report_strong_event", "report_strong_event"],
            "copy_number_estimate": [2.55, 2.55],
            "a_zscore": [35.0, 34.0],
        }
    )
    manifest_tsv = tmp_path / "S1.plot_event_manifest.tsv"
    output_svg = tmp_path / "S1.final_cnv.svg"

    build_cnv_plot_svg(
        sample_id="S1",
        bins_df=bins,
        branch_b_events_df=events,
        a_branch_df=pd.DataFrame(),
        ref_bins_df=ref_bins,
        output_svg=output_svg,
        output_bins_tsv=tmp_path / "S1.plot_bins.tsv",
        output_copy_number_svg=tmp_path / "S1.final_cnv_cn.svg",
        output_copy_number_bins_tsv=tmp_path / "S1.plot_bins_cn.tsv",
        output_copy_number_event_support_tsv=tmp_path / "S1.plot_event_support.tsv",
        output_event_manifest_tsv=manifest_tsv,
    )

    manifest = pd.read_csv(manifest_tsv, sep="\t")
    svg = output_svg.read_text(encoding="utf-8")

    assert len(manifest) == 1
    row = manifest.iloc[0]
    assert row["start"] == 0
    assert row["end"] == 5_000_000
    assert row["merged_source_event_ids"] == "left_arm;right_arm"
    assert row["plot_visibility"] == "final_report_plot"
    assert svg.count('class="report-region"') == 2
    assert svg.count('class="report-ref-z-trend"') == 2


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
    gender = pd.DataFrame(
        {
            "sample_id": ["XY1"],
            "sex_call": ["XY"],
            "predict_gender": ["M"],
            "bam_y_relative_depth": [0.32],
            "bam_y_to_x_ratio": [0.53],
            "bam_inferred_sex": ["XY"],
        }
    )
    output_cn_bins = tmp_path / "XY1.plot_bins_cn.tsv"
    output_cn_svg = tmp_path / "XY1.final_cnv_cn.svg"

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
    cn_svg = output_cn_svg.read_text(encoding="utf-8")
    chr_y = cn_bins.loc[cn_bins["chrom"].eq("chrY")].iloc[0]
    assert pd.isna(chr_y["copy_number"])
    assert chr_y["cn_scatter_state_sex_aware"] == "neutral"
    assert chr_y["report_state"] == "neutral"
    assert chr_y["copy_number_interpretation_status"] == "sex_chrom_cn_unavailable"
    assert chr_y["copy_number_source"] == "sex_chrom_ref_ratio_not_interpretable"
    assert chr_y["chrY_display_mode"] == "sex_chrom_cn_unavailable"
    assert 'chry-presence-guide' not in cn_svg
    assert 'data-chrom="chrY"' not in cn_svg


def test_copy_number_scatter_uses_valid_xy_chry_bins_when_chrom_median_ref_is_zero(tmp_path):
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
            ("chr1", 0, 0, 1_000_000, 100, 0.1, "mixed"),
            ("chrY", 0, 0, 1_000_000, 0, 0.1, "XY"),
            ("chrY", 1, 1_000_000, 2_000_000, 0, 0.1, "XY"),
            ("chrY", 2, 2_000_000, 3_000_000, 5, 0.1, "XY"),
        ]
    )
    gender = pd.DataFrame(
        {
            "sample_id": ["XY1"],
            "sex_call": ["XY"],
            "predict_gender": ["M"],
            "bam_y_relative_depth": [0.31],
            "bam_y_to_x_ratio": [0.52],
            "bam_inferred_sex": ["XY"],
        }
    )
    output_cn_bins = tmp_path / "XY1.plot_bins_cn.tsv"
    output_cn_svg = tmp_path / "XY1.final_cnv_cn.svg"

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
    chr_y = cn_bins.loc[cn_bins["chrom"].eq("chrY")]
    invalid_bins = chr_y.loc[chr_y["start"].isin([0, 1_000_000])]
    valid_bin = chr_y.loc[chr_y["start"].eq(2_000_000)].iloc[0]

    assert invalid_bins["copy_number"].isna().all()
    assert set(invalid_bins["copy_number_interpretation_status"]) == {"sex_chrom_cn_unavailable"}
    assert set(invalid_bins["copy_number_source"]) == {"sex_chrom_ref_ratio_not_interpretable"}
    assert valid_bin["copy_number"] == pytest.approx(10.0, abs=0.02)
    assert valid_bin["copy_number_interpretation_status"] == "sex_aware_interpretable"
    assert valid_bin["copy_number_source"] == "normalized_signal_ref_median_log2r_y_haploid_sex_specific_XY"
    assert valid_bin["chrY_display_mode"] == "xy_y_ratio_interpretable"
    assert 'data-chrom="chrY"' in output_cn_svg.read_text(encoding="utf-8")


def test_branch_s_cn_trend_skips_uninterpretable_chry_reference(tmp_path):
    bins = pd.DataFrame(
        {
            "chrom": ["chrY", "chrY"],
            "bin_index": [0, 1],
            "start": [0, 1_000_000],
            "end": [1_000_000, 2_000_000],
            "calibrated_z": [0.1, 0.2],
            "normalized_signal": [_norm_signal(50), _norm_signal(50)],
        }
    )
    ref_bins = _ref_bins(
        [
            ("chrY", 0, 0, 1_000_000, 0),
            ("chrY", 1, 1_000_000, 2_000_000, 0),
        ]
    )
    gender = pd.DataFrame({"sample_id": ["XY1"], "sex_call": ["XY"], "predict_gender": ["M"]})
    branch_s_summary = pd.DataFrame(
        {
            "sample_id": ["XY1"],
            "sex_call": ["XY"],
            "sca_candidate_state": ["Y_GAIN"],
            "sca_report_layer_class": ["sca_report_review_event"],
        }
    )
    branch_s_scores = pd.DataFrame({"sample_id": ["XY1"], "sca_state": ["Y_GAIN"], "state_score": [18.0]})
    branch_s_evidence = pd.DataFrame(
        {
            "sample_id": ["XY1"],
            "region_class": ["Y_NONPAR"],
            "chrom": ["chrY"],
            "start": [0],
            "end": [2_000_000],
        }
    )
    output_cn_svg = tmp_path / "XY1.final_cnv_cn.svg"

    build_cnv_plot_svg(
        sample_id="XY1",
        bins_df=bins,
        branch_b_events_df=pd.DataFrame(),
        a_branch_df=pd.DataFrame(),
        gender_df=gender,
        branch_s_summary_df=branch_s_summary,
        branch_s_scores_df=branch_s_scores,
        branch_s_evidence_df=branch_s_evidence,
        ref_bins_df=ref_bins,
        output_svg=tmp_path / "XY1.final_cnv.svg",
        output_bins_tsv=tmp_path / "XY1.plot_bins.tsv",
        output_copy_number_svg=output_cn_svg,
        output_copy_number_bins_tsv=tmp_path / "XY1.plot_bins_cn.tsv",
    )

    assert 'class="branch-s-cn-region"' in output_cn_svg.read_text(encoding="utf-8")
    assert 'class="branch-s-cn-trend"' not in output_cn_svg.read_text(encoding="utf-8")


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
    output_event_support = tmp_path / "G3.plot_event_support.tsv"

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
        output_copy_number_event_support_tsv=output_event_support,
    )

    svg = output_svg.read_text(encoding="utf-8")
    cn_svg = output_cn_svg.read_text(encoding="utf-8")
    cn_bins = pd.read_csv(output_cn_bins, sep="\t")
    support = pd.read_csv(output_event_support, sep="\t")
    assert 'class="branch-s-review-region"' in svg
    assert 'class="branch-s-cn-region"' in cn_svg
    assert 'class="branch-s-cn-trend"' in cn_svg
    assert "event-level score=24.160" in svg
    assert set(cn_bins["event_layer"]) == {"branch_s_review"}
    assert len(support) == 1
    row = support.iloc[0]
    assert row["event_layer"] == "branch_s_review"
    assert row["branch_s_state"] == "X_LOSS"
    assert row["sex_call"] == "XX"
    assert row["sex_chrom_region_class"] == "X_NONPAR"
    assert row["expected_copy_number"] == pytest.approx(2.0)
    assert row["same_direction_ref_z_bin_count"] == 3
    assert row["support_interpretation_status"] == "Z_DIRECTION_SUPPORTED"


def test_branch_s_ref_z_support_is_not_flipped_by_autosomal_centering(tmp_path):
    bins = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1", "chrX", "chrX", "chrX"],
            "bin_index": [0, 1, 2, 0, 1, 2],
            "start": [0, 1_000_000, 2_000_000, 0, 1_000_000, 2_000_000],
            "end": [1_000_000, 2_000_000, 3_000_000, 1_000_000, 2_000_000, 3_000_000],
            "calibrated_z": [-0.1, -0.2, -0.1, -0.1, -0.2, -0.1],
            "normalized_signal": [
                _norm_signal(90),
                _norm_signal(91),
                _norm_signal(89),
                _norm_signal(85),
                _norm_signal(84),
                _norm_signal(86),
            ],
        }
    )
    ref_bins = _ref_bins(
        [
            ("chr1", 0, 0, 1_000_000, 100),
            ("chr1", 1, 1_000_000, 2_000_000, 100),
            ("chr1", 2, 2_000_000, 3_000_000, 100),
            ("chrX", 0, 0, 1_000_000, 100),
            ("chrX", 1, 1_000_000, 2_000_000, 100),
            ("chrX", 2, 2_000_000, 3_000_000, 100),
        ]
    )
    gender = pd.DataFrame({"sample_id": ["G5"], "sex_call": ["XX"], "predict_gender": ["F"]})
    branch_s_summary = pd.DataFrame(
        {
            "sample_id": ["G5"],
            "sex_call": ["XX"],
            "sca_candidate_state": ["X_LOSS"],
            "sca_report_layer_class": ["sca_report_review_event"],
            "sca_report_layer_reason": ["sca_review_weak_report_visible_with_branch_a_support"],
        }
    )
    branch_s_scores = pd.DataFrame(
        {
            "sample_id": ["G5"],
            "sca_state": ["X_LOSS"],
            "state_score": [31.0],
            "state_score_reason": ["branch_a_candidate_zscore_sex_call_compatible_review"],
        }
    )
    branch_s_evidence = pd.DataFrame(
        {
            "sample_id": ["G5"],
            "region_class": ["X_NONPAR"],
            "chrom": ["chrX"],
            "start": [0],
            "end": [3_000_000],
            "bin_count": [3],
            "median_calibrated_z": [-0.1],
        }
    )
    output_event_support = tmp_path / "G5.plot_event_support.tsv"
    output_bins = tmp_path / "G5.plot_bins.tsv"

    build_cnv_plot_svg(
        sample_id="G5",
        bins_df=bins,
        branch_b_events_df=pd.DataFrame(),
        a_branch_df=pd.DataFrame(),
        gender_df=gender,
        branch_s_summary_df=branch_s_summary,
        branch_s_scores_df=branch_s_scores,
        branch_s_evidence_df=branch_s_evidence,
        ref_bins_df=ref_bins,
        output_svg=tmp_path / "G5.final_cnv.svg",
        output_bins_tsv=output_bins,
        output_copy_number_svg=tmp_path / "G5.final_cnv_cn.svg",
        output_copy_number_bins_tsv=tmp_path / "G5.plot_bins_cn.tsv",
        output_copy_number_event_support_tsv=output_event_support,
    )

    plot_bins = pd.read_csv(output_bins, sep="\t")
    autosome_bins = plot_bins[plot_bins["chrom"] == "chr1"]
    x_bins = plot_bins[plot_bins["chrom"] == "chrX"]
    assert autosome_bins["display_ref_z"].median() == pytest.approx(0.0)
    assert (x_bins["display_ref_z"] == x_bins["branch_a_ref_z"]).all()
    assert set(x_bins["display_z_source"]) == {"branch_a_ref_z_sex_chrom_raw_no_autosome_centering"}

    support = pd.read_csv(output_event_support, sep="\t")
    row = support.iloc[0]
    assert row["event_layer"] == "branch_s_review"
    assert row["branch_s_state"] == "X_LOSS"
    assert row["same_direction_ref_z_bin_count"] == 3
    assert row["support_interpretation_status"] == "Z_DIRECTION_SUPPORTED"
