#!/biosoftware/miniconda/envs/snakemake_env/bin/python
import argparse
import html
import math
from pathlib import Path

import numpy as np
import pandas as pd

from pgta.core.logging import setup_logger


AUTOSOME_ORDER = {f"chr{index}": index for index in range(1, 23)}
CHROM_ORDER = {**AUTOSOME_ORDER, "chrX": 23, "chrY": 24}
STATE_COLOR = {"gain": "#facc15", "loss": "#2563eb", "neutral": "#64748b"}
REPORT_STATE_COLOR = {"dup": "#facc15", "del": "#2563eb", "neutral": "#64748b"}
CN_REPORT_STATE_COLOR = {"dup": "#1d4ed8", "del": "#ef4444", "neutral": "#64748b"}
NEUTRAL_COLOR = "#64748b"
TREND_COLOR = "#dc2626"
CN_CHROM_GAP_BP = 20_000_000
CN_SCATTER_RADIUS = 2.0
HG19_CENTROMERES = {
    "chr1": (121_535_434, 124_535_434),
    "chr2": (92_326_171, 95_326_171),
    "chr3": (90_504_854, 93_504_854),
    "chr4": (49_660_117, 52_660_117),
    "chr5": (46_405_641, 49_405_641),
    "chr6": (58_830_166, 61_830_166),
    "chr7": (58_054_331, 61_054_331),
    "chr8": (43_838_887, 46_838_887),
    "chr9": (47_367_679, 50_367_679),
    "chr10": (39_254_935, 42_254_935),
    "chr11": (51_644_205, 54_644_205),
    "chr12": (34_856_694, 37_856_694),
    "chr13": (16_000_000, 19_000_000),
    "chr14": (16_000_000, 19_000_000),
    "chr15": (17_000_000, 20_000_000),
    "chr16": (35_335_801, 38_335_801),
    "chr17": (22_263_006, 25_263_006),
    "chr18": (15_460_898, 18_460_898),
    "chr19": (24_681_782, 27_681_782),
    "chr20": (26_369_569, 29_369_569),
    "chr21": (11_288_129, 14_288_129),
    "chr22": (13_000_000, 16_000_000),
    "chrX": (58_632_012, 61_632_012),
    "chrY": (10_104_553, 13_104_553),
}

# Plot idiom reference: cnvpro/cnvseqpipe/CNVcalling.R keeps neutral points subdued
# and highlights final gain/loss segments. The runtime implementation stays here
# to avoid coupling PGT-A reports to cnvpro HDF5/R schema.


def parse_args():
    parser = argparse.ArgumentParser(description="Render a per-sample final CNV gain/loss profile SVG.")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--input-bins", required=True)
    parser.add_argument("--input-events", required=True)
    parser.add_argument("--input-a-branch", default="")
    parser.add_argument("--output-svg", required=True)
    parser.add_argument("--output-bins-tsv", default="")
    parser.add_argument("--output-copy-number-svg", default="")
    parser.add_argument("--output-copy-number-bins-tsv", default="")
    parser.add_argument("--max-points", type=int, default=8000)
    parser.add_argument("--log", default="")
    return parser.parse_args()


def read_table(path_value, empty_ok=True):
    if not path_value:
        return pd.DataFrame()
    path = Path(path_value)
    if not path.exists():
        if empty_ok:
            return pd.DataFrame()
        raise FileNotFoundError(path)
    if path.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def normalize_chrom(value):
    token = str(value).strip()
    if not token:
        return ""
    if token.lower().startswith("chr"):
        token = token[3:]
    token = token.upper()
    if token == "23":
        token = "X"
    elif token == "24":
        token = "Y"
    return f"chr{token}"


def chrom_sort_key(chrom):
    normalized = normalize_chrom(chrom)
    return (CHROM_ORDER.get(normalized, 999), normalized)


def choose_signal_column(bins_df):
    if "calibrated_z" not in bins_df.columns:
        raise ValueError("Input bins must contain calibrated_z; robust_z/raw_robust_z/normalized_signal fallback is disabled.")
    return "calibrated_z"


def coerce_bins(bins_df):
    if bins_df.empty:
        raise ValueError("Input bins table is empty.")
    required = {"chrom", "bin_index", "start", "end"}
    missing = required - set(bins_df.columns)
    if missing:
        raise ValueError(f"Input bins table is missing required columns: {sorted(missing)}")
    frame = bins_df.copy()
    frame["chrom"] = frame["chrom"].map(normalize_chrom)
    for column in ("bin_index", "start", "end"):
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    signal_column = choose_signal_column(frame)
    frame["z"] = pd.to_numeric(frame[signal_column], errors="coerce")
    frame = frame.dropna(subset=["bin_index", "start", "end", "z"]).copy()
    frame = frame[np.isfinite(frame["z"].astype(float))].copy()
    frame["plot_signal"] = frame["z"].clip(lower=-12.0, upper=12.0)
    frame["bin_index"] = frame["bin_index"].astype(int)
    frame["start"] = frame["start"].astype(int)
    frame["end"] = frame["end"].astype(int)
    frame = frame[frame["chrom"].isin(CHROM_ORDER)].copy()
    if frame.empty:
        raise ValueError("No supported chromosomes found in input bins table.")
    return frame.sort_values(["chrom"], key=lambda series: series.map(chrom_sort_key)).sort_values(
        ["chrom", "bin_index"]
    )


def normalize_report_state(value):
    text = str(value or "").strip().lower()
    if text in {"gain", "dup", "duplication"}:
        return "dup"
    if text in {"loss", "del", "deletion"}:
        return "del"
    return ""


def coerce_final_events(events_df, sample_id=""):
    if events_df.empty:
        return events_df
    frame = events_df.copy()
    if sample_id and "sample_id" in frame.columns:
        frame = frame[frame["sample_id"].astype(str).eq(str(sample_id))].copy()
    frame["chrom"] = frame["chrom"].map(normalize_chrom)
    frame = frame[frame["chrom"].isin(CHROM_ORDER)].copy()
    for column in ("start", "end", "keep_event", "priority_score"):
        if column in frame.columns:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
    for column, default in {
        "state": "",
        "artifact_status": "unreviewed",
        "technical_confidence": "",
        "artifact_flags": "",
    }.items():
        if column not in frame.columns:
            frame[column] = default
        frame[column] = frame[column].fillna(default).astype(str)
    frame = frame.dropna(subset=["start", "end"]).copy()
    has_v2_report_contract = "v2_report_layer_class" in frame.columns or "v2_report_visibility" in frame.columns
    if has_v2_report_contract:
        report_class = frame.get("v2_report_layer_class", pd.Series("", index=frame.index)).fillna("").astype(str)
        visibility = frame.get("v2_report_visibility", pd.Series("", index=frame.index)).fillna("").astype(str)
        frame = frame[
            report_class.eq("report_event")
            | visibility.isin({"report_strong_event", "report_weak_event", "final_report"})
        ].copy()
    elif "keep_event" in frame.columns:
        frame = frame[pd.to_numeric(frame["keep_event"], errors="coerce").fillna(0).astype(int).eq(1)].copy()
    if not has_v2_report_contract and "artifact_status" in frame.columns:
        frame = frame[~frame["artifact_status"].str.lower().eq("artifact")].copy()
    frame["report_state"] = frame["state"].map(normalize_report_state)
    frame = frame[frame["report_state"].isin({"dup", "del"})].copy()
    for column in ("priority_score", "a_abs_zscore"):
        if column in frame.columns:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
    for column in ("copy_number_estimate", "sex_adjusted_copy_number", "a_ratio"):
        if column in frame.columns:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
    sort_columns = [column for column in ["priority_score", "a_abs_zscore", "end"] if column in frame.columns]
    if sort_columns:
        frame = frame.sort_values(sort_columns, ascending=[False] * len(sort_columns)).copy()
    return frame


def annotate_report_states(bins_df, final_events):
    frame = bins_df.copy()
    frame["report_state"] = "neutral"
    if final_events.empty:
        return frame
    for row in final_events.itertuples(index=False):
        mask = (
            frame["chrom"].astype(str).eq(str(row.chrom))
            & frame["end"].astype(int).gt(int(row.start))
            & frame["start"].astype(int).lt(int(row.end))
            & frame["report_state"].eq("neutral")
        )
        frame.loc[mask, "report_state"] = str(row.report_state)
    return frame


def event_copy_number(row_dict):
    for column in ("copy_number_estimate", "sex_adjusted_copy_number"):
        value = row_dict.get(column, np.nan)
        try:
            numeric = float(value)
        except (TypeError, ValueError):
            numeric = np.nan
        if np.isfinite(numeric):
            return numeric, column
    value = row_dict.get("a_ratio", np.nan)
    try:
        ratio = float(value)
    except (TypeError, ValueError):
        ratio = np.nan
    if np.isfinite(ratio):
        return float(2.0 * (1.0 + ratio)), "a_ratio"
    raise ValueError("Final report event is missing copy_number_estimate, sex_adjusted_copy_number, and a_ratio.")


def annotate_copy_number_bins(bins_df, final_events):
    frame = bins_df.copy()
    frame["copy_number"] = 2.0
    frame["copy_number_source"] = "neutral_diploid_baseline"
    frame["is_structure_gap_blank"] = structural_gap_mask(frame)
    frame.loc[frame["is_structure_gap_blank"], "copy_number"] = np.nan
    frame.loc[frame["is_structure_gap_blank"], "copy_number_source"] = "structure_gap_blank"
    if final_events.empty:
        return frame
    for row in final_events.itertuples(index=False):
        row_dict = row._asdict()
        copy_number, source = event_copy_number(row_dict)
        mask = (
            frame["chrom"].astype(str).eq(str(row_dict["chrom"]))
            & frame["end"].astype(int).gt(int(row_dict["start"]))
            & frame["start"].astype(int).lt(int(row_dict["end"]))
            & ~frame["is_structure_gap_blank"]
        )
        event_bins = frame.loc[mask].copy()
        if event_bins.empty:
            continue
        event_z = pd.to_numeric(event_bins["z"], errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
        median_abs_z = float(event_z.abs().median()) if not event_z.empty else np.nan
        frame.loc[mask, "report_state"] = str(row_dict["report_state"])
        if np.isfinite(median_abs_z) and median_abs_z > 0.0:
            denominator = max(median_abs_z, 0.25)
            amplitude = abs(copy_number - 2.0)
            scaled = 2.0 + (pd.to_numeric(frame.loc[mask, "z"], errors="coerce") * amplitude / denominator)
            frame.loc[mask, "copy_number"] = scaled.astype(float)
            frame.loc[mask, "copy_number_source"] = "event_calibrated_z_scaled_to_cn_proxy"
        else:
            frame.loc[mask, "copy_number"] = copy_number
            frame.loc[mask, "copy_number_source"] = "event_cn_uniform_median_z_uninformative"
    return frame


def boolean_like_series(series):
    text = series.fillna("").astype(str).str.strip().str.lower()
    numeric = pd.to_numeric(series, errors="coerce").fillna(0.0)
    return text.isin({"true", "t", "yes", "y"}) | numeric.gt(0.0)


def structural_gap_mask(bins_df):
    if bins_df.empty:
        return pd.Series(dtype=bool, index=bins_df.index)
    centromere_mask = pd.Series(False, index=bins_df.index)
    centromere_columns_seen = False
    for column in ("is_near_centromere", "near_centromere", "is_centromere", "is_centromere_bin"):
        if column in bins_df.columns:
            centromere_columns_seen = True
            centromere_mask = centromere_mask | boolean_like_series(bins_df[column])
    for column in ("centromere_overlap_fraction", "centromere_fraction", "centromere_bin_fraction"):
        if column in bins_df.columns:
            centromere_columns_seen = True
            overlap = pd.to_numeric(bins_df[column], errors="coerce").fillna(0.0)
            centromere_mask = centromere_mask | overlap.ge(0.5)
    if "nearest_centromere_distance_bp" in bins_df.columns:
        centromere_columns_seen = True
        distance = pd.to_numeric(bins_df["nearest_centromere_distance_bp"], errors="coerce")
        centromere_mask = centromere_mask | distance.le(5_000_000)
    if centromere_columns_seen:
        return centromere_mask

    starts = pd.to_numeric(bins_df.get("start", pd.Series(index=bins_df.index, dtype=float)), errors="coerce")
    ends = pd.to_numeric(bins_df.get("end", pd.Series(index=bins_df.index, dtype=float)), errors="coerce")
    chroms = bins_df.get("chrom", pd.Series(index=bins_df.index, dtype=str)).map(normalize_chrom)
    for chrom, (centromere_start, centromere_end) in HG19_CENTROMERES.items():
        centromere_mask = centromere_mask | (
            chroms.eq(chrom) & starts.lt(centromere_end) & ends.gt(centromere_start)
        )
    return centromere_mask


def build_chrom_layout(bins_df, gap_bp=2_000_000):
    layout = {}
    cursor = 0
    for chrom in sorted(bins_df["chrom"].unique(), key=chrom_sort_key):
        chrom_df = bins_df[bins_df["chrom"] == chrom]
        chrom_start = int(chrom_df["start"].min())
        chrom_end = int(chrom_df["end"].max())
        span = max(chrom_end - chrom_start, 1)
        layout[chrom] = {"offset": cursor, "start": chrom_start, "end": chrom_end, "span": span}
        cursor += span + gap_bp
    return layout, max(cursor - gap_bp, 1)


def genome_position(chrom, pos, layout):
    item = layout[chrom]
    bounded = min(max(int(pos), item["start"]), item["end"])
    return item["offset"] + (bounded - item["start"])


def scale_x(genome_pos, total_span, left, plot_width):
    return left + (float(genome_pos) / float(total_span)) * plot_width


def scale_y(value, mid_y, half_height):
    clipped = max(min(float(value), 12.0), -12.0)
    return mid_y - (clipped / 12.0) * half_height


def scale_copy_number_y(value, mid_y, half_height):
    clipped = max(min(float(value), 4.0), 0.0)
    return mid_y - ((clipped - 2.0) / 2.0) * half_height


def downsample_bins(frame, max_points):
    if max_points <= 0 or len(frame) <= max_points:
        return frame
    stride = int(math.ceil(len(frame) / float(max_points)))
    return frame.iloc[::stride].copy()


def svg_text(x, y, text, size=12, fill="#0f172a", weight="normal", anchor="start"):
    escaped = html.escape(str(text))
    return (
        f'<text x="{x:.2f}" y="{y:.2f}" font-size="{size}" font-family="Arial,sans-serif" '
        f'font-weight="{weight}" text-anchor="{anchor}" fill="{fill}">{escaped}</text>'
    )


def render_event_region(row, layout, total_span, left, plot_width, y, height, color, label):
    chrom = str(row["chrom"])
    if chrom not in layout:
        return ""
    x1 = scale_x(genome_position(chrom, row["start"], layout), total_span, left, plot_width)
    x2 = scale_x(genome_position(chrom, row["end"], layout), total_span, left, plot_width)
    width = max(x2 - x1, 1.5)
    title = html.escape(label)
    return (
        f'<rect x="{x1:.2f}" y="{y:.2f}" width="{width:.2f}" height="{height:.2f}" '
        f'fill="{color}" opacity="0.14" stroke="{color}" stroke-width="1.25">'
        f"<title>{title}</title></rect>"
    )


def render_report_event_trend_lines(bins_df, final_events, layout, total_span, left, plot_width, mid_y, half_h):
    chunks = []
    if final_events.empty:
        return chunks
    for row in final_events.itertuples(index=False):
        row_dict = row._asdict()
        chrom = str(row_dict.get("chrom", ""))
        if chrom not in layout:
            continue
        event_bins = bins_df[
            bins_df["chrom"].astype(str).eq(chrom)
            & bins_df["end"].astype(int).gt(int(row_dict.get("start")))
            & bins_df["start"].astype(int).lt(int(row_dict.get("end")))
        ].copy()
        if event_bins.empty:
            continue
        trend_value = float(event_bins["plot_signal"].median())
        x1 = scale_x(genome_position(chrom, int(row_dict.get("start")), layout), total_span, left, plot_width)
        x2 = scale_x(genome_position(chrom, int(row_dict.get("end")), layout), total_span, left, plot_width)
        if x2 <= x1:
            continue
        y = scale_y(trend_value, mid_y, half_h)
        label = (
            f"report z trend {row_dict.get('report_state')} "
            f"{chrom}:{int(row_dict.get('start'))}-{int(row_dict.get('end'))}; "
            f"median z={trend_value:.3f}"
        )
        chunks.append(
            f'<line class="report-z-trend" x1="{x1:.2f}" y1="{y:.2f}" '
            f'x2="{x2:.2f}" y2="{y:.2f}" stroke="{TREND_COLOR}" '
            f'stroke-width="2.4" opacity="0.94" stroke-linecap="round">'
            f"<title>{html.escape(label)}</title></line>"
        )
    return chunks


def render_report_event_cn_trend_lines(bins_df, final_events, layout, total_span, left, plot_width, mid_y, half_h):
    chunks = []
    if final_events.empty:
        return chunks
    for row in final_events.itertuples(index=False):
        row_dict = row._asdict()
        chrom = str(row_dict.get("chrom", ""))
        if chrom not in layout:
            continue
        trend_value, _source = event_copy_number(row_dict)
        event_bins = bins_df[
            bins_df["chrom"].astype(str).eq(chrom)
            & bins_df["end"].astype(int).gt(int(row_dict.get("start")))
            & bins_df["start"].astype(int).lt(int(row_dict.get("end")))
            & ~bins_df.get("is_structure_gap_blank", pd.Series(False, index=bins_df.index)).astype(bool)
        ].copy()
        if event_bins.empty:
            continue
        event_bins = event_bins.sort_values(["chrom", "start"]).copy()
        current = []
        previous_end = None
        for bin_row in event_bins.itertuples(index=False):
            if previous_end is not None and int(bin_row.start) > int(previous_end):
                chunks.extend(render_cn_trend_chunk(current, row_dict, trend_value, layout, total_span, left, plot_width, mid_y, half_h))
                current = []
            current.append(bin_row)
            previous_end = int(bin_row.end)
        chunks.extend(render_cn_trend_chunk(current, row_dict, trend_value, layout, total_span, left, plot_width, mid_y, half_h))
    return chunks


def render_cn_trend_chunk(chunk_rows, row_dict, trend_value, layout, total_span, left, plot_width, mid_y, half_h):
    if not chunk_rows:
        return []
    chrom = str(row_dict.get("chrom", ""))
    x1 = scale_x(genome_position(chrom, int(chunk_rows[0].start), layout), total_span, left, plot_width)
    x2 = scale_x(genome_position(chrom, int(chunk_rows[-1].end), layout), total_span, left, plot_width)
    if x2 <= x1:
        return []
    y = scale_copy_number_y(trend_value, mid_y, half_h)
    label = (
        f"report CN trend {row_dict.get('report_state')} "
        f"{chrom}:{int(chunk_rows[0].start)}-{int(chunk_rows[-1].end)}; "
        f"event CN={trend_value:.3f}"
    )
    state = str(row_dict.get("report_state", "")).lower()
    trend_color = CN_REPORT_STATE_COLOR.get(state, TREND_COLOR)
    return [
        f'<line class="report-cn-trend" x1="{x1:.2f}" y1="{y:.2f}" '
        f'x2="{x2:.2f}" y2="{y:.2f}" stroke="{trend_color}" '
        f'stroke-width="2.4" opacity="0.94" stroke-linecap="round">'
        f"<title>{html.escape(label)}</title></line>"
    ]


def write_plot_bins_tsv(path_value, bins_df, layout):
    if not path_value:
        return
    output = bins_df.copy()
    genome_positions = []
    for row in output.itertuples(index=False):
        center = int(row.start + ((row.end - row.start) / 2))
        genome_positions.append(int(genome_position(row.chrom, center, layout)))
    output["genome_pos"] = genome_positions
    output = output[["chrom", "start", "end", "genome_pos", "z", "report_state"]].copy()
    path = Path(path_value)
    path.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(path, sep="\t", index=False)


def write_copy_number_bins_tsv(path_value, bins_df, layout):
    if not path_value:
        return
    output = bins_df.copy()
    genome_positions = []
    for row in output.itertuples(index=False):
        center = int(row.start + ((row.end - row.start) / 2))
        genome_positions.append(int(genome_position(row.chrom, center, layout)))
    output["genome_pos"] = genome_positions
    output = output[["chrom", "start", "end", "genome_pos", "z", "copy_number", "report_state", "copy_number_source"]].copy()
    path = Path(path_value)
    path.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(path, sep="\t", index=False)


def build_copy_number_plot_svg(sample_id, bins, final_events, layout, total_span, output_svg, output_bins_tsv="", max_points=8000):
    cn_bins = annotate_copy_number_bins(bins, final_events)
    cn_layout, cn_total_span = build_chrom_layout(cn_bins, gap_bp=CN_CHROM_GAP_BP)
    write_copy_number_bins_tsv(output_bins_tsv, cn_bins, cn_layout)

    width = 2560
    height = 620
    left = 82
    right = 42
    top = 74
    plot_width = width - left - right
    signal_top = top + 72
    signal_height = 350
    mid_y = signal_top + signal_height / 2.0
    half_h = signal_height / 2.0

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        svg_text(24, 34, f"CNV final copy-number profile - {sample_id}", size=24, weight="bold"),
        svg_text(24, 58, "Final reported dup/del regions over event-level copy-number proxy", size=13, fill="#475569"),
        svg_text(20, signal_top - 18, "Copy number", size=12, fill="#334155"),
        f'<rect x="{left:.2f}" y="{signal_top:.2f}" width="{plot_width:.2f}" height="{signal_height:.2f}" fill="#ffffff"/>',
    ]

    for idx, chrom in enumerate(sorted(cn_layout, key=chrom_sort_key)):
        item = cn_layout[chrom]
        x1 = scale_x(item["offset"], cn_total_span, left, plot_width)
        x2 = scale_x(item["offset"] + item["span"], cn_total_span, left, plot_width)
        svg.append(
            f'<line class="chrom-separator" x1="{x1:.2f}" y1="{signal_top:.2f}" '
            f'x2="{x1:.2f}" y2="{signal_top + signal_height:.2f}" '
            f'stroke="#94a3b8" stroke-width="1.1" opacity="0.95"/>'
        )
        svg.append(svg_text((x1 + x2) / 2.0, signal_top + signal_height + 18, chrom, size=10, fill="#334155", anchor="middle"))
        tick = ((int(item["start"]) // 50_000_000) + 1) * 50_000_000
        while tick < int(item["end"]):
            tick_x = scale_x(genome_position(chrom, tick, cn_layout), cn_total_span, left, plot_width)
            svg.append(
                f'<line class="chrom-50mb-tick" x1="{tick_x:.2f}" y1="{signal_top + signal_height:.2f}" '
                f'x2="{tick_x:.2f}" y2="{signal_top + signal_height + 6:.2f}" stroke="#64748b" stroke-width="0.7"/>'
            )
            svg.append(svg_text(tick_x, signal_top + signal_height + 31, f"{int(tick / 1_000_000)}Mb", size=8, fill="#64748b", anchor="middle"))
            tick += 50_000_000

    gap_bins = cn_bins[cn_bins["is_structure_gap_blank"].astype(bool)].copy()
    for row in gap_bins.itertuples(index=False):
        if str(row.chrom) not in cn_layout:
            continue
        x1 = scale_x(genome_position(row.chrom, int(row.start), cn_layout), cn_total_span, left, plot_width)
        x2 = scale_x(genome_position(row.chrom, int(row.end), cn_layout), cn_total_span, left, plot_width)
        svg.append(
            f'<rect class="structure-gap-blank" x="{x1:.2f}" y="{signal_top:.2f}" '
            f'width="{max(x2 - x1, 1):.2f}" height="{signal_height:.2f}" fill="#cbd5e1" opacity="0.82"/>'
        )

    for cn in (1, 2, 3):
        y = scale_copy_number_y(cn, mid_y, half_h)
        color = "#94a3b8" if cn == 2 else "#cbd5e1"
        svg.append(f'<line x1="{left}" y1="{y:.2f}" x2="{left + plot_width}" y2="{y:.2f}" stroke="{color}" stroke-width="1" stroke-dasharray="4,4"/>')
        svg.append(svg_text(18, y + 4, f"CN={cn}", size=11, fill="#64748b"))

    for row in final_events.itertuples(index=False):
        row_dict = row._asdict()
        state = str(row_dict.get("report_state", "")).lower()
        color = CN_REPORT_STATE_COLOR.get(state, NEUTRAL_COLOR)
        label = f"{state} {row_dict.get('chrom')}:{int(row_dict.get('start'))}-{int(row_dict.get('end'))}"
        svg.append(render_event_region(row_dict, cn_layout, cn_total_span, left, plot_width, signal_top, signal_height, color, label))

    plot_bins = downsample_bins(cn_bins, max_points=max_points)
    for row in plot_bins.itertuples(index=False):
        if getattr(row, "is_structure_gap_blank", False) or not np.isfinite(float(row.copy_number)):
            continue
        x = scale_x(genome_position(row.chrom, int(row.start + ((row.end - row.start) / 2)), cn_layout), cn_total_span, left, plot_width)
        y = scale_copy_number_y(row.copy_number, mid_y, half_h)
        color = CN_REPORT_STATE_COLOR.get(str(row.report_state), NEUTRAL_COLOR)
        opacity = 0.82 if row.report_state in {"dup", "del"} else 0.62
        svg.append(
            f'<circle class="cn-bin-scatter" cx="{x:.2f}" cy="{y:.2f}" '
            f'r="{CN_SCATTER_RADIUS:.2f}" fill="{color}" opacity="{opacity:.2f}" '
            f'stroke="#ffffff" stroke-width="0.35"/>'
        )
    svg.extend(render_report_event_cn_trend_lines(cn_bins, final_events, cn_layout, cn_total_span, left, plot_width, mid_y, half_h))

    legend_x = left
    legend_y = height - 62
    legend_items = [
        ("dup", CN_REPORT_STATE_COLOR["dup"]),
        ("del", CN_REPORT_STATE_COLOR["del"]),
    ]
    for idx, (label, color) in enumerate(legend_items):
        x = legend_x + idx * 180
        svg.append(f'<rect x="{x:.2f}" y="{legend_y:.2f}" width="14" height="10" fill="{color}" opacity="0.65"/>')
        svg.append(svg_text(x + 20, legend_y + 10, label, size=12, fill="#334155"))

    svg.append("</svg>")

    output_path = Path(output_svg)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(svg) + "\n", encoding="utf-8")


def build_cnv_plot_svg(
    sample_id,
    bins_df,
    branch_b_events_df,
    a_branch_df,
    output_svg,
    output_bins_tsv="",
    output_copy_number_svg="",
    output_copy_number_bins_tsv="",
    max_points=8000,
):
    bins = coerce_bins(bins_df)
    final_events = coerce_final_events(branch_b_events_df, sample_id=sample_id)
    layout, total_span = build_chrom_layout(bins)
    bins = annotate_report_states(bins, final_events)
    write_plot_bins_tsv(output_bins_tsv, bins, layout)

    width = 1280
    height = 620
    left = 70
    right = 30
    top = 74
    plot_width = width - left - right
    signal_top = top + 72
    signal_height = 350
    mid_y = signal_top + signal_height / 2.0
    half_h = signal_height / 2.0

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        svg_text(24, 34, f"CNV final profile - {sample_id}", size=24, weight="bold"),
        svg_text(24, 58, "Final reported dup/del regions over genome-wide calibrated z signal", size=13, fill="#475569"),
    ]

    # Chromosome background and labels.
    for idx, chrom in enumerate(sorted(layout, key=chrom_sort_key)):
        item = layout[chrom]
        x1 = scale_x(item["offset"], total_span, left, plot_width)
        x2 = scale_x(item["offset"] + item["span"], total_span, left, plot_width)
        fill = "#f8fafc" if idx % 2 == 0 else "#eef2f7"
        svg.append(
            f'<rect x="{x1:.2f}" y="{signal_top:.2f}" width="{max(x2 - x1, 1):.2f}" height="{signal_height:.2f}" fill="{fill}"/>'
        )
        svg.append(f'<line x1="{x1:.2f}" y1="{signal_top:.2f}" x2="{x1:.2f}" y2="{signal_top + signal_height:.2f}" stroke="#cbd5e1" stroke-width="0.5"/>')
        svg.append(svg_text((x1 + x2) / 2.0, signal_top + signal_height + 18, chrom, size=10, fill="#334155", anchor="middle"))

    for z in (-6, 0, 6):
        y = scale_y(z, mid_y, half_h)
        color = "#94a3b8" if z == 0 else "#cbd5e1"
        svg.append(f'<line x1="{left}" y1="{y:.2f}" x2="{left + plot_width}" y2="{y:.2f}" stroke="{color}" stroke-width="1" stroke-dasharray="4,4"/>')
        svg.append(svg_text(20, y + 4, f"z={z}", size=11, fill="#64748b"))

    for row in final_events.itertuples(index=False):
        row_dict = row._asdict()
        state = str(row_dict.get("report_state", "")).lower()
        color = STATE_COLOR.get(state, NEUTRAL_COLOR)
        color = REPORT_STATE_COLOR.get(state, NEUTRAL_COLOR)
        label = f"{state} {row_dict.get('chrom')}:{int(row_dict.get('start'))}-{int(row_dict.get('end'))}"
        svg.append(render_event_region(row_dict, layout, total_span, left, plot_width, signal_top, signal_height, color, label))

    plot_bins = downsample_bins(bins, max_points=max_points)
    point_chunks = []
    for row in plot_bins.itertuples(index=False):
        x = scale_x(genome_position(row.chrom, int(row.start + ((row.end - row.start) / 2)), layout), total_span, left, plot_width)
        y = scale_y(row.plot_signal, mid_y, half_h)
        color = REPORT_STATE_COLOR.get(str(row.report_state), NEUTRAL_COLOR)
        opacity = 0.82 if row.report_state in {"dup", "del"} else 0.62
        point_chunks.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="1.25" fill="{color}" opacity="{opacity:.2f}"/>')
    svg.extend(point_chunks)
    svg.extend(render_report_event_trend_lines(bins, final_events, layout, total_span, left, plot_width, mid_y, half_h))

    legend_x = left
    legend_y = height - 62
    legend_items = [
        ("dup", REPORT_STATE_COLOR["dup"]),
        ("del", REPORT_STATE_COLOR["del"]),
        ("neutral bin", NEUTRAL_COLOR),
        ("report z trend", TREND_COLOR),
    ]
    for idx, (label, color) in enumerate(legend_items):
        x = legend_x + idx * 180
        svg.append(f'<rect x="{x:.2f}" y="{legend_y:.2f}" width="14" height="10" fill="{color}" opacity="0.65"/>')
        svg.append(svg_text(x + 20, legend_y + 10, label, size=12, fill="#334155"))

    svg.append("</svg>")

    output_path = Path(output_svg)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(svg) + "\n", encoding="utf-8")
    if output_copy_number_svg or output_copy_number_bins_tsv:
        build_copy_number_plot_svg(
            sample_id=sample_id,
            bins=bins,
            final_events=final_events,
            layout=layout,
            total_span=total_span,
            output_svg=output_copy_number_svg,
            output_bins_tsv=output_copy_number_bins_tsv,
            max_points=max_points,
        )


def main():
    args = parse_args()
    logger = setup_logger("cnv_branch_b_plot", args.log or None)
    bins_df = read_table(args.input_bins, empty_ok=False)
    events_df = read_table(args.input_events, empty_ok=True)
    a_branch_df = read_table(args.input_a_branch, empty_ok=True)
    build_cnv_plot_svg(
        sample_id=args.sample_id,
        bins_df=bins_df,
        branch_b_events_df=events_df,
        a_branch_df=a_branch_df,
        output_svg=args.output_svg,
        output_bins_tsv=args.output_bins_tsv,
        output_copy_number_svg=args.output_copy_number_svg,
        output_copy_number_bins_tsv=args.output_copy_number_bins_tsv,
        max_points=args.max_points,
    )
    logger.info(
        "CNV plot written: sample=%s bins=%d branch_b_events=%d a_branch_events=%d output=%s",
        args.sample_id,
        len(bins_df),
        len(events_df),
        len(a_branch_df),
        args.output_svg,
    )


if __name__ == "__main__":
    main()
