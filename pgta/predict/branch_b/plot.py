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
STATE_COLOR = {"gain": "#dc2626", "loss": "#2563eb", "neutral": "#64748b"}
REPORT_STATE_COLOR = {"dup": "#dc2626", "del": "#2563eb", "neutral": "#64748b"}
NEUTRAL_COLOR = "#64748b"
TREND_COLOR = "#b91c1c"

# Plot idiom reference: cnvpro/cnvseqpipe/CNVcalling.R keeps neutral points subdued,
# highlights final gain/loss segments, and overlays a smoothed trend. The runtime
# implementation stays here to avoid coupling PGT-A reports to cnvpro HDF5/R schema.


def parse_args():
    parser = argparse.ArgumentParser(description="Render a per-sample final CNV gain/loss profile SVG.")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--input-bins", required=True)
    parser.add_argument("--input-events", required=True)
    parser.add_argument("--input-a-branch", default="")
    parser.add_argument("--output-svg", required=True)
    parser.add_argument("--output-bins-tsv", default="")
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


def build_chrom_layout(bins_df):
    layout = {}
    cursor = 0
    gap = 2_000_000
    for chrom in sorted(bins_df["chrom"].unique(), key=chrom_sort_key):
        chrom_df = bins_df[bins_df["chrom"] == chrom]
        chrom_start = int(chrom_df["start"].min())
        chrom_end = int(chrom_df["end"].max())
        span = max(chrom_end - chrom_start, 1)
        layout[chrom] = {"offset": cursor, "start": chrom_start, "end": chrom_end, "span": span}
        cursor += span + gap
    return layout, max(cursor - gap, 1)


def genome_position(chrom, pos, layout):
    item = layout[chrom]
    bounded = min(max(int(pos), item["start"]), item["end"])
    return item["offset"] + (bounded - item["start"])


def scale_x(genome_pos, total_span, left, plot_width):
    return left + (float(genome_pos) / float(total_span)) * plot_width


def scale_y(value, mid_y, half_height):
    clipped = max(min(float(value), 12.0), -12.0)
    return mid_y - (clipped / 12.0) * half_height


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


def smooth_signal_by_chrom(frame):
    smoothed_frames = []
    for chrom, chrom_df in frame.groupby("chrom", sort=False):
        ordered = chrom_df.sort_values("bin_index").copy()
        n_bins = len(ordered)
        if n_bins == 0:
            continue
        window = max(3, int(round(n_bins / 40.0)))
        if window % 2 == 0:
            window += 1
        window = min(window, 101)
        ordered["smooth_signal"] = (
            ordered["plot_signal"]
            .rolling(window=window, center=True, min_periods=1)
            .median()
            .clip(lower=-12.0, upper=12.0)
        )
        smoothed_frames.append(ordered)
    if not smoothed_frames:
        return frame.assign(smooth_signal=frame["plot_signal"])
    return pd.concat(smoothed_frames, ignore_index=True)


def render_smooth_polylines(frame, layout, total_span, left, plot_width, mid_y, half_h):
    chunks = []
    for chrom, chrom_df in frame.groupby("chrom", sort=False):
        if chrom not in layout or chrom_df.empty:
            continue
        points = []
        for row in chrom_df.itertuples(index=False):
            center = int(row.start + ((row.end - row.start) / 2))
            x = scale_x(genome_position(row.chrom, center, layout), total_span, left, plot_width)
            y = scale_y(row.smooth_signal, mid_y, half_h)
            points.append(f"{x:.2f},{y:.2f}")
        if len(points) >= 2:
            chunks.append(
                '<polyline points="'
                + " ".join(points)
                + f'" fill="none" stroke="{TREND_COLOR}" stroke-width="2.0" opacity="0.92">'
                + f"<title>Smoothed signal trend {html.escape(str(chrom))}</title></polyline>"
            )
    return chunks


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


def build_cnv_plot_svg(
    sample_id,
    bins_df,
    branch_b_events_df,
    a_branch_df,
    output_svg,
    output_bins_tsv="",
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
    smooth_bins = smooth_signal_by_chrom(plot_bins)
    point_chunks = []
    for row in plot_bins.itertuples(index=False):
        x = scale_x(genome_position(row.chrom, int(row.start + ((row.end - row.start) / 2)), layout), total_span, left, plot_width)
        y = scale_y(row.plot_signal, mid_y, half_h)
        color = REPORT_STATE_COLOR.get(str(row.report_state), NEUTRAL_COLOR)
        opacity = 0.82 if row.report_state in {"dup", "del"} else 0.62
        point_chunks.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="1.25" fill="{color}" opacity="{opacity:.2f}"/>')
    svg.extend(point_chunks)
    svg.extend(render_smooth_polylines(smooth_bins, layout, total_span, left, plot_width, mid_y, half_h))

    legend_x = left
    legend_y = height - 62
    legend_items = [
        ("dup", REPORT_STATE_COLOR["dup"]),
        ("del", REPORT_STATE_COLOR["del"]),
        ("neutral bin", NEUTRAL_COLOR),
        ("smooth z trend", TREND_COLOR),
    ]
    for idx, (label, color) in enumerate(legend_items):
        x = legend_x + idx * 180
        svg.append(f'<rect x="{x:.2f}" y="{legend_y:.2f}" width="14" height="10" fill="{color}" opacity="0.65"/>')
        svg.append(svg_text(x + 20, legend_y + 10, label, size=12, fill="#334155"))

    svg.append("</svg>")

    output_path = Path(output_svg)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(svg) + "\n", encoding="utf-8")


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
