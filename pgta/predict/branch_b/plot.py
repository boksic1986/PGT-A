#!/biosoftware/miniconda/envs/snakemake_env/bin/python
import argparse
import html
import json
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
CN_SCATTER_RADIUS = 1.35
CN_NEUTRAL_LOWER = 1.7
CN_NEUTRAL_UPPER = 2.3
CN_HAPLOID_NEUTRAL_LOWER = 0.85
CN_HAPLOID_NEUTRAL_UPPER = 1.15
SEX_CHROM_REF_MIN_CPM = 1.0
CN_COPY_NUMBER_MIN = 0.0
CN_COPY_NUMBER_MAX = 4.0
CN_LOG2R_PSEUDOCOUNT = 1e-3
CN_BIN_COPY_NUMBER_SOURCE = "normalized_signal_ref_median_log2r_autosome_centered"
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
    parser.add_argument("--input-ref-bins", default="")
    parser.add_argument("--gender-tsv", default="")
    parser.add_argument("--branch-s-summary", default="")
    parser.add_argument("--branch-s-scores", default="")
    parser.add_argument("--branch-s-evidence", default="")
    parser.add_argument("--output-svg", required=True)
    parser.add_argument("--output-bins-tsv", default="")
    parser.add_argument("--output-copy-number-svg", default="")
    parser.add_argument("--output-copy-number-bins-tsv", default="")
    parser.add_argument("--output-copy-number-event-support-tsv", default="")
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
    if path.suffix.lower() == ".json":
        payload = json.loads(path.read_text(encoding="utf-8"))
        if isinstance(payload, list):
            return pd.DataFrame(payload)
        if isinstance(payload, dict):
            return pd.DataFrame([payload])
        raise ValueError(f"Unsupported JSON payload for table input: {path}")
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


def coerce_bins(bins_df):
    if bins_df.empty:
        raise ValueError("Input bins table is empty.")
    required = {"chrom", "bin_index", "start", "end", "normalized_signal"}
    missing = required - set(bins_df.columns)
    if missing:
        raise ValueError(f"Input bins table is missing required columns: {sorted(missing)}")
    frame = bins_df.copy()
    frame["chrom"] = frame["chrom"].map(normalize_chrom)
    for column in ("bin_index", "start", "end"):
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame["normalized_signal"] = pd.to_numeric(frame["normalized_signal"], errors="coerce")
    if "calibrated_z" in frame.columns:
        frame["residual_calibrated_z"] = pd.to_numeric(frame["calibrated_z"], errors="coerce")
    else:
        frame["residual_calibrated_z"] = np.nan
    frame = frame.dropna(subset=["bin_index", "start", "end", "normalized_signal"]).copy()
    frame = frame[np.isfinite(frame["normalized_signal"].astype(float))].copy()
    frame["bin_index"] = frame["bin_index"].astype(int)
    frame["start"] = frame["start"].astype(int)
    frame["end"] = frame["end"].astype(int)
    frame = frame[frame["chrom"].isin(CHROM_ORDER)].copy()
    if frame.empty:
        raise ValueError("No supported chromosomes found in input bins table.")
    return frame.sort_values(["chrom"], key=lambda series: series.map(chrom_sort_key)).sort_values(
        ["chrom", "bin_index"]
    )


def validate_ref_bins_for_ref_z(ref_bins_df):
    if ref_bins_df is None or ref_bins_df.empty:
        raise ValueError("Input reference bin stability is required for branch_a_ref_z.")
    required = {"chrom", "bin_index", "ref_median", "ref_mad"}
    missing = required - set(ref_bins_df.columns)
    if missing:
        raise ValueError(f"Input reference bin stability is missing required columns: {sorted(missing)}")


def hard_mask_series(frame):
    if "mask_label" not in frame.columns:
        return pd.Series(False, index=frame.index)
    label = frame["mask_label"].fillna("").astype(str).str.strip().str.lower()
    return label.isin({"hard", "hard_mask", "hard_exclude", "exclude", "excluded"})


def annotate_branch_a_ref_z_bins(bins_df, ref_bins_df):
    validate_ref_bins_for_ref_z(ref_bins_df)
    frame = bins_df.copy()
    ref = ref_bins_df.copy()
    ref["chrom"] = ref["chrom"].map(normalize_chrom)
    ref["bin_index"] = pd.to_numeric(ref["bin_index"], errors="coerce")
    ref["ref_median"] = pd.to_numeric(ref["ref_median"], errors="coerce")
    ref["ref_mad"] = pd.to_numeric(ref["ref_mad"], errors="coerce")
    ref = ref.dropna(subset=["chrom", "bin_index"]).copy()
    ref["bin_index"] = ref["bin_index"].astype(int)
    ref = ref.sort_values(["chrom", "bin_index"]).drop_duplicates(["chrom", "bin_index"], keep="first")
    frame = frame.merge(ref[["chrom", "bin_index", "ref_median", "ref_mad"]], on=["chrom", "bin_index"], how="left")

    raw_scale = 1.4826 * pd.to_numeric(frame["ref_mad"], errors="coerce")
    structure_or_hard = structural_gap_mask(frame) | hard_mask_series(frame)
    autosome_valid_scale = (
        frame["chrom"].isin(AUTOSOME_ORDER)
        & ~structure_or_hard
        & np.isfinite(raw_scale)
        & raw_scale.gt(0.0)
    )
    scale_pool = raw_scale.loc[autosome_valid_scale].dropna()
    scale_source_label = "autosomal_non_gap_mad_x1.4826_floor_p10"
    if scale_pool.empty:
        all_valid_scale = ~structure_or_hard & np.isfinite(raw_scale) & raw_scale.gt(0.0)
        scale_pool = raw_scale.loc[all_valid_scale].dropna()
        scale_source_label = "all_non_gap_mad_x1.4826_floor_p10"
    if scale_pool.empty:
        raise ValueError("No valid autosomal reference MAD scale is available for branch_a_ref_z.")
    min_ref_scale = float(scale_pool.quantile(0.10))
    if not np.isfinite(min_ref_scale) or min_ref_scale <= 0.0:
        raise ValueError("Invalid reference scale floor for branch_a_ref_z.")

    ref_scale = raw_scale.copy()
    ref_scale = ref_scale.where(ref_scale.ge(min_ref_scale), min_ref_scale)
    valid = (
        ~structure_or_hard
        & np.isfinite(frame["normalized_signal"])
        & np.isfinite(frame["ref_median"])
        & np.isfinite(ref_scale)
        & ref_scale.gt(0.0)
    )
    frame["ref_z_scale"] = np.where(valid, ref_scale, np.nan)
    frame["ref_z_scale_source"] = np.where(
        valid,
        f"{scale_source_label}={min_ref_scale:.6g}",
        "ref_z_scale_unavailable",
    )
    frame["branch_a_ref_z"] = np.nan
    frame.loc[valid, "branch_a_ref_z"] = (
        frame.loc[valid, "normalized_signal"] - frame.loc[valid, "ref_median"]
    ) / ref_scale.loc[valid]
    frame["z"] = frame["branch_a_ref_z"]
    frame["plot_signal"] = pd.to_numeric(frame["z"], errors="coerce").clip(lower=-12.0, upper=12.0)
    frame["z_source"] = "ref_z_unavailable_invalid_ref"
    frame.loc[valid, "z_source"] = "branch_a_ref_median_mad_z"
    frame.loc[structure_or_hard, "z_source"] = "ref_z_unavailable_masked_or_structure_gap"
    return frame


def normalize_report_state(value):
    text = str(value or "").strip().lower()
    if text in {"gain", "dup", "duplication"}:
        return "dup"
    if text in {"loss", "del", "deletion"}:
        return "del"
    return ""


def sample_sex_call(gender_df, branch_s_summary_df=None, sample_id=""):
    for frame in (gender_df, branch_s_summary_df):
        if frame is None or frame.empty or "sex_call" not in frame.columns:
            continue
        rows = frame
        if sample_id and "sample_id" in frame.columns:
            rows = frame[frame["sample_id"].astype(str).eq(str(sample_id))]
        if rows.empty:
            continue
        value = str(rows.iloc[0].get("sex_call", "") or "").strip().upper()
        if value in {"XX", "XY"}:
            return value
    return ""


def sex_chrom_region_class(chrom, is_par=False):
    normalized = normalize_chrom(chrom)
    if normalized == "chrX":
        return "X_PAR" if bool(is_par) else "X_NONPAR"
    if normalized == "chrY":
        return "Y_PAR" if bool(is_par) else "Y_NONPAR"
    return "AUTOSOME"


def expected_copy_number_for_bin(chrom, sex_call, is_par=False):
    normalized = normalize_chrom(chrom)
    sex = str(sex_call or "").strip().upper()
    if normalized in AUTOSOME_ORDER:
        return 2.0
    if bool(is_par):
        return np.nan
    if normalized == "chrX":
        if sex == "XX":
            return 2.0
        if sex == "XY":
            return 1.0
    if normalized == "chrY":
        if sex == "XX":
            return 0.0
        if sex == "XY":
            return 1.0
    return np.nan


def branch_s_state_to_report_state(state):
    text = str(state or "").strip().upper()
    if text.endswith("_GAIN"):
        return "dup"
    if text.endswith("_LOSS"):
        return "del"
    return ""


def branch_s_axis_from_state(state):
    text = str(state or "").strip().upper()
    if text.startswith("X_"):
        return "X"
    if text.startswith("Y_"):
        return "Y"
    return ""


def coerce_branch_s_review_events(summary_df=None, scores_df=None, evidence_df=None, sample_id=""):
    if summary_df is None or summary_df.empty:
        return pd.DataFrame()
    frame = summary_df.copy()
    if sample_id and "sample_id" in frame.columns:
        frame = frame[frame["sample_id"].astype(str).eq(str(sample_id))].copy()
    if frame.empty or "sca_report_layer_class" not in frame.columns:
        return pd.DataFrame()
    frame = frame[frame["sca_report_layer_class"].astype(str).eq("sca_report_review_event")].copy()
    if frame.empty:
        return pd.DataFrame()
    scores = scores_df.copy() if scores_df is not None else pd.DataFrame()
    if not scores.empty and sample_id and "sample_id" in scores.columns:
        scores = scores[scores["sample_id"].astype(str).eq(str(sample_id))].copy()
    evidence = evidence_df.copy() if evidence_df is not None else pd.DataFrame()
    if not evidence.empty and sample_id and "sample_id" in evidence.columns:
        evidence = evidence[evidence["sample_id"].astype(str).eq(str(sample_id))].copy()
    rows = []
    for _, summary in frame.iterrows():
        state = str(summary.get("sca_candidate_state", "") or "")
        report_state = branch_s_state_to_report_state(state)
        axis = branch_s_axis_from_state(state)
        if not report_state or not axis:
            continue
        region = pd.DataFrame()
        if not evidence.empty and "region_class" in evidence.columns:
            region = evidence[evidence["region_class"].astype(str).eq(f"{axis}_NONPAR")].copy()
        if region.empty:
            continue
        start = int(pd.to_numeric(region["start"], errors="coerce").min())
        end = int(pd.to_numeric(region["end"], errors="coerce").max())
        score = np.nan
        reason = ""
        if not scores.empty and "sca_state" in scores.columns:
            matched = scores[scores["sca_state"].astype(str).eq(state)].copy()
            if not matched.empty:
                score = pd.to_numeric(matched.iloc[0].get("state_score", np.nan), errors="coerce")
                reason = str(matched.iloc[0].get("state_score_reason", "") or "")
        rows.append(
            {
                "sample_id": str(summary.get("sample_id", sample_id) or sample_id),
                "chrom": f"chr{axis}",
                "start": start,
                "end": end,
                "state": "gain" if report_state == "dup" else "loss",
                "report_state": report_state,
                "event_layer": "branch_s_review",
                "branch_s_state": state,
                "branch_s_score": float(score) if np.isfinite(score) else np.nan,
                "branch_s_score_reason": reason,
                "branch_s_report_layer_reason": str(summary.get("sca_report_layer_reason", "") or ""),
            }
        )
    return pd.DataFrame(rows)


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


def coerce_ref_bins(ref_bins_df):
    if ref_bins_df is None or ref_bins_df.empty:
        return pd.DataFrame()
    required = {"chrom", "bin_index", "ref_median"}
    missing = required - set(ref_bins_df.columns)
    if missing:
        raise ValueError(f"Reference bin table is missing required columns: {sorted(missing)}")
    frame = ref_bins_df.copy()
    frame["chrom"] = frame["chrom"].map(normalize_chrom)
    frame["bin_index"] = pd.to_numeric(frame["bin_index"], errors="coerce")
    frame["ref_median"] = pd.to_numeric(frame["ref_median"], errors="coerce")
    frame = frame.dropna(subset=["chrom", "bin_index", "ref_median"]).copy()
    frame["bin_index"] = frame["bin_index"].astype(int)
    keep_columns = ["chrom", "bin_index", "ref_median"]
    for column in ("ref_sample_count", "ref_mad", "ref_mad_z", "ref_stability_label"):
        if column in frame.columns:
            keep_columns.append(column)
    return frame[keep_columns].drop_duplicates(["chrom", "bin_index"], keep="first")


def _expm1_series(series):
    values = pd.to_numeric(series, errors="coerce").astype(float)
    return np.expm1(values)


def derive_ratio_copy_number(frame, ref_bins_df=None):
    result = frame.copy()
    input_copy_number = pd.to_numeric(result["copy_number"], errors="coerce") if "copy_number" in result.columns else None
    input_log2r = pd.to_numeric(result["log2r"], errors="coerce") if "log2r" in result.columns else None
    result["log2r"] = np.nan
    result["raw_log2r"] = np.nan
    result["raw_ratio"] = np.nan
    result["ref_cpm"] = np.nan
    result["copy_number"] = np.nan
    result["copy_number_centering_log2_shift"] = np.nan
    result["copy_number_source"] = ""

    if input_copy_number is not None:
        copy_number = input_copy_number
        valid = np.isfinite(copy_number)
        result.loc[valid, "copy_number"] = copy_number[valid]
        result.loc[valid, "log2r"] = np.log2(copy_number[valid] / 2.0)
        result.loc[valid, "raw_log2r"] = result.loc[valid, "log2r"]
        result.loc[valid, "raw_ratio"] = copy_number[valid] / 2.0
        result.loc[valid, "copy_number_centering_log2_shift"] = 0.0
        result.loc[valid, "copy_number_source"] = "explicit_copy_number"
        if valid.any():
            return result

    if input_log2r is not None:
        log2r = input_log2r
        valid = np.isfinite(log2r)
        result.loc[valid, "log2r"] = log2r[valid]
        result.loc[valid, "raw_log2r"] = log2r[valid]
        result.loc[valid, "raw_ratio"] = np.power(2.0, log2r[valid])
        result.loc[valid, "copy_number"] = 2.0 * np.power(2.0, log2r[valid])
        result.loc[valid, "copy_number_centering_log2_shift"] = 0.0
        result.loc[valid, "copy_number_source"] = "explicit_log2r"
        if valid.any():
            return result

    if "normalized_signal" not in result.columns:
        raise ValueError("Copy-number plot requires log2r, copy_number, or normalized_signal plus reference bin medians.")

    ref_bins = coerce_ref_bins(ref_bins_df)
    if ref_bins.empty:
        raise ValueError("Copy-number plot requires reference bin medians when log2r/copy_number are not present.")

    ref_columns = [column for column in ("ref_median", "ref_mad", "ref_mad_z", "ref_stability_label") if column in result.columns]
    if ref_columns:
        result = result.drop(columns=ref_columns)
    result = result.merge(ref_bins, on=["chrom", "bin_index"], how="left")
    sample_cpm = _expm1_series(result["normalized_signal"])
    ref_cpm = _expm1_series(result["ref_median"])
    result.loc[np.isfinite(ref_cpm), "ref_cpm"] = ref_cpm.loc[np.isfinite(ref_cpm)]
    valid = np.isfinite(sample_cpm) & np.isfinite(ref_cpm) & (ref_cpm > 0.0)
    ratio = pd.Series(np.nan, index=result.index, dtype=float)
    ratio.loc[valid] = (sample_cpm.loc[valid] + CN_LOG2R_PSEUDOCOUNT) / (
        ref_cpm.loc[valid] + CN_LOG2R_PSEUDOCOUNT
    )
    valid_ratio = np.isfinite(ratio) & (ratio > 0.0)
    result.loc[valid_ratio, "raw_ratio"] = ratio.loc[valid_ratio]
    raw_log2r = pd.Series(np.nan, index=result.index, dtype=float)
    raw_log2r.loc[valid_ratio] = np.log2(ratio.loc[valid_ratio])
    centering_mask = (
        valid_ratio
        & result["chrom"].astype(str).isin(set(AUTOSOME_ORDER))
        & ~structural_gap_mask(result)
    )
    center_shift = float(raw_log2r.loc[centering_mask].median()) if centering_mask.any() else 0.0
    if not np.isfinite(center_shift):
        center_shift = 0.0
    centered_log2r = raw_log2r - center_shift
    result.loc[valid_ratio, "raw_log2r"] = raw_log2r.loc[valid_ratio]
    result.loc[valid_ratio, "log2r"] = centered_log2r.loc[valid_ratio]
    result.loc[valid_ratio, "copy_number"] = 2.0 * np.power(2.0, centered_log2r.loc[valid_ratio])
    result.loc[valid_ratio, "copy_number_centering_log2_shift"] = center_shift
    result.loc[valid_ratio, "copy_number_source"] = CN_BIN_COPY_NUMBER_SOURCE
    result.loc[~valid_ratio, "copy_number_source"] = "ref_median_unavailable"
    return result


def classify_copy_number_state(value):
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return "neutral"
    if not np.isfinite(numeric):
        return "neutral"
    if numeric < CN_NEUTRAL_LOWER:
        return "del"
    if numeric > CN_NEUTRAL_UPPER:
        return "dup"
    return "neutral"


def classify_sex_aware_copy_number_state(
    copy_number,
    expected_copy_number,
    chrom,
    sex_call,
    ref_cpm=np.nan,
    chrom_ref_cpm_median=np.nan,
):
    normalized = normalize_chrom(chrom)
    sex = str(sex_call or "").strip().upper()
    try:
        expected = float(expected_copy_number)
        value = float(copy_number)
    except (TypeError, ValueError):
        return "neutral", "copy_number_unavailable"
    if normalized == "chrY" and sex == "XY":
        try:
            ref_value = float(ref_cpm)
        except (TypeError, ValueError):
            ref_value = np.nan
        try:
            chrom_ref_value = float(chrom_ref_cpm_median)
        except (TypeError, ValueError):
            chrom_ref_value = np.nan
        if (
            not np.isfinite(ref_value)
            or ref_value < SEX_CHROM_REF_MIN_CPM
            or not np.isfinite(chrom_ref_value)
            or chrom_ref_value < SEX_CHROM_REF_MIN_CPM
        ):
            return "neutral", "sex_chrom_ref_ratio_not_interpretable"
    if not np.isfinite(expected):
        return "neutral", "sex_chrom_context_only"
    if expected == 0.0:
        return "neutral", "sex_aware_absent_expected"
    if not np.isfinite(value):
        return "neutral", "copy_number_unavailable"
    if expected <= 1.0:
        lower = expected * (CN_HAPLOID_NEUTRAL_LOWER / 1.0)
        upper = expected * (CN_HAPLOID_NEUTRAL_UPPER / 1.0)
    else:
        lower = expected * (CN_NEUTRAL_LOWER / 2.0)
        upper = expected * (CN_NEUTRAL_UPPER / 2.0)
    if value < lower:
        return "del", "sex_aware_interpretable"
    if value > upper:
        return "dup", "sex_aware_interpretable"
    return "neutral", "sex_aware_interpretable"


def annotate_copy_number_bins(bins_df, final_events, branch_s_events=None, ref_bins_df=None, sex_call=""):
    frame = derive_ratio_copy_number(bins_df.copy(), ref_bins_df=ref_bins_df)
    frame["cn_scatter_state"] = frame["copy_number"].map(classify_copy_number_state)
    if "is_par_region" in frame.columns:
        par_flag = boolean_like_series(frame["is_par_region"])
    elif "is_PAR" in frame.columns:
        par_flag = boolean_like_series(frame["is_PAR"])
    elif "par_overlap_fraction" in frame.columns:
        par_flag = pd.to_numeric(frame["par_overlap_fraction"], errors="coerce").fillna(0.0).gt(0.0)
    else:
        par_flag = pd.Series(False, index=frame.index)
    frame["sex_call"] = str(sex_call or "")
    frame["sex_chrom_region_class"] = [
        sex_chrom_region_class(chrom, is_par)
        for chrom, is_par in zip(frame["chrom"], par_flag)
    ]
    frame["expected_copy_number"] = [
        expected_copy_number_for_bin(chrom, sex_call, is_par)
        for chrom, is_par in zip(frame["chrom"], par_flag)
    ]
    raw_ratio = pd.to_numeric(frame.get("raw_ratio", pd.Series(np.nan, index=frame.index)), errors="coerce")
    ref_cpm = pd.to_numeric(frame.get("ref_cpm", pd.Series(np.nan, index=frame.index)), errors="coerce")
    frame["chrom_ref_cpm_median"] = frame.groupby("chrom")["ref_cpm"].transform(
        lambda values: pd.to_numeric(values, errors="coerce").dropna().median()
    )
    chrom_ref_cpm_median = pd.to_numeric(frame["chrom_ref_cpm_median"], errors="coerce")
    expected_numeric = pd.to_numeric(frame["expected_copy_number"], errors="coerce")
    y_haploid_interpretable = (
        frame["chrom"].astype(str).eq("chrY")
        & expected_numeric.eq(1.0)
        & raw_ratio.notna()
        & ref_cpm.ge(SEX_CHROM_REF_MIN_CPM)
        & chrom_ref_cpm_median.ge(SEX_CHROM_REF_MIN_CPM)
    )
    frame.loc[y_haploid_interpretable, "copy_number"] = raw_ratio.loc[y_haploid_interpretable] * expected_numeric.loc[
        y_haploid_interpretable
    ]
    frame.loc[y_haploid_interpretable, "log2r"] = np.log2(
        frame.loc[y_haploid_interpretable, "copy_number"] / expected_numeric.loc[y_haploid_interpretable]
    )
    frame.loc[y_haploid_interpretable, "copy_number_source"] = (
        "normalized_signal_ref_median_log2r_y_haploid_centered"
    )
    frame["copy_number_delta"] = pd.to_numeric(frame["copy_number"], errors="coerce") - pd.to_numeric(
        frame["expected_copy_number"], errors="coerce"
    )
    sex_aware = [
        classify_sex_aware_copy_number_state(cn, expected, chrom, sex_call, ref_cpm, chrom_ref)
        for cn, expected, chrom, ref_cpm, chrom_ref in zip(
            frame["copy_number"],
            frame["expected_copy_number"],
            frame["chrom"],
            frame.get("ref_cpm", pd.Series(np.nan, index=frame.index)),
            frame.get("chrom_ref_cpm_median", pd.Series(np.nan, index=frame.index)),
        )
    ]
    frame["cn_scatter_state_sex_aware"] = [item[0] for item in sex_aware]
    frame["copy_number_interpretation_status"] = [item[1] for item in sex_aware]
    not_interpretable = frame["copy_number_interpretation_status"].eq("sex_chrom_ref_ratio_not_interpretable")
    frame.loc[not_interpretable, "copy_number"] = np.nan
    frame.loc[not_interpretable, "log2r"] = np.nan
    frame.loc[not_interpretable, "copy_number_delta"] = np.nan
    frame.loc[not_interpretable, "copy_number_source"] = "sex_chrom_ref_ratio_not_interpretable"
    frame["report_state"] = frame["cn_scatter_state_sex_aware"]
    frame["event_report_state"] = "neutral"
    frame["event_layer"] = "neutral"
    missing_source = frame["copy_number_source"].fillna("").astype(str).eq("")
    frame.loc[~np.isfinite(frame["copy_number"]) & missing_source, "copy_number_source"] = "copy_number_unavailable"
    frame["is_structure_gap_blank"] = structural_gap_mask(frame)
    frame.loc[frame["is_structure_gap_blank"], "copy_number"] = np.nan
    frame.loc[frame["is_structure_gap_blank"], "log2r"] = np.nan
    frame.loc[frame["is_structure_gap_blank"], "cn_scatter_state"] = "neutral"
    frame.loc[frame["is_structure_gap_blank"], "cn_scatter_state_sex_aware"] = "neutral"
    frame.loc[frame["is_structure_gap_blank"], "report_state"] = "neutral"
    frame.loc[frame["is_structure_gap_blank"], "event_layer"] = "neutral"
    frame.loc[frame["is_structure_gap_blank"], "copy_number_source"] = "structure_gap_blank"
    event_frames = []
    if final_events is not None and not final_events.empty:
        event_frames.append(final_events.assign(event_layer="autosomal_report"))
    if branch_s_events is not None and not branch_s_events.empty:
        event_frames.append(branch_s_events)
    if not event_frames:
        return frame
    events = pd.concat(event_frames, ignore_index=True, sort=False)
    for row in events.itertuples(index=False):
        row_dict = row._asdict()
        mask = (
            frame["chrom"].astype(str).eq(str(row_dict["chrom"]))
            & frame["end"].astype(int).gt(int(row_dict["start"]))
            & frame["start"].astype(int).lt(int(row_dict["end"]))
            & ~frame["is_structure_gap_blank"]
            & frame["event_report_state"].eq("neutral")
        )
        frame.loc[mask, "event_report_state"] = str(row_dict["report_state"])
        frame.loc[mask, "event_layer"] = str(row_dict.get("event_layer", "autosomal_report"))
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


def render_event_region(row, layout, total_span, left, plot_width, y, height, color, label, class_name=""):
    chrom = str(row["chrom"])
    if chrom not in layout:
        return ""
    x1 = scale_x(genome_position(chrom, row["start"], layout), total_span, left, plot_width)
    x2 = scale_x(genome_position(chrom, row["end"], layout), total_span, left, plot_width)
    width = max(x2 - x1, 1.5)
    title = html.escape(label)
    class_attr = f' class="{html.escape(str(class_name))}"' if class_name else ""
    return (
        f'<rect{class_attr} x="{x1:.2f}" y="{y:.2f}" width="{width:.2f}" height="{height:.2f}" '
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
        valid_z = pd.to_numeric(event_bins.get("branch_a_ref_z", pd.Series(np.nan, index=event_bins.index)), errors="coerce")
        valid_z = valid_z[np.isfinite(valid_z)]
        if valid_z.empty:
            continue
        state = str(row_dict.get("report_state", "")).lower()
        if state == "dup":
            quantile = 0.75
            quantile_label = "Q75"
        elif state == "del":
            quantile = 0.25
            quantile_label = "Q25"
        else:
            continue
        trend_value = float(valid_z.quantile(quantile))
        x1 = scale_x(genome_position(chrom, int(row_dict.get("start")), layout), total_span, left, plot_width)
        x2 = scale_x(genome_position(chrom, int(row_dict.get("end")), layout), total_span, left, plot_width)
        if x2 <= x1:
            continue
        y = scale_y(trend_value, mid_y, half_h)
        score = first_finite_value(
            row_dict,
            ("a_zscore", "branch_a_zscore", "a_abs_zscore", "max_abs_zscore", "priority_score"),
        )
        score_text = f"; a_zscore={score:.3f}" if score is not None else ""
        label = (
            f"branch_a_ref_z {quantile_label} {row_dict.get('report_state')} "
            f"{chrom}:{int(row_dict.get('start'))}-{int(row_dict.get('end'))}; "
            f"{quantile_label} branch_a_ref_z={trend_value:.3f}; "
            f"valid bins={len(valid_z)}{score_text}"
        )
        chunks.append(
            f'<line class="report-ref-z-trend" x1="{x1:.2f}" y1="{y:.2f}" '
            f'x2="{x2:.2f}" y2="{y:.2f}" stroke="{TREND_COLOR}" '
            f'stroke-width="2.4" opacity="0.94" stroke-linecap="round">'
            f"<title>{html.escape(label)}</title></line>"
        )
    return chunks


def first_finite_value(row_dict, keys):
    for key in keys:
        if key not in row_dict:
            continue
        value = pd.to_numeric(pd.Series([row_dict.get(key)]), errors="coerce").iloc[0]
        if np.isfinite(value):
            return float(value)
    return None


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


def render_branch_s_cn_trend_lines(bins_df, branch_s_events, layout, total_span, left, plot_width, mid_y, half_h):
    chunks = []
    if branch_s_events is None or branch_s_events.empty:
        return chunks
    for row in branch_s_events.itertuples(index=False):
        row_dict = row._asdict()
        chrom = str(row_dict.get("chrom", ""))
        if chrom not in layout:
            continue
        event_bins = bins_df[
            bins_df["chrom"].astype(str).eq(chrom)
            & bins_df["end"].astype(int).gt(int(row_dict.get("start")))
            & bins_df["start"].astype(int).lt(int(row_dict.get("end")))
            & ~bins_df.get("is_structure_gap_blank", pd.Series(False, index=bins_df.index)).astype(bool)
            & bins_df.get("copy_number_interpretation_status", pd.Series("", index=bins_df.index)).astype(str).eq(
                "sex_aware_interpretable"
            )
        ].copy()
        if event_bins.empty:
            continue
        cn_values = pd.to_numeric(event_bins["copy_number"], errors="coerce")
        expected_values = pd.to_numeric(event_bins["expected_copy_number"], errors="coerce")
        event_bins = event_bins[np.isfinite(cn_values) & np.isfinite(expected_values)].copy()
        if event_bins.empty:
            continue
        trend_value = float(pd.to_numeric(event_bins["copy_number"], errors="coerce").median())
        expected = float(pd.to_numeric(event_bins["expected_copy_number"], errors="coerce").median())
        if not np.isfinite(trend_value) or not np.isfinite(expected) or expected <= 0.0:
            continue
        if math.isclose(expected, 2.0, rel_tol=0.0, abs_tol=0.05):
            threshold = 0.10
        elif math.isclose(expected, 1.0, rel_tol=0.0, abs_tol=0.05):
            threshold = 0.05
        else:
            continue
        if abs(trend_value - expected) < threshold:
            continue
        event_bins = event_bins.sort_values(["chrom", "start"]).copy()
        current = []
        previous_end = None
        for bin_row in event_bins.itertuples(index=False):
            if previous_end is not None and int(bin_row.start) > int(previous_end):
                chunks.extend(
                    render_branch_s_cn_trend_chunk(
                        current,
                        row_dict,
                        trend_value,
                        expected,
                        threshold,
                        layout,
                        total_span,
                        left,
                        plot_width,
                        mid_y,
                        half_h,
                    )
                )
                current = []
            current.append(bin_row)
            previous_end = int(bin_row.end)
        chunks.extend(
            render_branch_s_cn_trend_chunk(
                current,
                row_dict,
                trend_value,
                expected,
                threshold,
                layout,
                total_span,
                left,
                plot_width,
                mid_y,
                half_h,
            )
        )
    return chunks


def render_branch_s_cn_trend_chunk(
    chunk_rows,
    row_dict,
    trend_value,
    expected,
    threshold,
    layout,
    total_span,
    left,
    plot_width,
    mid_y,
    half_h,
):
    if not chunk_rows:
        return []
    chrom = str(row_dict.get("chrom", ""))
    x1 = scale_x(genome_position(chrom, int(chunk_rows[0].start), layout), total_span, left, plot_width)
    x2 = scale_x(genome_position(chrom, int(chunk_rows[-1].end), layout), total_span, left, plot_width)
    if x2 <= x1:
        return []
    y = scale_copy_number_y(trend_value, mid_y, half_h)
    state = str(row_dict.get("report_state", "")).lower()
    trend_color = CN_REPORT_STATE_COLOR.get(state, TREND_COLOR)
    label = (
        f"Branch S CN trend {row_dict.get('branch_s_state')} "
        f"{chrom}:{int(chunk_rows[0].start)}-{int(chunk_rows[-1].end)}; "
        f"median CN={trend_value:.3f}; expected CN={expected:.3f}; "
        f"threshold={threshold:.3f}"
    )
    return [
        f'<line class="branch-s-cn-trend" x1="{x1:.2f}" y1="{y:.2f}" '
        f'x2="{x2:.2f}" y2="{y:.2f}" stroke="{trend_color}" '
        f'stroke-width="2.4" opacity="0.94" stroke-linecap="round">'
        f"<title>{html.escape(label)}</title></line>"
    ]


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


def cn_support_mask(values, state):
    if str(state).lower() == "dup":
        return values > CN_NEUTRAL_UPPER
    if str(state).lower() == "del":
        return values < CN_NEUTRAL_LOWER
    return pd.Series(False, index=values.index)


def z_support_mask(values, state):
    if str(state).lower() == "dup":
        return values > 0.0
    if str(state).lower() == "del":
        return values < 0.0
    return pd.Series(False, index=values.index)


def summarize_copy_number_event_support(cn_bins, final_events):
    columns = [
        "sample_id",
        "chrom",
        "start",
        "end",
        "state",
        "event_report_state",
        "valid_bin_count",
        "cn_support_bin_count",
        "cn_same_direction_fraction",
        "median_bin_cn",
        "mean_bin_cn",
        "median_log2r",
        "median_calibrated_z",
        "z_support_bin_count",
        "centromere_gap_bin_count",
        "cn_direction_consistency_status",
    ]
    if final_events.empty:
        return pd.DataFrame(columns=columns)
    rows = []
    for event in final_events.itertuples(index=False):
        row = event._asdict()
        chrom = str(row.get("chrom", ""))
        start = int(row.get("start"))
        end = int(row.get("end"))
        state = str(row.get("report_state", ""))
        event_bins = cn_bins[
            cn_bins["chrom"].astype(str).eq(chrom)
            & cn_bins["end"].astype(int).gt(start)
            & cn_bins["start"].astype(int).lt(end)
        ].copy()
        gap_count = int(event_bins.get("is_structure_gap_blank", pd.Series(False, index=event_bins.index)).astype(bool).sum())
        valid = event_bins[
            ~event_bins.get("is_structure_gap_blank", pd.Series(False, index=event_bins.index)).astype(bool)
            & pd.to_numeric(event_bins["copy_number"], errors="coerce").notna()
        ].copy()
        cn_values = pd.to_numeric(valid["copy_number"], errors="coerce")
        log2r_values = pd.to_numeric(valid.get("log2r", pd.Series(dtype=float)), errors="coerce")
        z_values = pd.to_numeric(valid.get("z", pd.Series(dtype=float)), errors="coerce")
        cn_support = cn_support_mask(cn_values, state)
        z_support = z_support_mask(z_values, state)
        valid_count = int(cn_values.notna().sum())
        support_count = int(cn_support.fillna(False).sum())
        support_fraction = float(support_count / valid_count) if valid_count else np.nan
        if valid_count == 0:
            consistency = "CN_NO_VALID_BINS"
        elif support_fraction >= 0.50:
            consistency = "CN_DIRECTION_SUPPORTED"
        elif support_count > 0:
            consistency = "CN_DIRECTION_WEAK_OR_MIXED"
        else:
            consistency = "CN_DIRECTION_NOT_SUPPORTED"
        rows.append(
            {
                "sample_id": str(row.get("sample_id", "")),
                "chrom": chrom,
                "start": start,
                "end": end,
                "state": str(row.get("state", "")),
                "event_report_state": state,
                "valid_bin_count": valid_count,
                "cn_support_bin_count": support_count,
                "cn_same_direction_fraction": round(support_fraction, 6) if np.isfinite(support_fraction) else np.nan,
                "median_bin_cn": round(float(cn_values.median()), 6) if valid_count else np.nan,
                "mean_bin_cn": round(float(cn_values.mean()), 6) if valid_count else np.nan,
                "median_log2r": round(float(log2r_values.median()), 6) if log2r_values.notna().any() else np.nan,
                "median_calibrated_z": round(float(z_values.median()), 6) if z_values.notna().any() else np.nan,
                "z_support_bin_count": int(z_support.fillna(False).sum()),
                "centromere_gap_bin_count": gap_count,
                "cn_direction_consistency_status": consistency,
            }
        )
    return pd.DataFrame(rows, columns=columns)


def write_copy_number_event_support_tsv(path_value, support_df):
    if not path_value:
        return
    path = Path(path_value)
    path.parent.mkdir(parents=True, exist_ok=True)
    support_df.to_csv(path, sep="\t", index=False)


def write_plot_bins_tsv(path_value, bins_df, layout):
    if not path_value:
        return
    output = bins_df.copy()
    genome_positions = []
    for row in output.itertuples(index=False):
        center = int(row.start + ((row.end - row.start) / 2))
        genome_positions.append(int(genome_position(row.chrom, center, layout)))
    output["genome_pos"] = genome_positions
    columns = [
        "chrom",
        "start",
        "end",
        "genome_pos",
        "z",
        "branch_a_ref_z",
        "residual_calibrated_z",
        "z_source",
        "ref_z_scale",
        "ref_z_scale_source",
        "report_state",
    ]
    for column in columns:
        if column not in output.columns:
            output[column] = np.nan
    output = output[columns].copy()
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
    output = output[
        [
            "chrom",
            "start",
            "end",
            "genome_pos",
            "z",
            "raw_log2r",
            "log2r",
            "raw_ratio",
            "ref_cpm",
            "chrom_ref_cpm_median",
            "copy_number",
            "copy_number_centering_log2_shift",
            "cn_scatter_state",
            "cn_scatter_state_sex_aware",
            "report_state",
            "event_report_state",
            "event_layer",
            "sex_call",
            "expected_copy_number",
            "copy_number_delta",
            "sex_chrom_region_class",
            "copy_number_interpretation_status",
            "is_structure_gap_blank",
            "copy_number_source",
        ]
    ].copy()
    path = Path(path_value)
    path.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(path, sep="\t", index=False)


def build_copy_number_plot_svg(
    sample_id,
    bins,
    final_events,
    branch_s_events,
    layout,
    total_span,
    output_svg,
    output_bins_tsv="",
    output_event_support_tsv="",
    ref_bins_df=None,
    sex_call="",
    max_points=8000,
):
    cn_bins = annotate_copy_number_bins(
        bins,
        final_events,
        branch_s_events=branch_s_events,
        ref_bins_df=ref_bins_df,
        sex_call=sex_call,
    )
    cn_layout, cn_total_span = build_chrom_layout(cn_bins, gap_bp=CN_CHROM_GAP_BP)
    write_copy_number_bins_tsv(output_bins_tsv, cn_bins, cn_layout)
    write_copy_number_event_support_tsv(
        output_event_support_tsv,
        summarize_copy_number_event_support(cn_bins, final_events),
    )

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
        svg_text(24, 58, "Final reported dup/del regions over ratio-derived bin copy number", size=13, fill="#475569"),
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

    for row in final_events.to_dict("records"):
        color = CN_REPORT_STATE_COLOR.get(str(row.get("report_state", "")), NEUTRAL_COLOR)
        region = render_event_region(
            row,
            cn_layout,
            cn_total_span,
            left,
            plot_width,
            signal_top,
            signal_height,
            color,
            f"report CN region {row.get('report_state')} {row.get('chrom')}:{int(row.get('start'))}-{int(row.get('end'))}",
            class_name="report-cn-region",
        )
        if region:
            svg.append(region)
    if branch_s_events is not None and not branch_s_events.empty:
        for row in branch_s_events.to_dict("records"):
            color = CN_REPORT_STATE_COLOR.get(str(row.get("report_state", "")), NEUTRAL_COLOR)
            score = row.get("branch_s_score", np.nan)
            score_text = f"; event-level score={float(score):.3f}" if np.isfinite(float(score)) else ""
            region = render_event_region(
                row,
                cn_layout,
                cn_total_span,
                left,
                plot_width,
                signal_top,
                signal_height,
                color,
                (
                    f"Branch S review {row.get('branch_s_state')} "
                    f"{row.get('chrom')}:{int(row.get('start'))}-{int(row.get('end'))}"
                    f"{score_text}; not a bin-level CN value"
                ),
                class_name="branch-s-cn-region",
            )
            if region:
                svg.append(region)

    plot_bins = downsample_bins(cn_bins, max_points=max_points)
    for row in plot_bins.itertuples(index=False):
        if getattr(row, "is_structure_gap_blank", False) or not np.isfinite(float(row.copy_number)):
            continue
        x = scale_x(genome_position(row.chrom, int(row.start + ((row.end - row.start) / 2)), cn_layout), cn_total_span, left, plot_width)
        y = scale_copy_number_y(row.copy_number, mid_y, half_h)
        color = CN_REPORT_STATE_COLOR.get(str(row.cn_scatter_state_sex_aware), NEUTRAL_COLOR)
        opacity = 0.78 if row.cn_scatter_state_sex_aware in {"dup", "del"} else 0.58
        out_of_range_attr = (
            ' data-copy-number-out-of-range="true"'
            if float(row.copy_number) < CN_COPY_NUMBER_MIN or float(row.copy_number) > CN_COPY_NUMBER_MAX
            else ""
        )
        svg.append(
            f'<circle class="cn-bin-scatter"{out_of_range_attr} cx="{x:.2f}" cy="{y:.2f}" '
            f'r="{CN_SCATTER_RADIUS:.2f}" fill="{color}" opacity="{opacity:.2f}" '
            f'stroke="#ffffff" stroke-width="0.35"/>'
        )
    svg.extend(render_report_event_cn_trend_lines(cn_bins, final_events, cn_layout, cn_total_span, left, plot_width, mid_y, half_h))
    svg.extend(render_branch_s_cn_trend_lines(cn_bins, branch_s_events, cn_layout, cn_total_span, left, plot_width, mid_y, half_h))

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
    gender_df=None,
    branch_s_summary_df=None,
    branch_s_scores_df=None,
    branch_s_evidence_df=None,
    ref_bins_df=None,
    output_bins_tsv="",
    output_copy_number_svg="",
    output_copy_number_bins_tsv="",
    output_copy_number_event_support_tsv="",
    max_points=8000,
):
    bins = coerce_bins(bins_df)
    bins = annotate_branch_a_ref_z_bins(bins, ref_bins_df)
    final_events = coerce_final_events(branch_b_events_df, sample_id=sample_id)
    sex_call = sample_sex_call(gender_df, branch_s_summary_df, sample_id=sample_id)
    branch_s_events = coerce_branch_s_review_events(
        branch_s_summary_df,
        scores_df=branch_s_scores_df,
        evidence_df=branch_s_evidence_df,
        sample_id=sample_id,
    )
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
        svg_text(24, 58, "Final reported dup/del regions over Branch A ref-normalized z signal", size=13, fill="#475569"),
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
    if not branch_s_events.empty:
        for row in branch_s_events.to_dict("records"):
            color = REPORT_STATE_COLOR.get(str(row.get("report_state", "")), NEUTRAL_COLOR)
            score = row.get("branch_s_score", np.nan)
            score_text = f"; event-level score={float(score):.3f}" if np.isfinite(float(score)) else ""
            label = (
                f"Branch S review {row.get('branch_s_state')} "
                f"{row.get('chrom')}:{int(row.get('start'))}-{int(row.get('end'))}"
                f"{score_text}; branch_a_ref_z points remain bin-level"
            )
            svg.append(
                render_event_region(
                    row,
                    layout,
                    total_span,
                    left,
                    plot_width,
                    signal_top,
                    signal_height,
                    color,
                    label,
                    class_name="branch-s-review-region",
                )
            )

    plot_bins = downsample_bins(bins, max_points=max_points)
    point_chunks = []
    for row in plot_bins.itertuples(index=False):
        if not np.isfinite(float(row.z)):
            continue
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
        ("event ref-z trend", TREND_COLOR),
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
            branch_s_events=branch_s_events,
            layout=layout,
            total_span=total_span,
            output_svg=output_copy_number_svg,
            output_bins_tsv=output_copy_number_bins_tsv,
            output_event_support_tsv=output_copy_number_event_support_tsv,
            ref_bins_df=ref_bins_df,
            sex_call=sex_call,
            max_points=max_points,
        )


def main():
    args = parse_args()
    logger = setup_logger("cnv_branch_b_plot", args.log or None)
    bins_df = read_table(args.input_bins, empty_ok=False)
    events_df = read_table(args.input_events, empty_ok=True)
    a_branch_df = read_table(args.input_a_branch, empty_ok=True)
    ref_bins_df = read_table(args.input_ref_bins, empty_ok=True)
    gender_df = read_table(args.gender_tsv, empty_ok=True)
    branch_s_summary_df = read_table(args.branch_s_summary, empty_ok=True)
    branch_s_scores_df = read_table(args.branch_s_scores, empty_ok=True)
    branch_s_evidence_df = read_table(args.branch_s_evidence, empty_ok=True)
    build_cnv_plot_svg(
        sample_id=args.sample_id,
        bins_df=bins_df,
        branch_b_events_df=events_df,
        a_branch_df=a_branch_df,
        gender_df=gender_df,
        branch_s_summary_df=branch_s_summary_df,
        branch_s_scores_df=branch_s_scores_df,
        branch_s_evidence_df=branch_s_evidence_df,
        ref_bins_df=ref_bins_df,
        output_svg=args.output_svg,
        output_bins_tsv=args.output_bins_tsv,
        output_copy_number_svg=args.output_copy_number_svg,
        output_copy_number_bins_tsv=args.output_copy_number_bins_tsv,
        output_copy_number_event_support_tsv=args.output_copy_number_event_support_tsv,
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
