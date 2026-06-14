from __future__ import annotations

import math


BROAD_GAIN_CALLERS = {"chromosome_dosage_detector", "raw_chromosome_dosage_detector"}


def normalize_chrom(value):
    chrom = str(value).strip()
    if not chrom:
        return chrom
    return chrom if chrom.startswith("chr") else f"chr{chrom}"


def normalize_state(value):
    text = str(value).strip().lower()
    mapping = {
        "gain": "gain",
        "dup": "gain",
        "duplication": "gain",
        "amplification": "gain",
        "+": "gain",
        "loss": "loss",
        "del": "loss",
        "deletion": "loss",
        "-": "loss",
    }
    return mapping.get(text, text)


def filter_truth_to_samples(truth_df, sample_ids):
    if truth_df is None or truth_df.empty or "sample_id" not in truth_df.columns:
        return truth_df
    sample_set = {str(sample_id) for sample_id in sample_ids if str(sample_id)}
    if not sample_set:
        return truth_df.iloc[0:0].copy()
    return truth_df[truth_df["sample_id"].astype(str).isin(sample_set)].copy()


def _is_missing(value):
    try:
        return bool(math.isnan(float(value)))
    except (TypeError, ValueError):
        return value is None


def overlap_fraction(left_start, left_end, right_start, right_end):
    if any(_is_missing(value) for value in [left_start, left_end, right_start, right_end]):
        return float("nan")
    overlap = max(0.0, min(float(left_end), float(right_end)) - max(float(left_start), float(right_start)) + 1.0)
    length = max(float(left_end) - float(left_start) + 1.0, 1.0)
    return overlap / length


def event_support_z(frame):
    import numpy as np
    import pandas as pd

    support_columns = [
        column for column in ["calibrated_mean_z", "calibrated_median_z", "event_corr_adjusted_z"] if column in frame.columns
    ]
    if frame.empty:
        return pd.Series(np.nan, index=frame.index, dtype=float)
    support_frame = (
        frame[support_columns].apply(pd.to_numeric, errors="coerce").abs()
        if support_columns
        else pd.DataFrame(index=frame.index)
    )
    if "a_abs_zscore" in frame.columns:
        a_support = pd.to_numeric(frame["a_abs_zscore"], errors="coerce").abs()
        caller = _text_column(frame, "caller")
        support_frame["a_abs_zscore"] = a_support.where(caller.eq("wisecondorx_a_branch"))
    if support_frame.empty:
        return pd.Series(np.nan, index=frame.index, dtype=float)
    return support_frame.max(axis=1, skipna=True)


def _numeric_column(frame, column, default):
    import pandas as pd

    if column not in frame.columns:
        return pd.Series(default, index=frame.index, dtype=float)
    return pd.to_numeric(frame[column], errors="coerce").fillna(default)


def _text_column(frame, column):
    import pandas as pd

    if column not in frame.columns:
        return pd.Series("", index=frame.index, dtype=str)
    return frame[column].fillna("").astype(str)


def _autosome_mask(chrom_series):
    return chrom_series.astype(str).str.match(r"^chr([1-9]|1[0-9]|2[0-2])$")


def broad_low_level_gain_mask(
    frame,
    broad_gain_qvalue=0.10,
    broad_gain_min_effective_bins=10.0,
    broad_gain_min_support_fraction=0.35,
):
    import pandas as pd

    if frame.empty:
        return pd.Series(False, index=frame.index, dtype=bool)

    chrom = _text_column(frame, "chrom").map(normalize_chrom)
    state = _text_column(frame, "state").map(normalize_state)
    caller = _text_column(frame, "caller")
    qvalue = _numeric_column(frame, "empirical_qvalue", math.nan)
    effective_bins = _numeric_column(frame, "effective_bin_count", math.nan)
    n_bins = _numeric_column(frame, "n_bins", math.nan)
    clean_fraction = _numeric_column(frame, "clean_bin_fraction", 0.0)
    high_fraction = _numeric_column(frame, "high_risk_bin_fraction", 0.0)
    moderate_fraction = _numeric_column(frame, "moderate_risk_bin_fraction", math.nan)
    moderate_fraction = moderate_fraction.where(
        moderate_fraction.notna(),
        (1.0 - clean_fraction - high_fraction).clip(lower=0.0),
    )
    weighted_non_high_support = clean_fraction + (0.5 * moderate_fraction.clip(lower=0.0))
    keep_event = _numeric_column(frame, "keep_event", 1.0)
    broad_flag = _text_column(frame, "artifact_flags").str.contains("broad_chrom_fraction", regex=False)

    return (
        keep_event.astype(float).eq(1.0)
        & _autosome_mask(chrom)
        & state.eq("gain")
        & caller.isin(BROAD_GAIN_CALLERS)
        & (broad_flag | caller.isin(BROAD_GAIN_CALLERS))
        & (effective_bins.fillna(n_bins) >= float(broad_gain_min_effective_bins))
        & (weighted_non_high_support >= float(broad_gain_min_support_fraction))
        & (high_fraction < 0.50)
        & (~qvalue.notna() | (qvalue <= float(broad_gain_qvalue)))
    )


def event_detection_threshold(
    frame,
    standard_z_threshold,
    broad_gain_z_threshold=1.8,
    broad_gain_qvalue=0.10,
):
    import pandas as pd

    threshold = pd.Series(float(standard_z_threshold), index=frame.index, dtype=float)
    broad_mask = broad_low_level_gain_mask(frame, broad_gain_qvalue=broad_gain_qvalue)
    threshold.loc[broad_mask] = min(float(standard_z_threshold), float(broad_gain_z_threshold))
    return threshold


def event_detected(
    frame,
    standard_z_threshold,
    broad_gain_z_threshold=1.8,
    broad_gain_qvalue=0.10,
):
    support = event_support_z(frame)
    threshold = event_detection_threshold(
        frame,
        standard_z_threshold=standard_z_threshold,
        broad_gain_z_threshold=broad_gain_z_threshold,
        broad_gain_qvalue=broad_gain_qvalue,
    )
    return support.ge(threshold)
