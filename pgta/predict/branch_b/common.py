#!/biosoftware/miniconda/envs/snakemake_env/bin/python
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd


CHROM_ORDER = tuple(str(index) for index in range(1, 23)) + ("X", "Y")
CHROM_TO_LABEL = {chrom: f"chr{chrom}" for chrom in CHROM_ORDER}
SEX_CHROM_NUMERIC_ALIASES = {"X": ("23", 23), "Y": ("24", 24)}
STATE_TO_SVTYPE = {"loss": "DEL", "gain": "DUP", "neutral": "NEU"}
OPTIONAL_REGION_FRACTION_COLUMNS = {
    "par_overlap_fraction": 0.0,
    "xtr_overlap_fraction": 0.0,
    "sex_homology_overlap_fraction": 0.0,
    "segmental_duplication_overlap_fraction": 0.0,
    "low_mappability_overlap_fraction": 0.0,
    "gap_centromere_telomere_overlap_fraction": 0.0,
    "repeat_rich_overlap_fraction": 0.0,
    "blacklist_overlap_fraction": 0.0,
    "ambiguous_alignment_overlap_fraction": 0.0,
}
OPTIONAL_REGION_BOOL_COLUMNS = {
    "is_PAR": 0,
    "is_XTR": 0,
    "is_sex_homology": 0,
    "is_segmental_duplication": 0,
    "is_low_mappability": 0,
    "is_gap_centromere_telomere": 0,
    "is_repeat_rich": 0,
    "is_blacklist_region": 0,
    "is_ambiguous_alignment_region": 0,
}


@dataclass
class CandidateEvent:
    event_id: str
    sample_id: str
    branch: str
    correction_model: str
    caller: str
    caller_stage: str
    chrom: str
    start: int
    end: int
    start_bin: int
    end_bin: int
    n_bins: int
    state: str
    svtype: str
    segment_id: str
    segment_weight: float
    segment_mean_signal: float
    segment_median_signal: float
    segment_mean_robust_z: float
    segment_median_robust_z: float
    segment_abs_max_robust_z: float
    calibrated_mean_z: float = np.nan
    calibrated_median_z: float = np.nan
    empirical_pvalue: float = np.nan
    empirical_qvalue: float = np.nan
    artifact_status: str = "unreviewed"
    artifact_flags: str = ""
    artifact_explanations: str = ""
    keep_event: int = 1

    def to_dict(self):
        payload = asdict(self)
        for key, value in list(payload.items()):
            if isinstance(value, float) and not np.isfinite(value):
                payload[key] = None
        return payload


def _to_numeric_vector(raw, npz_path, chrom):
    if raw is None:
        return None
    arr = np.asarray(raw)
    if arr.size == 0:
        return None
    if not np.issubdtype(arr.dtype, np.number):
        raise ValueError(f"Non-numeric counts for chr{chrom} in {npz_path}")
    vec = np.nan_to_num(arr.astype(np.float64).reshape(-1), nan=0.0, posinf=0.0, neginf=0.0)
    return vec if vec.size > 0 else None


def _coerce_quality_metric(raw):
    if raw is None:
        return np.nan
    arr = np.asarray(raw)
    if arr.shape == () and np.issubdtype(arr.dtype, np.number):
        return float(arr.item())
    return np.nan


def _sample_vector_for_chrom(sample_obj, chrom, npz_path):
    candidate_keys = [chrom]
    if str(chrom).isdigit():
        candidate_keys.append(int(chrom))
    candidate_keys.extend(SEX_CHROM_NUMERIC_ALIASES.get(str(chrom), ()))
    for key in candidate_keys:
        if key not in sample_obj:
            continue
        vec = _to_numeric_vector(sample_obj.get(key), npz_path, chrom)
        if vec is not None:
            return vec
    return None


def load_sample_bins(npz_path):
    npz_path = Path(npz_path)
    with np.load(npz_path, allow_pickle=True, encoding="latin1") as data:
        if "sample" not in data.files:
            raise ValueError(f"Missing 'sample' key in {npz_path}")
        sample_obj = data["sample"].item()
        if not isinstance(sample_obj, dict):
            raise ValueError(f"'sample' is not dict-like in {npz_path}")
        binsize = int(np.asarray(data["binsize"]).item()) if "binsize" in data.files else 0
        quality = _coerce_quality_metric(data["quality"]) if "quality" in data.files else np.nan
        rows = []
        for chrom in CHROM_ORDER:
            vec = _sample_vector_for_chrom(sample_obj, chrom, npz_path)
            if vec is None:
                continue
            label = CHROM_TO_LABEL[chrom]
            for bin_index, count in enumerate(vec):
                start = int(bin_index * binsize)
                rows.append(
                    {
                        "chrom": label,
                        "chrom_key": chrom,
                        "bin_index": int(bin_index),
                        "start": start,
                        "end": int(start + binsize),
                        "raw_count": float(count),
                    }
                )
    if not rows:
        raise ValueError(f"No usable chromosome vectors found in {npz_path}")
    bins_df = pd.DataFrame(rows)
    total = float(bins_df["raw_count"].sum())
    total = total if total > 0.0 else 1.0
    bins_df["normalized_signal"] = np.log1p((bins_df["raw_count"] / total) * 1e6)
    bins_df["bin_weight"] = np.sqrt(np.clip(bins_df["raw_count"], a_min=1.0, a_max=None))
    return bins_df, binsize, quality


def robust_z(values):
    values = np.asarray(values, dtype=np.float64)
    if values.size == 0:
        return np.zeros(0, dtype=np.float64)
    median = np.median(values)
    mad = np.median(np.abs(values - median))
    if not np.isfinite(mad) or mad <= 0.0:
        std = np.std(values)
        if not np.isfinite(std) or std <= 0.0:
            return np.zeros_like(values)
        return (values - median) / std
    return (values - median) / (1.4826 * mad)


def interval_overlap(start, end, left, right):
    return max(0, min(int(end), int(right)) - max(int(start), int(left)))


def ensure_region_annotation_columns(df):
    frame = df.copy()
    for column, default in OPTIONAL_REGION_FRACTION_COLUMNS.items():
        if column not in frame.columns:
            frame[column] = default
        frame[column] = pd.to_numeric(frame[column], errors="coerce").fillna(default).clip(lower=0.0, upper=1.0)
    for column, default in OPTIONAL_REGION_BOOL_COLUMNS.items():
        if column not in frame.columns:
            frame[column] = default
        frame[column] = frame[column].fillna(default).astype(int)
    if "mask_label" not in frame.columns:
        frame["mask_label"] = "pass"
    else:
        frame["mask_label"] = frame["mask_label"].fillna("pass").astype(str)
    if "is_autosome" not in frame.columns:
        frame["is_autosome"] = frame["chrom"].astype(str).str.match(r"^chr([1-9]|1[0-9]|2[0-2])$")
    frame["is_autosome"] = frame["is_autosome"].fillna(0).astype(int)
    if "mappability_score" not in frame.columns:
        if "mappability_proxy" in frame.columns:
            source = frame["mappability_proxy"]
        elif "atgc_fraction" in frame.columns:
            source = frame["atgc_fraction"]
        else:
            source = pd.Series(1.0, index=frame.index)
        frame["mappability_score"] = pd.to_numeric(source, errors="coerce").fillna(1.0).clip(lower=0.0, upper=1.0)
    else:
        frame["mappability_score"] = pd.to_numeric(frame["mappability_score"], errors="coerce").fillna(1.0).clip(
            lower=0.0, upper=1.0
        )
    return frame


def annotate_region_risk(df):
    frame = ensure_region_annotation_columns(df)

    mapping_component = np.maximum(
        frame["low_mappability_overlap_fraction"].to_numpy(dtype=np.float64),
        np.clip((0.85 - frame["mappability_score"].to_numpy(dtype=np.float64)) / 0.50, 0.0, 1.0),
    )
    repeat_component = np.maximum.reduce(
        [
            frame["segmental_duplication_overlap_fraction"].to_numpy(dtype=np.float64),
            frame["repeat_rich_overlap_fraction"].to_numpy(dtype=np.float64) * 0.85,
            frame["blacklist_overlap_fraction"].to_numpy(dtype=np.float64),
            frame["gap_centromere_telomere_overlap_fraction"].to_numpy(dtype=np.float64),
        ]
    )
    homology_component = np.maximum.reduce(
        [
            frame["par_overlap_fraction"].to_numpy(dtype=np.float64),
            frame["xtr_overlap_fraction"].to_numpy(dtype=np.float64),
            frame["sex_homology_overlap_fraction"].to_numpy(dtype=np.float64),
            frame["ambiguous_alignment_overlap_fraction"].to_numpy(dtype=np.float64),
        ]
    )
    mask_penalty = frame["mask_label"].map({"pass": 0.0, "soft": 0.20, "dynamic": 0.55, "hard": 0.95}).fillna(0.0)
    risk_score = (
        0.35 * mapping_component
        + 0.30 * repeat_component
        + 0.25 * homology_component
        + 0.10 * mask_penalty.to_numpy(dtype=np.float64)
    )
    risk_score = np.maximum(
        risk_score,
        np.maximum.reduce(
            [
                frame["blacklist_overlap_fraction"].to_numpy(dtype=np.float64),
                frame["gap_centromere_telomere_overlap_fraction"].to_numpy(dtype=np.float64),
                frame["ambiguous_alignment_overlap_fraction"].to_numpy(dtype=np.float64) * 0.95,
            ]
        ),
    )
    risk_score = np.clip(risk_score, 0.0, 1.0)

    risk_class = np.full(len(frame), "clean", dtype=object)
    moderate_condition = (
        (risk_score >= 0.35)
        | (homology_component > 0.0)
        | (repeat_component > 0.0)
        | (mapping_component > 0.0)
        | frame["mask_label"].isin({"soft", "dynamic"}).to_numpy(dtype=bool)
    )
    high_condition = (
        (risk_score >= 0.75)
        | (frame["blacklist_overlap_fraction"].to_numpy(dtype=np.float64) > 0.0)
        | (frame["gap_centromere_telomere_overlap_fraction"].to_numpy(dtype=np.float64) > 0.0)
        | (frame["ambiguous_alignment_overlap_fraction"].to_numpy(dtype=np.float64) >= 0.25)
        | (frame["mask_label"].isin({"hard", "dynamic"}).to_numpy(dtype=bool))
    )
    risk_class[moderate_condition] = "moderate"
    risk_class[high_condition] = "high"

    region_risk_weight = np.clip(np.exp(-2.0 * risk_score), 0.10, 1.0)
    variance_inflation = 1.0 + (2.5 * risk_score) + np.where(high_condition, 0.5, 0.0)

    fit_block = (
        (frame["blacklist_overlap_fraction"].to_numpy(dtype=np.float64) > 0.0)
        | (frame["gap_centromere_telomere_overlap_fraction"].to_numpy(dtype=np.float64) > 0.0)
        | (frame["ambiguous_alignment_overlap_fraction"].to_numpy(dtype=np.float64) > 0.0)
        | (mapping_component >= 0.75)
        | (frame["segmental_duplication_overlap_fraction"].to_numpy(dtype=np.float64) >= 0.75)
        | (homology_component >= 0.25)
    )
    seed_block = (
        high_condition
        | (frame["gap_centromere_telomere_overlap_fraction"].to_numpy(dtype=np.float64) > 0.0)
        | (frame["blacklist_overlap_fraction"].to_numpy(dtype=np.float64) > 0.0)
        | (frame["ambiguous_alignment_overlap_fraction"].to_numpy(dtype=np.float64) > 0.0)
    )
    null_block = (
        ~frame["is_autosome"].to_numpy(dtype=bool)
        | high_condition
        | moderate_condition
        | frame["mask_label"].isin({"soft", "dynamic", "hard"}).to_numpy(dtype=bool)
    )

    frame["region_risk_score"] = risk_score
    frame["region_risk_class"] = risk_class
    frame["region_risk_weight"] = region_risk_weight
    frame["variance_inflation"] = variance_inflation
    frame["correction_fit_eligible"] = (~fit_block & ~frame["mask_label"].isin({"hard", "dynamic"})).astype(int)
    frame["calling_seed_eligible"] = (~seed_block & ~frame["mask_label"].isin({"hard"})).astype(int)
    frame["calibration_null_eligible"] = (~null_block).astype(int)
    return frame


def weighted_mean(values, weights):
    values = np.asarray(values, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    valid = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    if not np.any(valid):
        return np.nan
    return float(np.average(values[valid], weights=weights[valid]))


def effective_sample_size(weights):
    weights = np.asarray(weights, dtype=np.float64)
    valid = np.isfinite(weights) & (weights > 0.0)
    if not np.any(valid):
        return 0.0
    numer = np.square(weights[valid].sum())
    denom = np.square(weights[valid]).sum()
    if denom <= 0.0:
        return 0.0
    return float(numer / denom)


def stouffer_weighted_z(values, weights):
    values = np.asarray(values, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    valid = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    if not np.any(valid):
        return np.nan
    denom = np.sqrt(np.square(weights[valid]).sum())
    if denom <= 0.0:
        return np.nan
    return float(np.dot(weights[valid], values[valid]) / denom)


def estimate_autocorrelation_by_chrom(df, value_column, max_lag=5):
    if df.empty or value_column not in df.columns:
        return {}
    correlations = {lag: [] for lag in range(1, max(1, int(max_lag)) + 1)}
    for _, chrom_df in df.groupby("chrom", sort=False):
        series = pd.to_numeric(chrom_df[value_column], errors="coerce").to_numpy(dtype=np.float64)
        series = series[np.isfinite(series)]
        if series.size < 3:
            continue
        series = series - np.mean(series)
        variance = np.var(series)
        if not np.isfinite(variance) or variance <= 0.0:
            continue
        for lag in correlations:
            if series.size <= lag:
                continue
            left = series[:-lag]
            right = series[lag:]
            corr = float(np.dot(left, right) / (len(left) * variance))
            if np.isfinite(corr):
                correlations[lag].append(corr)
    return {
        lag: float(np.clip(np.mean(values), -0.95, 0.95))
        for lag, values in correlations.items()
        if values
    }


def autocorrelation_inflation(event_length, lag_correlations):
    n = max(int(round(event_length)), 1)
    if n <= 1 or not lag_correlations:
        return 1.0
    inflation = 1.0
    for lag, rho in lag_correlations.items():
        if lag >= n:
            continue
        inflation += 2.0 * float(rho) * (1.0 - (float(lag) / float(n)))
    return float(max(inflation, 1.0))


def benjamini_hochberg(pvalues):
    pvalues = np.asarray(pvalues, dtype=np.float64)
    result = np.full_like(pvalues, np.nan)
    valid = np.isfinite(pvalues)
    if not np.any(valid):
        return result
    order = np.argsort(pvalues[valid])
    ranked = pvalues[valid][order]
    n = ranked.size
    adjusted = np.empty(n, dtype=np.float64)
    running = 1.0
    for index in range(n - 1, -1, -1):
        value = ranked[index] * n / float(index + 1)
        running = min(running, value)
        adjusted[index] = running
    valid_indices = np.where(valid)[0][order]
    result[valid_indices] = np.clip(adjusted, 0.0, 1.0)
    return result


def write_table(path_value, df):
    path = Path(path_value)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def read_table(path_value, required_columns=None, empty_ok=True):
    df = pd.read_csv(path_value, sep="\t")
    missing = [column for column in (required_columns or []) if column not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in {path_value}: {','.join(missing)}")
    if not empty_ok and df.empty:
        raise ValueError(f"Input table is empty: {path_value}")
    return df


def read_bins_and_candidates(
    bins_path,
    candidates_path,
    bins_required_columns=None,
    candidate_required_columns=None,
    empty_candidates_ok=True,
):
    bins_df = read_table(bins_path, required_columns=bins_required_columns, empty_ok=False)
    candidates_df = read_table(
        candidates_path,
        required_columns=candidate_required_columns,
        empty_ok=empty_candidates_ok,
    )
    return bins_df, candidates_df


def write_json(path_value, payload):
    path = Path(path_value)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
