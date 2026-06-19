#!/biosoftware/miniconda/envs/snakemake_env/bin/python
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from pgta.core.logging import setup_logger
from pgta.predict.branch_b.common import read_table, write_json, write_table


BRANCH_S_VERSION = "branch_s_shadow_v1"
EVIDENCE_COLUMNS = [
    "branch_s_version",
    "sample_id",
    "sex_call",
    "predict_gender",
    "sex_call_source",
    "region_class",
    "chrom",
    "is_par_region",
    "bin_count",
    "start",
    "end",
    "span_bp",
    "a_candidate_count",
    "mean_calibrated_z",
    "median_calibrated_z",
    "mean_robust_z",
    "median_robust_z",
    "abs_mean_calibrated_z",
    "par_bin_fraction",
    "branch_s_action",
    "final_report_impact",
]
SCORE_COLUMNS = [
    "branch_s_version",
    "sample_id",
    "sex_call",
    "predict_gender",
    "sca_state",
    "source_region_class",
    "state_score",
    "state_score_status",
    "state_score_reason",
    "branch_s_action",
    "final_report_impact",
]


def parse_args():
    parser = argparse.ArgumentParser(description="Build Branch S sex-chromosome shadow evidence.")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--input-bins", required=True)
    parser.add_argument("--a-candidates", required=True)
    parser.add_argument("--gender-tsv", default="")
    parser.add_argument("--output-evidence", required=True)
    parser.add_argument("--output-scores", required=True)
    parser.add_argument("--output-summary", required=True)
    parser.add_argument("--version", default=BRANCH_S_VERSION)
    parser.add_argument("--log", default="")
    return parser.parse_args()


def _empty_table(columns):
    return pd.DataFrame(columns=list(columns))


def read_optional_table(path_value, columns=None):
    if not path_value:
        return _empty_table(columns or [])
    path = Path(path_value)
    if not path.exists() or path.stat().st_size == 0:
        return _empty_table(columns or [])
    return read_table(path, empty_ok=True)


def normalize_chrom(value):
    token = str(value).strip()
    if token.lower().startswith("chr"):
        token = token[3:]
    token = token.upper()
    if token in {"23", "X"}:
        return "chrX"
    if token in {"24", "Y"}:
        return "chrY"
    return f"chr{token}" if token else ""


def _coerce_bool_series(series):
    if series is None:
        return None
    if series.dtype == bool:
        return series.fillna(False)
    text = series.fillna("").astype(str).str.strip().str.lower()
    return text.isin({"1", "true", "t", "yes", "y"})


def prepare_bins(bins):
    frame = bins.copy()
    if frame.empty:
        for column, default in {"is_par_region": False, "calibrated_z": np.nan, "robust_z": np.nan}.items():
            if column not in frame.columns:
                frame[column] = default
        return frame
    required = {"chrom", "start", "end"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Branch S input bins missing required columns: {','.join(sorted(missing))}")
    frame["chrom"] = frame["chrom"].map(normalize_chrom)
    frame["start"] = pd.to_numeric(frame["start"], errors="coerce")
    frame["end"] = pd.to_numeric(frame["end"], errors="coerce")
    frame = frame[frame["chrom"].isin({"chrX", "chrY"}) & frame["start"].notna() & frame["end"].notna()].copy()
    if frame.empty:
        for column, default in {"is_par_region": False, "calibrated_z": np.nan, "robust_z": np.nan}.items():
            if column not in frame.columns:
                frame[column] = default
        return frame
    if "is_PAR" in frame.columns:
        par_flag = _coerce_bool_series(frame["is_PAR"])
    else:
        par_flag = pd.Series(False, index=frame.index)
    if "par_overlap_fraction" in frame.columns:
        par_overlap = pd.to_numeric(frame["par_overlap_fraction"], errors="coerce").fillna(0.0)
        par_flag = par_flag | (par_overlap > 0.0)
    frame["is_par_region"] = par_flag.astype(bool)
    for column in ["calibrated_z", "robust_z"]:
        if column not in frame.columns:
            frame[column] = np.nan
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame["start"] = frame["start"].astype(int)
    frame["end"] = frame["end"].astype(int)
    return frame


def prepare_candidates(candidates):
    frame = candidates.copy()
    if frame.empty:
        return frame
    required = {"chrom", "start", "end"}
    if not required.issubset(frame.columns):
        return pd.DataFrame(columns=["chrom", "start", "end"])
    frame["chrom"] = frame["chrom"].map(normalize_chrom)
    frame["start"] = pd.to_numeric(frame["start"], errors="coerce")
    frame["end"] = pd.to_numeric(frame["end"], errors="coerce")
    if "state" not in frame.columns and "cnv_type" in frame.columns:
        frame["state"] = frame["cnv_type"]
    if "state" in frame.columns:
        frame["state"] = frame["state"].map(normalize_candidate_state)
    for column in ["a_zscore", "a_abs_zscore"]:
        if column in frame.columns:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
    return frame[frame["chrom"].isin({"chrX", "chrY"}) & frame["start"].notna() & frame["end"].notna()].copy()


def normalize_candidate_state(value):
    text = str(value).strip().lower()
    if text in {"gain", "dup", "duplication"}:
        return "gain"
    if text in {"loss", "del", "deletion"}:
        return "loss"
    return ""


def gender_context(gender):
    if gender is None or gender.empty:
        return {"sex_call": "", "predict_gender": "", "sex_call_source": ""}
    row = gender.iloc[0]
    return {
        "sex_call": str(row.get("sex_call", "") or ""),
        "predict_gender": str(row.get("predict_gender", "") or ""),
        "sex_call_source": str(row.get("sex_call_source", "") or ""),
    }


def finite_mean(values):
    vector = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(dtype=np.float64)
    vector = vector[np.isfinite(vector)]
    return float(np.mean(vector)) if vector.size else np.nan


def finite_median(values):
    vector = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(dtype=np.float64)
    vector = vector[np.isfinite(vector)]
    return float(np.median(vector)) if vector.size else np.nan


def count_overlapping_candidates(region_bins, candidates):
    return int(len(overlapping_candidate_rows(region_bins, candidates)))


def overlapping_candidate_rows(region_bins, candidates):
    if region_bins.empty or candidates is None or candidates.empty:
        return pd.DataFrame()
    chrom = str(region_bins["chrom"].iat[0])
    region_candidates = candidates[candidates["chrom"].astype(str).eq(chrom)].copy()
    if region_candidates.empty:
        return pd.DataFrame()
    matched = np.zeros(len(region_candidates), dtype=bool)
    candidate_starts = pd.to_numeric(region_candidates["start"], errors="coerce").to_numpy(dtype=np.float64)
    candidate_ends = pd.to_numeric(region_candidates["end"], errors="coerce").to_numpy(dtype=np.float64)
    for _, row in region_bins.iterrows():
        matched |= (candidate_starts < int(row["end"])) & (candidate_ends > int(row["start"]))
    return region_candidates.loc[matched].copy()


def region_class(chrom, is_par):
    axis = "X" if str(chrom) == "chrX" else "Y"
    return f"{axis}_PAR" if bool(is_par) else f"{axis}_NONPAR"


def summarize_region(sample_id, region_bins, candidates, context, version):
    chrom = str(region_bins["chrom"].iat[0])
    is_par = bool(region_bins["is_par_region"].iat[0])
    mean_calibrated_z = finite_mean(region_bins["calibrated_z"])
    start = int(region_bins["start"].min())
    end = int(region_bins["end"].max())
    return {
        "branch_s_version": str(version),
        "sample_id": str(sample_id),
        "sex_call": context["sex_call"],
        "predict_gender": context["predict_gender"],
        "sex_call_source": context["sex_call_source"],
        "region_class": region_class(chrom, is_par),
        "chrom": chrom,
        "is_par_region": int(is_par),
        "bin_count": int(len(region_bins)),
        "start": start,
        "end": end,
        "span_bp": int(max(end - start, 0)),
        "a_candidate_count": count_overlapping_candidates(region_bins, candidates),
        "mean_calibrated_z": mean_calibrated_z,
        "median_calibrated_z": finite_median(region_bins["calibrated_z"]),
        "mean_robust_z": finite_mean(region_bins["robust_z"]),
        "median_robust_z": finite_median(region_bins["robust_z"]),
        "abs_mean_calibrated_z": float(abs(mean_calibrated_z)) if np.isfinite(mean_calibrated_z) else np.nan,
        "par_bin_fraction": float(region_bins["is_par_region"].mean()),
        "branch_s_action": "SHADOW_ONLY",
        "final_report_impact": "none_shadow_only",
    }


def state_direction(state):
    text = str(state).strip().upper()
    if text.endswith("_GAIN"):
        return 1.0
    if text.endswith("_LOSS"):
        return -1.0
    return 0.0


def candidate_abs_zscore(rows):
    if rows is None or rows.empty:
        return np.nan
    if "a_abs_zscore" in rows.columns:
        values = pd.to_numeric(rows["a_abs_zscore"], errors="coerce").to_numpy(dtype=np.float64)
    elif "a_zscore" in rows.columns:
        values = np.abs(pd.to_numeric(rows["a_zscore"], errors="coerce").to_numpy(dtype=np.float64))
    else:
        values = np.array([], dtype=np.float64)
    values = values[np.isfinite(values)]
    return float(np.max(values)) if values.size else np.nan


def candidate_state_score(region_bins, candidates, state):
    overlaps = overlapping_candidate_rows(region_bins, candidates)
    if overlaps.empty or "state" not in overlaps.columns:
        return np.nan, ""
    target_state = "gain" if str(state).upper().endswith("_GAIN") else "loss"
    same = overlaps[overlaps["state"].astype(str).eq(target_state)].copy()
    opposite_state = "loss" if target_state == "gain" else "gain"
    opposite = overlaps[overlaps["state"].astype(str).eq(opposite_state)].copy()
    same_z = candidate_abs_zscore(same)
    opposite_z = candidate_abs_zscore(opposite)
    if np.isfinite(same_z) or np.isfinite(opposite_z):
        same_value = same_z if np.isfinite(same_z) else 0.0
        opposite_value = opposite_z if np.isfinite(opposite_z) else 0.0
        if same_value >= opposite_value:
            return float(same_value), "branch_a_candidate_zscore"
        return float(-opposite_value), "branch_a_candidate_zscore"
    if not same.empty or not opposite.empty:
        return (1.0 if len(same) >= len(opposite) else -1.0), "branch_a_candidate_overlap"
    return np.nan, ""


def make_state_score(sample_id, state, source_region_class, source_value, context, version, candidate_score=np.nan, candidate_reason=""):
    if not np.isfinite(source_value):
        score = np.nan
        status = "INSUFFICIENT_EVIDENCE"
        reason = "missing_nonpar_region"
    elif np.isfinite(candidate_score):
        score = float(candidate_score)
        status = "AVAILABLE"
        reason = str(candidate_reason or "branch_a_candidate_support")
    else:
        score = float(state_direction(state) * source_value)
        status = "AVAILABLE"
        reason = "nonpar_mean_calibrated_z"
    return {
        "branch_s_version": str(version),
        "sample_id": str(sample_id),
        "sex_call": context["sex_call"],
        "predict_gender": context["predict_gender"],
        "sca_state": state,
        "source_region_class": source_region_class,
        "state_score": score,
        "state_score_status": status,
        "state_score_reason": reason,
        "branch_s_action": "SHADOW_ONLY",
        "final_report_impact": "none_shadow_only",
    }


def build_branch_s_shadow(sample_id, bins, a_candidates, gender=None, version=BRANCH_S_VERSION):
    context = gender_context(gender)
    sex_bins = prepare_bins(bins)
    candidates = prepare_candidates(a_candidates)
    evidence_rows = []
    for chrom in ["chrX", "chrY"]:
        for is_par in [False, True]:
            region_bins = sex_bins[
                sex_bins["chrom"].astype(str).eq(chrom) & sex_bins["is_par_region"].astype(bool).eq(is_par)
            ].copy()
            if region_bins.empty:
                continue
            evidence_rows.append(summarize_region(sample_id, region_bins, candidates, context, version))
    evidence = pd.DataFrame(evidence_rows, columns=EVIDENCE_COLUMNS)

    means = {
        row["region_class"]: row["mean_calibrated_z"]
        for _, row in evidence.iterrows()
        if row.get("region_class") in {"X_NONPAR", "Y_NONPAR"}
    }
    region_bins_by_class = {}
    for chrom in ["chrX", "chrY"]:
        region_bins = sex_bins[
            sex_bins["chrom"].astype(str).eq(chrom) & sex_bins["is_par_region"].astype(bool).eq(False)
        ].copy()
        if not region_bins.empty:
            region_bins_by_class[f"{'X' if chrom == 'chrX' else 'Y'}_NONPAR"] = region_bins
    candidate_scores = {}
    for state, source_region in [
        ("X_GAIN", "X_NONPAR"),
        ("X_LOSS", "X_NONPAR"),
        ("Y_GAIN", "Y_NONPAR"),
        ("Y_LOSS", "Y_NONPAR"),
    ]:
        candidate_scores[state] = candidate_state_score(region_bins_by_class.get(source_region, pd.DataFrame()), candidates, state)
    scores = pd.DataFrame(
        [
            make_state_score(
                sample_id,
                "X_GAIN",
                "X_NONPAR",
                means.get("X_NONPAR", np.nan),
                context,
                version,
                *candidate_scores["X_GAIN"],
            ),
            make_state_score(
                sample_id,
                "X_LOSS",
                "X_NONPAR",
                means.get("X_NONPAR", np.nan),
                context,
                version,
                *candidate_scores["X_LOSS"],
            ),
            make_state_score(
                sample_id,
                "Y_GAIN",
                "Y_NONPAR",
                means.get("Y_NONPAR", np.nan),
                context,
                version,
                *candidate_scores["Y_GAIN"],
            ),
            make_state_score(
                sample_id,
                "Y_LOSS",
                "Y_NONPAR",
                means.get("Y_NONPAR", np.nan),
                context,
                version,
                *candidate_scores["Y_LOSS"],
            ),
        ],
        columns=SCORE_COLUMNS,
    )
    return evidence, scores


def summarize_branch_s_shadow(sample_id, evidence, scores, gender=None, version=BRANCH_S_VERSION):
    context = gender_context(gender)
    available_scores = (
        scores["state_score_status"].fillna("").astype(str).eq("AVAILABLE").sum()
        if scores is not None and "state_score_status" in scores.columns
        else 0
    )
    return {
        "sample_id": str(sample_id),
        "branch_s_version": str(version),
        "sex_call": context["sex_call"],
        "predict_gender": context["predict_gender"],
        "sex_call_source": context["sex_call_source"],
        "region_count": int(len(evidence) if evidence is not None else 0),
        "score_count": int(len(scores) if scores is not None else 0),
        "available_score_count": int(available_scores),
        "final_report_impact": "none_shadow_only",
        "replaces_current_sex_calling": False,
        "replaces_final_report": False,
    }


def main():
    args = parse_args()
    logger = setup_logger("branch_s_shadow", args.log or None)
    bins = read_table(args.input_bins, required_columns=["chrom", "start", "end"], empty_ok=False)
    candidates = read_optional_table(args.a_candidates)
    gender = read_optional_table(args.gender_tsv)
    evidence, scores = build_branch_s_shadow(
        sample_id=args.sample_id,
        bins=bins,
        a_candidates=candidates,
        gender=gender,
        version=args.version,
    )
    write_table(args.output_evidence, evidence)
    write_table(args.output_scores, scores)
    write_json(args.output_summary, summarize_branch_s_shadow(args.sample_id, evidence, scores, gender, args.version))
    logger.info(
        "wrote Branch S shadow evidence rows=%d score_rows=%d to %s",
        len(evidence),
        len(scores),
        args.output_evidence,
    )


if __name__ == "__main__":
    main()
