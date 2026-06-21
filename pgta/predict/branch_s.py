#!/biosoftware/miniconda/envs/snakemake_env/bin/python
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from pgta.core.logging import setup_logger
from pgta.predict.branch_b.common import read_table, write_json, write_table


BRANCH_S_VERSION = "branch_s_shadow_v1"
NONPAR_MEDIAN_CALIBRATED_MIN_ABS_Z = 2.0
NONPAR_MEDIAN_ROBUST_MIN_ABS_Z = 5.0
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


def nonpar_directional_support(region_bins, state):
    if region_bins is None or region_bins.empty:
        return "missing"
    direction = state_direction(state)
    if direction == 0.0:
        return "neutral"
    calibrated_median = finite_median(region_bins.get("calibrated_z", pd.Series(dtype=float)))
    robust_median = finite_median(region_bins.get("robust_z", pd.Series(dtype=float)))
    calibrated_signed = direction * calibrated_median if np.isfinite(calibrated_median) else np.nan
    robust_signed = direction * robust_median if np.isfinite(robust_median) else np.nan
    if (
        np.isfinite(calibrated_signed)
        and calibrated_signed >= NONPAR_MEDIAN_CALIBRATED_MIN_ABS_Z
    ) or (
        np.isfinite(robust_signed)
        and robust_signed >= NONPAR_MEDIAN_ROBUST_MIN_ABS_Z
    ):
        return "supported"
    if (
        np.isfinite(calibrated_signed)
        and calibrated_signed <= -NONPAR_MEDIAN_CALIBRATED_MIN_ABS_Z
    ) or (
        np.isfinite(robust_signed)
        and robust_signed <= -NONPAR_MEDIAN_ROBUST_MIN_ABS_Z
    ):
        return "opposite"
    return "neutral"


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


def sex_call_compatible_uncorroborated_review(context, state):
    sex = str((context or {}).get("sex_call", "") or "").strip().upper()
    state_text = str(state).strip().upper()
    return sex == "XX" and state_text == "X_LOSS"


def candidate_state_score(region_bins, candidates, state, context=None):
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
            support = nonpar_directional_support(region_bins, state)
            if support == "supported":
                return float(same_value), "branch_a_candidate_zscore_nonpar_corroborated"
            if support == "neutral" and sex_call_compatible_uncorroborated_review(context, state):
                return float(same_value), "branch_a_candidate_zscore_sex_call_compatible_uncorroborated_review"
            if support == "opposite":
                return float(-same_value), "branch_a_only_uncorroborated_by_nonpar_median"
            return 0.0, "branch_a_only_uncorroborated_by_nonpar_median"
        if nonpar_directional_support(region_bins, state) == "supported":
            return np.nan, ""
        return float(-opposite_value), "branch_a_candidate_zscore_opposite_direction"
    if not same.empty or not opposite.empty:
        support = nonpar_directional_support(region_bins, state)
        if len(same) >= len(opposite) and support == "supported":
            return 1.0, "branch_a_candidate_overlap_nonpar_corroborated"
        if len(same) >= len(opposite):
            return 0.0, "branch_a_only_uncorroborated_by_nonpar_median"
        if support == "supported":
            return np.nan, ""
        return -1.0, "branch_a_candidate_overlap_opposite_direction"
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
        reason = "nonpar_median_calibrated_z"
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


def expected_sex_ploidy(sex_call):
    sex = str(sex_call or "").strip().upper()
    if sex == "XX":
        return 2, 0
    if sex == "XY":
        return 1, 1
    return np.nan, np.nan


def _score_value(scores, state):
    if scores is None or scores.empty or "sca_state" not in scores.columns:
        return np.nan
    matched = scores[scores["sca_state"].astype(str).eq(str(state))]
    if matched.empty or "state_score" not in matched.columns:
        return np.nan
    return float(pd.to_numeric(matched["state_score"], errors="coerce").iloc[0])


def axis_direction_from_scores(scores, axis, min_abs_score=5.0):
    gain_score = _score_value(scores, f"{axis}_GAIN")
    loss_score = _score_value(scores, f"{axis}_LOSS")
    gain_value = gain_score if np.isfinite(gain_score) else -np.inf
    loss_value = loss_score if np.isfinite(loss_score) else -np.inf
    best = max(gain_value, loss_value)
    if not np.isfinite(best) or best < float(min_abs_score):
        return "neutral_or_uncertain"
    return "gain" if gain_value >= loss_value else "loss"


def dominant_sca_state(scores, min_score=5.0):
    if scores is None or scores.empty or "state_score" not in scores.columns:
        return "none_detected", np.nan
    frame = scores.copy()
    frame["state_score"] = pd.to_numeric(frame["state_score"], errors="coerce")
    frame = frame[frame["state_score"].notna()].copy()
    if frame.empty:
        return "none_detected", np.nan
    idx = frame["state_score"].idxmax()
    best_score = float(frame.loc[idx, "state_score"])
    if best_score < float(min_score):
        return "none_detected", best_score
    return str(frame.loc[idx, "sca_state"]), best_score


def region_context(evidence, region_class_name):
    if evidence is None or evidence.empty or "region_class" not in evidence.columns:
        return "not_available"
    return "available" if evidence["region_class"].astype(str).eq(region_class_name).any() else "not_available"


def branch_a_axis_support(evidence, axis):
    if evidence is None or evidence.empty or "region_class" not in evidence.columns or "a_candidate_count" not in evidence.columns:
        return "unknown"
    frame = evidence[evidence["region_class"].astype(str).str.startswith(f"{axis}_")].copy()
    if frame.empty:
        return "unknown"
    counts = pd.to_numeric(frame["a_candidate_count"], errors="coerce").fillna(0)
    return "present" if int(counts.sum()) > 0 else "absent"


def sca_confidence_tier(sca_state, score):
    if str(sca_state) == "none_detected" or not np.isfinite(score):
        return "SCA_NO_CALL"
    if score >= 30.0:
        return "SCA_REVIEW_STRONG"
    return "SCA_REVIEW_WEAK"


def sca_report_layer_class(sca_state, confidence_tier, scores, evidence):
    if str(sca_state) == "none_detected":
        if branch_a_axis_support(evidence, "X") == "present" or branch_a_axis_support(evidence, "Y") == "present":
            if has_uncorroborated_branch_a_support(scores):
                return "sca_filtered_or_sex_consistent_event"
        return "sca_no_call"
    if str(confidence_tier) == "SCA_REVIEW_STRONG":
        return "sca_report_review_event"
    return "sca_internal_review_event"


def sca_report_layer_reason(sca_state, confidence_tier, scores, evidence):
    if str(sca_state) == "none_detected":
        if has_uncorroborated_branch_a_support(scores):
            return "branch_a_only_uncorroborated_by_nonpar_median"
        return "insufficient_sca_evidence"
    score_reason = dominant_state_score_reason(scores, sca_state)
    if "sex_call_compatible_uncorroborated_review" in score_reason:
        return f"{str(confidence_tier).lower()}_with_sex_call_compatible_branch_a_support"
    return f"{str(confidence_tier).lower()}_with_nonpar_corroboration"


def dominant_state_score_reason(scores, sca_state):
    if scores is None or scores.empty or "sca_state" not in scores.columns:
        return ""
    matches = scores[scores["sca_state"].fillna("").astype(str).eq(str(sca_state))]
    if matches.empty or "state_score_reason" not in matches.columns:
        return ""
    return str(matches.iloc[0].get("state_score_reason", "") or "")


def has_uncorroborated_branch_a_support(scores):
    if scores is None or scores.empty or "state_score_reason" not in scores.columns:
        return False
    return scores["state_score_reason"].fillna("").astype(str).str.contains(
        "branch_a_only_uncorroborated", regex=False
    ).any()


def sca_uncertainty_reason(sca_state, context):
    reasons = ["locked_sca_truth_incomplete", "branch_s_not_final_validated"]
    if not str(context.get("sex_call", "")).strip():
        reasons.append("sex_call_missing")
    if str(sca_state) == "none_detected":
        reasons.append("insufficient_sca_evidence")
    return ";".join(reasons)


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
        row["region_class"]: row["median_calibrated_z"]
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
        candidate_scores[state] = candidate_state_score(
            region_bins_by_class.get(source_region, pd.DataFrame()),
            candidates,
            state,
            context,
        )
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
    expected_x_ploidy, expected_y_ploidy = expected_sex_ploidy(context["sex_call"])
    sca_state, sca_score = dominant_sca_state(scores)
    confidence_tier = sca_confidence_tier(sca_state, sca_score)
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
        "sample": str(sample_id),
        "expected_x_ploidy": expected_x_ploidy,
        "expected_y_ploidy": expected_y_ploidy,
        "x_nonpar_direction": axis_direction_from_scores(scores, "X"),
        "x_par_context": region_context(evidence, "X_PAR"),
        "y_nonpar_direction": axis_direction_from_scores(scores, "Y"),
        "y_par_or_homology_context": region_context(evidence, "Y_PAR"),
        "branch_a_x_support": branch_a_axis_support(evidence, "X"),
        "branch_a_y_support": branch_a_axis_support(evidence, "Y"),
        "sca_candidate_state": sca_state,
        "sca_confidence_tier": confidence_tier,
        "sca_report_layer_class": sca_report_layer_class(sca_state, confidence_tier, scores, evidence),
        "sca_report_layer_reason": sca_report_layer_reason(sca_state, confidence_tier, scores, evidence),
        "sca_output_mode": "review_development_only",
        "sca_uncertainty_reason": sca_uncertainty_reason(sca_state, context),
        "report_text_status": "development_only_not_final_reportable",
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
