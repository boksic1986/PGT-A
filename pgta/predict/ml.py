#!/biosoftware/miniconda/envs/snakemake_env/bin/python
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from pgta.core.logging import setup_logger


def parse_args():
    parser = argparse.ArgumentParser(description="Train an event-level candidate classifier with optional XGBoost/LightGBM and sklearn fallback.")
    parser.add_argument("--event-tsv", action="append", default=[], help="Repeat for multiple event TSVs")
    parser.add_argument("--labels-tsv", default="")
    parser.add_argument("--output-features", required=True)
    parser.add_argument("--output-cv-metrics", required=True)
    parser.add_argument("--output-calibration", required=True)
    parser.add_argument("--output-importance", required=True)
    parser.add_argument("--output-predictions", required=True)
    parser.add_argument("--output-summary", required=True)
    parser.add_argument("--backend", default="auto")
    parser.add_argument("--cv-folds", type=int, default=5)
    parser.add_argument("--random-seed", type=int, default=1)
    parser.add_argument("--log", default="")
    return parser.parse_args()


FEATURE_COLUMNS = [
    "n_bins",
    "segment_weight",
    "segment_mean_signal",
    "segment_median_signal",
    "segment_mean_robust_z",
    "segment_median_robust_z",
    "segment_abs_max_robust_z",
    "calibrated_mean_z",
    "calibrated_median_z",
    "empirical_pvalue",
    "empirical_qvalue",
    "priority_score",
    "chrom_fraction",
    "is_gain",
    "is_loss",
    "is_autosome",
    "is_sex_chromosome",
    "flag_count",
    "has_edge_flag",
    "has_weak_support_flag",
    "has_par_flag",
]


def ensure_parent(path_value):
    path = Path(path_value)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def load_event_tables(paths):
    frames = []
    for path_value in paths:
        path = Path(path_value)
        if not path.exists():
            continue
        df = pd.read_csv(path, sep="\t")
        if not df.empty:
            frames.append(df)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def build_feature_table(events_df):
    if events_df.empty:
        return pd.DataFrame(columns=["event_id", "sample_id", "chrom", "state"] + FEATURE_COLUMNS)
    feature_df = events_df.copy()
    feature_df["event_id"] = feature_df["event_id"].astype(str)
    feature_df["sample_id"] = feature_df["sample_id"].astype(str)
    feature_df["chrom"] = feature_df["chrom"].astype(str)
    feature_df["state"] = feature_df["state"].astype(str)
    feature_df["chrom_fraction"] = (
        (feature_df["end_bin"].astype(float) - feature_df["start_bin"].astype(float) + 1.0)
        / np.clip(feature_df["end_bin"].astype(float) + 1.0, a_min=1.0, a_max=None)
    )
    feature_df["is_gain"] = (feature_df["state"] == "gain").astype(int)
    feature_df["is_loss"] = (feature_df["state"] == "loss").astype(int)
    feature_df["is_sex_chromosome"] = feature_df["chrom"].isin(["chrX", "chrY"]).astype(int)
    feature_df["is_autosome"] = (1 - feature_df["is_sex_chromosome"]).astype(int)
    feature_df["artifact_flags"] = feature_df.get("artifact_flags", "").fillna("").astype(str)
    feature_df["flag_count"] = feature_df["artifact_flags"].map(lambda value: 0 if not value else len([item for item in value.split(",") if item]))
    feature_df["has_edge_flag"] = feature_df["artifact_flags"].str.contains("edge_event", regex=False).astype(int)
    feature_df["has_weak_support_flag"] = feature_df["artifact_flags"].str.contains("weak_empirical_support", regex=False).astype(int)
    feature_df["has_par_flag"] = feature_df["artifact_flags"].str.contains("par_overlap", regex=False).astype(int)
    for column in FEATURE_COLUMNS:
        if column not in feature_df.columns:
            feature_df[column] = np.nan
    return feature_df[["event_id", "sample_id", "chrom", "state"] + FEATURE_COLUMNS].copy()


def merge_labels(feature_df, labels_tsv):
    if not labels_tsv:
        return feature_df, False
    path = Path(labels_tsv)
    if not path.exists():
        return feature_df, False
    labels_df = pd.read_csv(path, sep="\t")
    if labels_df.empty:
        return feature_df, False
    if "event_id" in labels_df.columns:
        merged = feature_df.merge(labels_df[["event_id", "label"]], on="event_id", how="left")
        return merged, merged["label"].notna().any()
    required = {"sample_id", "chrom", "state", "label"}
    if required.issubset(labels_df.columns) and {"chrom", "state"}.issubset(feature_df.columns):
        merged = feature_df.merge(labels_df[list(required)], on=["sample_id", "chrom", "state"], how="left")
        return merged, merged["label"].notna().any()
    return feature_df, False


def select_backend(name):
    backend = str(name).strip().lower()
    errors = []
    if backend in {"auto", "xgboost"}:
        try:
            from xgboost import XGBClassifier  # type: ignore

            return "xgboost", XGBClassifier(
                n_estimators=200,
                max_depth=4,
                learning_rate=0.05,
                subsample=0.8,
                colsample_bytree=0.8,
                objective="binary:logistic",
                eval_metric="logloss",
                random_state=1,
            )
        except Exception as exc:  # pragma: no cover
            errors.append(f"xgboost unavailable: {exc}")
    if backend in {"auto", "lightgbm"}:
        try:
            from lightgbm import LGBMClassifier  # type: ignore

            return "lightgbm", LGBMClassifier(
                n_estimators=200,
                learning_rate=0.05,
                num_leaves=31,
                subsample=0.8,
                colsample_bytree=0.8,
                random_state=1,
            )
        except Exception as exc:  # pragma: no cover
            errors.append(f"lightgbm unavailable: {exc}")
    try:
        from sklearn.ensemble import GradientBoostingClassifier  # type: ignore

        return "sklearn", GradientBoostingClassifier(random_state=1)
    except Exception as exc:  # pragma: no cover
        errors.append(f"sklearn unavailable: {exc}")
    raise RuntimeError("; ".join(errors))


def compute_feature_importance(model, feature_names):
    if hasattr(model, "feature_importances_"):
        return pd.DataFrame({"feature": feature_names, "importance": model.feature_importances_}).sort_values(
            "importance", ascending=False
        )
    return pd.DataFrame({"feature": feature_names, "importance": np.nan})


def main():
    args = parse_args()
    logger = setup_logger("cnv_ml", args.log or None)
    events_df = load_event_tables(args.event_tsv)
    feature_df = build_feature_table(events_df)
    ensure_parent(args.output_features)
    feature_df.to_csv(args.output_features, sep="\t", index=False)

    empty_cv = pd.DataFrame(columns=["fold", "backend", "accuracy", "precision", "recall", "f1", "roc_auc", "average_precision", "brier"])
    empty_cal = pd.DataFrame(columns=["bin_left", "bin_right", "mean_predicted", "observed_positive_rate", "count"])
    empty_imp = pd.DataFrame(columns=["feature", "importance"])
    empty_pred = pd.DataFrame(columns=["event_id", "sample_id", "y_true", "y_score", "y_pred"])

    feature_with_labels, has_labels = merge_labels(feature_df, args.labels_tsv)
    if not has_labels or feature_with_labels["label"].dropna().nunique() < 2:
        empty_cv.to_csv(args.output_cv_metrics, sep="\t", index=False)
        empty_cal.to_csv(args.output_calibration, sep="\t", index=False)
        empty_imp.to_csv(args.output_importance, sep="\t", index=False)
        empty_pred.to_csv(args.output_predictions, sep="\t", index=False)
        summary = {
            "status": "skipped",
            "reason": "labels_unavailable_or_single_class",
            "event_count": int(len(feature_df)),
        }
        ensure_parent(args.output_summary).write_text(json.dumps(summary, indent=2), encoding="utf-8")
        logger.info("ml stage skipped: %s", summary["reason"])
        return

    try:
        from sklearn.calibration import calibration_curve  # type: ignore
        from sklearn.metrics import accuracy_score, average_precision_score, brier_score_loss, f1_score, precision_score, recall_score, roc_auc_score  # type: ignore
        from sklearn.model_selection import StratifiedKFold  # type: ignore
    except Exception as exc:  # pragma: no cover
        empty_cv.to_csv(args.output_cv_metrics, sep="\t", index=False)
        empty_cal.to_csv(args.output_calibration, sep="\t", index=False)
        empty_imp.to_csv(args.output_importance, sep="\t", index=False)
        empty_pred.to_csv(args.output_predictions, sep="\t", index=False)
        summary = {
            "status": "skipped",
            "reason": f"sklearn_unavailable: {exc}",
            "event_count": int(len(feature_df)),
        }
        ensure_parent(args.output_summary).write_text(json.dumps(summary, indent=2), encoding="utf-8")
        logger.info("ml stage skipped: %s", summary["reason"])
        return

    modeling_df = feature_with_labels.dropna(subset=["label"]).copy()
    X = modeling_df[FEATURE_COLUMNS].fillna(0.0).to_numpy(dtype=np.float64)
    y = modeling_df["label"].astype(int).to_numpy(dtype=np.int64)

    backend_name, model = select_backend(args.backend)
    splitter = StratifiedKFold(n_splits=max(2, args.cv_folds), shuffle=True, random_state=args.random_seed)

    cv_rows = []
    prediction_rows = []
    for fold_index, (train_idx, test_idx) in enumerate(splitter.split(X, y), start=1):
        fold_model = select_backend(backend_name)[1]
        fold_model.fit(X[train_idx], y[train_idx])
        if hasattr(fold_model, "predict_proba"):
            y_score = fold_model.predict_proba(X[test_idx])[:, 1]
        else:
            y_score = fold_model.decision_function(X[test_idx])  # pragma: no cover
            y_score = 1.0 / (1.0 + np.exp(-y_score))
        y_pred = (y_score >= 0.5).astype(int)
        cv_rows.append(
            {
                "fold": fold_index,
                "backend": backend_name,
                "accuracy": accuracy_score(y[test_idx], y_pred),
                "precision": precision_score(y[test_idx], y_pred, zero_division=0),
                "recall": recall_score(y[test_idx], y_pred, zero_division=0),
                "f1": f1_score(y[test_idx], y_pred, zero_division=0),
                "roc_auc": roc_auc_score(y[test_idx], y_score),
                "average_precision": average_precision_score(y[test_idx], y_score),
                "brier": brier_score_loss(y[test_idx], y_score),
            }
        )
        prediction_rows.extend(
            {
                "event_id": modeling_df.iloc[test_index]["event_id"],
                "sample_id": modeling_df.iloc[test_index]["sample_id"],
                "y_true": int(y[test_pos]),
                "y_score": float(y_score[test_pos]),
                "y_pred": int(y_pred[test_pos]),
            }
            for test_pos, test_index in enumerate(test_idx)
        )

    model.fit(X, y)
    importance_df = compute_feature_importance(model, FEATURE_COLUMNS)
    prediction_df = pd.DataFrame(prediction_rows)
    calib_true, calib_pred = calibration_curve(prediction_df["y_true"], prediction_df["y_score"], n_bins=min(10, len(prediction_df)))
    calibration_df = pd.DataFrame(
        {
            "bin_left": np.linspace(0.0, 1.0, num=len(calib_true), endpoint=False),
            "bin_right": np.linspace(1.0 / max(len(calib_true), 1), 1.0, num=len(calib_true)),
            "mean_predicted": calib_pred,
            "observed_positive_rate": calib_true,
            "count": np.nan,
        }
    )
    cv_df = pd.DataFrame(cv_rows)

    ensure_parent(args.output_cv_metrics)
    cv_df.to_csv(args.output_cv_metrics, sep="\t", index=False)
    calibration_df.to_csv(args.output_calibration, sep="\t", index=False)
    importance_df.to_csv(args.output_importance, sep="\t", index=False)
    prediction_df.to_csv(args.output_predictions, sep="\t", index=False)
    summary = {
        "status": "completed",
        "backend": backend_name,
        "event_count": int(len(modeling_df)),
        "positive_event_count": int(y.sum()),
        "cv_folds": int(args.cv_folds),
        "mean_metrics": cv_df.mean(numeric_only=True).to_dict(),
    }
    ensure_parent(args.output_summary).write_text(json.dumps(summary, indent=2), encoding="utf-8")
    logger.info("ml stage completed: backend=%s events=%d", backend_name, len(modeling_df))


if __name__ == "__main__":
    main()
