#!/biosoftware/miniconda/envs/snakemake_env/bin/python
import argparse

import numpy as np
import pandas as pd

from pgta.predict.branch_b.common import (
    CandidateEvent,
    STATE_TO_SVTYPE,
    annotate_region_risk,
    effective_sample_size,
    load_sample_bins,
    robust_z,
    weighted_mean,
    write_json,
    write_table,
)
from pgta.core.logging import setup_logger


def parse_args():
    parser = argparse.ArgumentParser(description="Branch B CNV calling with weighted CBS + HMM refinement.")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--npz", default="")
    parser.add_argument("--input-bins", default="")
    parser.add_argument("--output-bins", required=True)
    parser.add_argument("--output-candidates", required=True)
    parser.add_argument("--output-summary", required=True)
    parser.add_argument("--branch", default="B")
    parser.add_argument("--correction-model", default="2d_loess_gc_mappability")
    parser.add_argument("--min-bins", type=int, default=5)
    parser.add_argument("--max-segments-per-chrom", type=int, default=12)
    parser.add_argument("--split-threshold", type=float, default=2.5)
    parser.add_argument("--hmm-state-shift", type=float, default=2.5)
    parser.add_argument("--hmm-stay-prob", type=float, default=0.995)
    parser.add_argument("--min-event-bins", type=int, default=3)
    parser.add_argument("--min-event-z", type=float, default=1.5)
    parser.add_argument("--chromosome-z-threshold", type=float, default=2.5)
    parser.add_argument("--chromosome-min-effective-bins", type=float, default=10.0)
    parser.add_argument("--chromosome-min-clean-fraction", type=float, default=0.30)
    parser.add_argument("--masked-gap-rescue-min-abs-local-z", type=float, default=1.8)
    parser.add_argument("--masked-gap-rescue-min-median-z", type=float, default=4.0)
    parser.add_argument("--masked-gap-rescue-min-chrom-shift-z", type=float, default=0.75)
    parser.add_argument("--fusion-iou-threshold", type=float, default=0.8)
    parser.add_argument("--log", default="")
    return parser.parse_args()


def load_calling_bins(args):
    if args.input_bins:
        bins_df = pd.read_csv(args.input_bins, sep="\t")
        if bins_df.empty:
            raise ValueError(f"input bins table is empty: {args.input_bins}")
        binsize = int(bins_df["end"].sub(bins_df["start"]).median()) if {"start", "end"}.issubset(bins_df.columns) else 0
        quality = float("nan")
        signal_column = "signal_for_calling" if "signal_for_calling" in bins_df.columns else "normalized_signal"
        if signal_column not in bins_df.columns:
            raise ValueError(f"input bins table missing signal column: {args.input_bins}")
        return bins_df, binsize, quality, signal_column
    if args.npz:
        bins_df, binsize, quality = load_sample_bins(args.npz)
        return bins_df, binsize, quality, "normalized_signal"
    raise ValueError("Either --npz or --input-bins is required")


def robust_center_scale(values):
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return 0.0, 1.0
    center = float(np.median(values))
    mad = float(np.median(np.abs(values - center)))
    scale = max(1.4826 * mad, float(np.std(values)), 1e-6)
    return center, scale


def select_background_signal(bins_df, signal_column):
    selector_candidates = [
        bins_df["is_autosome"].fillna(0).astype(int).eq(1) & bins_df["calling_seed_eligible"].fillna(1).astype(int).eq(1),
        bins_df["is_autosome"].fillna(0).astype(int).eq(1),
        bins_df["calling_seed_eligible"].fillna(1).astype(int).eq(1),
        pd.Series(True, index=bins_df.index),
    ]
    for selector in selector_candidates:
        values = pd.to_numeric(bins_df.loc[selector, signal_column], errors="coerce").to_numpy(dtype=np.float64)
        values = values[np.isfinite(values)]
        if values.size:
            return values
    return np.zeros(1, dtype=np.float64)


def annotate_calling_z_scores(bins_df, signal_column):
    frame = bins_df.copy()
    background_center, background_scale = robust_center_scale(select_background_signal(frame, signal_column))
    frame["global_robust_z"] = 0.0
    chrom_signal_means = {}
    chrom_is_autosome = {}
    for chrom, chrom_df in frame.groupby("chrom", sort=False):
        chrom_key = str(chrom)
        chrom_signal_means[chrom_key] = weighted_mean(
            pd.to_numeric(chrom_df[signal_column], errors="coerce").to_numpy(dtype=np.float64),
            chrom_df["calling_weight"].to_numpy(dtype=np.float64),
        )
        chrom_is_autosome[chrom_key] = bool(chrom_df["is_autosome"].fillna(0).astype(int).iloc[0])

    chrom_frames = []
    chrom_shift_summary = {}
    for chrom, chrom_df in frame.groupby("chrom", sort=False):
        chrom_df = chrom_df.copy()
        background_source = frame[frame["chrom"].astype(str) != str(chrom)].copy()
        if background_source.empty:
            background_source = frame
        chrom_background_center, chrom_background_scale = robust_center_scale(
            select_background_signal(background_source, signal_column)
        )
        chrom_key = str(chrom)
        chrom_df["global_robust_z"] = (
            pd.to_numeric(chrom_df[signal_column], errors="coerce").to_numpy(dtype=np.float64) - chrom_background_center
        ) / chrom_background_scale
        local_robust_z = robust_z(chrom_df[signal_column].to_numpy(dtype=np.float64))
        eligible = chrom_df["calling_seed_eligible"].fillna(1).astype(int).eq(1)
        weight_values = chrom_df["calling_weight"].to_numpy(dtype=np.float64)
        chrom_background_means = [
            mean_value
            for other_chrom, mean_value in chrom_signal_means.items()
            if other_chrom != chrom_key and np.isfinite(mean_value)
            and ((chrom_is_autosome.get(chrom_key, False) and chrom_is_autosome.get(other_chrom, False)) or not chrom_is_autosome.get(chrom_key, False))
        ]
        if not chrom_background_means:
            chrom_background_means = [
                mean_value for other_chrom, mean_value in chrom_signal_means.items() if other_chrom != chrom_key and np.isfinite(mean_value)
            ]
        shift_center, shift_scale = robust_center_scale(chrom_background_means)
        shift_scale = max(float(shift_scale), float(background_scale), 1e-6)
        chrom_mean_signal = chrom_signal_means.get(chrom_key, np.nan)
        shift_z = (chrom_mean_signal - shift_center) / shift_scale if np.isfinite(chrom_mean_signal) else np.nan
        shift_z = 0.0 if not np.isfinite(shift_z) else float(shift_z)
        chrom_df["raw_local_robust_z"] = local_robust_z
        chrom_df["chrom_shift_z"] = shift_z
        chrom_df["raw_robust_z"] = chrom_df["raw_local_robust_z"] + chrom_df["chrom_shift_z"]
        chrom_df["robust_z"] = chrom_df["raw_robust_z"] / chrom_df["variance_inflation"].to_numpy(dtype=np.float64)
        chrom_df.loc[chrom_df["calling_seed_eligible"] == 0, "calling_weight"] *= 0.5
        chrom_frames.append(chrom_df)
        chrom_shift_summary[str(chrom)] = shift_z

    return (
        pd.concat(chrom_frames, ignore_index=True),
        {
            "background_center": background_center,
            "background_scale": background_scale,
            "chrom_shift_z": chrom_shift_summary,
        },
    )


def build_chromosome_dosage_event(args, chrom_df):
    weights = np.clip(chrom_df["calling_weight"].to_numpy(dtype=np.float64), 1e-6, None)
    chrom_shift_z = float(chrom_df["chrom_shift_z"].iloc[0]) if "chrom_shift_z" in chrom_df.columns else np.nan
    if not np.isfinite(chrom_shift_z) or abs(chrom_shift_z) < float(args.chromosome_z_threshold):
        return []

    effective_bins = effective_sample_size(weights)
    if effective_bins < float(args.chromosome_min_effective_bins):
        return []

    clean_fraction = float(chrom_df["region_risk_class"].eq("clean").mean()) if "region_risk_class" in chrom_df.columns else 1.0
    if clean_fraction < float(args.chromosome_min_clean_fraction):
        return []

    state = "gain" if chrom_shift_z > 0 else "loss"
    mean_robust_z = float(np.average(chrom_df["robust_z"], weights=weights))
    median_robust_z = float(np.median(chrom_df["robust_z"]))
    signal_values = pd.to_numeric(chrom_df["normalized_signal"], errors="coerce").to_numpy(dtype=np.float64)
    chrom = str(chrom_df["chrom"].iloc[0])
    event = CandidateEvent(
        event_id=f"{args.sample_id}.{chrom}.chr_dosage.{state}",
        sample_id=args.sample_id,
        branch=args.branch,
        correction_model=args.correction_model,
        caller="chromosome_dosage_detector",
        caller_stage="chromosome_dosage",
        chrom=chrom,
        start=int(chrom_df["start"].iloc[0]),
        end=int(chrom_df["end"].iloc[-1]),
        start_bin=int(chrom_df["bin_index"].iloc[0]),
        end_bin=int(chrom_df["bin_index"].iloc[-1]),
        n_bins=int(len(chrom_df)),
        state=state,
        svtype=STATE_TO_SVTYPE[state],
        segment_id=f"{chrom}_chr_dosage",
        segment_weight=float(weights.sum()),
        segment_mean_signal=float(weighted_mean(signal_values, weights)),
        segment_median_signal=float(np.median(signal_values)),
        segment_mean_robust_z=mean_robust_z,
        segment_median_robust_z=median_robust_z,
        segment_abs_max_robust_z=float(np.max(np.abs(chrom_df["robust_z"].to_numpy(dtype=np.float64)))),
    )
    return [event.to_dict()]


def iter_boolean_blocks(mask_values, min_block_bins):
    blocks = []
    start = None
    for index, is_selected in enumerate(mask_values):
        if is_selected and start is None:
            start = index
            continue
        if (not is_selected) and start is not None:
            if (index - start) >= int(min_block_bins):
                blocks.append((start, index))
            start = None
    if start is not None and (len(mask_values) - start) >= int(min_block_bins):
        blocks.append((start, len(mask_values)))
    return blocks


def build_masked_gap_rescue_events(args, chrom_df):
    if chrom_df.empty:
        return []
    if "mask_label" not in chrom_df.columns or "gap_centromere_telomere_overlap_fraction" not in chrom_df.columns:
        return []
    chrom_shift_z = float(chrom_df["chrom_shift_z"].iloc[0]) if "chrom_shift_z" in chrom_df.columns else np.nan
    if not np.isfinite(chrom_shift_z) or abs(chrom_shift_z) < float(args.masked_gap_rescue_min_chrom_shift_z):
        return []

    rescue_state = "loss" if chrom_shift_z > 0 else "gain"
    rescue_direction = -1.0 if rescue_state == "loss" else 1.0
    normalized_local_z = robust_z(pd.to_numeric(chrom_df["normalized_signal"], errors="coerce").to_numpy(dtype=np.float64))
    hard_gap_mask = (
        chrom_df["mask_label"].fillna("").astype(str).eq("hard").to_numpy(dtype=bool)
        & (pd.to_numeric(chrom_df["gap_centromere_telomere_overlap_fraction"], errors="coerce").fillna(0.0).to_numpy(dtype=np.float64) > 0.0)
    )
    rescue_seed_mask = hard_gap_mask & ((rescue_direction * normalized_local_z) >= float(args.masked_gap_rescue_min_abs_local_z))
    rescue_blocks = iter_boolean_blocks(rescue_seed_mask, min_block_bins=args.min_event_bins)

    rescue_events = []
    for block_index, (left, right) in enumerate(rescue_blocks, start=1):
        segment = chrom_df.iloc[left:right].copy()
        local_z = normalized_local_z[left:right]
        if local_z.size < int(args.min_event_bins):
            continue
        if np.median(np.abs(local_z)) < float(args.masked_gap_rescue_min_median_z):
            continue
        weights = np.clip(segment["calling_weight"].to_numpy(dtype=np.float64), 1e-6, None)
        event = CandidateEvent(
            event_id=f"{args.sample_id}.{segment['chrom'].iloc[0]}.masked_gap_rescue_{block_index}.{rescue_state}",
            sample_id=args.sample_id,
            branch=args.branch,
            correction_model=args.correction_model,
            caller="masked_gap_rescue_detector",
            caller_stage="masked_gap_rescue",
            chrom=str(segment["chrom"].iloc[0]),
            start=int(segment["start"].iloc[0]),
            end=int(segment["end"].iloc[-1]),
            start_bin=int(segment["bin_index"].iloc[0]),
            end_bin=int(segment["bin_index"].iloc[-1]),
            n_bins=int(len(segment)),
            state=rescue_state,
            svtype=STATE_TO_SVTYPE[rescue_state],
            segment_id=f"{segment['chrom'].iloc[0]}_masked_gap_rescue_{block_index}",
            segment_weight=float(weights.sum()),
            segment_mean_signal=float(weighted_mean(segment["normalized_signal"], weights)),
            segment_median_signal=float(np.median(segment["normalized_signal"])),
            segment_mean_robust_z=float(np.average(local_z, weights=weights)),
            segment_median_robust_z=float(np.median(local_z)),
            segment_abs_max_robust_z=float(np.max(np.abs(local_z))),
        )
        rescue_events.append(event.to_dict())
    return rescue_events


def event_interval_iou(left_event, right_event):
    overlap_left = max(int(left_event.get("start_bin", -1)), int(right_event.get("start_bin", -1)))
    overlap_right = min(int(left_event.get("end_bin", -1)), int(right_event.get("end_bin", -1)))
    overlap = max(0, overlap_right - overlap_left + 1)
    left_length = max(1, int(left_event.get("end_bin", -1)) - int(left_event.get("start_bin", -1)) + 1)
    right_length = max(1, int(right_event.get("end_bin", -1)) - int(right_event.get("start_bin", -1)) + 1)
    union = max(left_length + right_length - overlap, 1)
    return float(overlap) / float(union)


def event_priority(event):
    caller = str(event.get("caller", ""))
    mean_z = abs(float(event.get("segment_mean_robust_z", 0.0) or 0.0))
    median_z = abs(float(event.get("segment_median_robust_z", 0.0) or 0.0))
    max_z = abs(float(event.get("segment_abs_max_robust_z", 0.0) or 0.0))
    support = max(mean_z, median_z, max_z)
    if caller == "chromosome_dosage_detector":
        return (2, support, int(event.get("n_bins", 0)))
    if caller == "segment_level_detector":
        return (1, support, -int(event.get("n_bins", 0)))
    return (0, support, -int(event.get("n_bins", 0)))


def fuse_candidate_events(segment_events, chromosome_events, iou_threshold=0.8):
    fused = []
    for event in list(chromosome_events) + list(segment_events):
        merged = False
        for index, existing in enumerate(fused):
            same_axis = (
                str(existing.get("chrom", "")) == str(event.get("chrom", ""))
                and str(existing.get("state", "")) == str(event.get("state", ""))
            )
            if not same_axis:
                continue
            if event_interval_iou(existing, event) < float(iou_threshold):
                continue
            if event_priority(event) > event_priority(existing):
                fused[index] = event
            merged = True
            break
        if not merged:
            fused.append(event)
    return fused


def segment_score(values, weights, split_index):
    left_weights = weights[:split_index]
    right_weights = weights[split_index:]
    if left_weights.sum() <= 0 or right_weights.sum() <= 0:
        return 0.0
    left_mean = np.average(values[:split_index], weights=left_weights)
    right_mean = np.average(values[split_index:], weights=right_weights)
    penalty = np.sqrt((left_weights.sum() * right_weights.sum()) / (left_weights.sum() + right_weights.sum()))
    return abs(left_mean - right_mean) * penalty


def recursive_segment(values, weights, min_bins, threshold, max_segments):
    segments = [(0, len(values))]
    pending = [(0, values, weights)]
    while pending and len(segments) < max_segments:
        origin, current_values, current_weights = pending.pop(0)
        if len(current_values) < (2 * min_bins):
            continue
        best_score = 0.0
        best_split = None
        for split in range(min_bins, len(current_values) - min_bins + 1):
            score = segment_score(current_values, current_weights, split)
            if score > best_score:
                best_score = score
                best_split = split
        if best_split is None or best_score < threshold:
            continue
        left = (origin, origin + best_split)
        right = (origin + best_split, origin + len(current_values))
        segments = [segment for segment in segments if segment != (origin, origin + len(current_values))]
        segments.extend([left, right])
        segments.sort()
        pending.append((left[0], current_values[:best_split], current_weights[:best_split]))
        pending.append((right[0], current_values[best_split:], current_weights[best_split:]))
    return sorted(segments)

def gaussian_logpdf(value, mean, sigma):
    sigma = max(float(sigma), 1e-6)
    return -0.5 * np.log(2.0 * np.pi * sigma * sigma) - ((value - mean) ** 2) / (2.0 * sigma * sigma)


def viterbi(segment_values, sigma, state_shift, stay_prob):
    states = ["loss", "neutral", "gain"]
    means = {"loss": -state_shift, "neutral": 0.0, "gain": state_shift}
    switch_prob = (1.0 - stay_prob) / 2.0
    transition = {
        source: {
            target: np.log(stay_prob if source == target else switch_prob)
            for target in states
        }
        for source in states
    }
    initial = {state: np.log(1.0 / len(states)) for state in states}
    scores = []
    paths = []
    for index, value in enumerate(segment_values):
        current_scores = {}
        current_paths = {}
        for state in states:
            emission = gaussian_logpdf(value, means[state], sigma)
            if index == 0:
                current_scores[state] = initial[state] + emission
                current_paths[state] = [state]
                continue
            best_prev = None
            best_score = None
            for prev_state in states:
                score = scores[-1][prev_state] + transition[prev_state][state] + emission
                if best_score is None or score > best_score:
                    best_score = score
                    best_prev = prev_state
            current_scores[state] = best_score
            current_paths[state] = paths[-1][best_prev] + [state]
        scores.append(current_scores)
        paths.append(current_paths)
    best_state = max(scores[-1], key=scores[-1].get)
    return paths[-1][best_state]


def infer_signal_state(mean_z, median_z, min_event_z):
    signed_support = float(mean_z) if abs(float(mean_z)) >= abs(float(median_z)) else float(median_z)
    if abs(signed_support) < float(min_event_z):
        return "neutral"
    return "gain" if signed_support > 0 else "loss"


def reconcile_segment_state(hmm_state, mean_z, median_z, min_event_z):
    signal_state = infer_signal_state(mean_z, median_z, min_event_z)
    if signal_state == "neutral":
        return "neutral"
    if hmm_state == "neutral":
        return signal_state
    if hmm_state != signal_state:
        return signal_state
    return hmm_state


def merge_adjacent_segment_calls(segment_calls):
    merged = []
    for call in segment_calls:
        if not merged:
            merged.append(dict(call))
            continue
        previous = merged[-1]
        same_state = str(previous.get("state", "")) == str(call.get("state", ""))
        same_block = int(previous.get("block_index", -1)) == int(call.get("block_index", -1))
        touching = int(previous.get("right", -1)) == int(call.get("left", -2))
        if same_state and same_block and touching:
            previous_width = max(int(previous["right"]) - int(previous["left"]), 1)
            current_width = max(int(call["right"]) - int(call["left"]), 1)
            previous["right"] = int(call["right"])
            previous["mean_z"] = float(
                np.average(
                    [float(previous["mean_z"]), float(call["mean_z"])],
                    weights=[previous_width, current_width],
                )
            )
            previous["median_z"] = float(
                np.average(
                    [float(previous["median_z"]), float(call["median_z"])],
                    weights=[previous_width, current_width],
                )
            )
            previous["source_segment_ids"] = list(previous.get("source_segment_ids", [])) + list(call.get("source_segment_ids", []))
            continue
        merged.append(dict(call))
    for merge_index, call in enumerate(merged, start=1):
        chrom = str(call.get("chrom", ""))
        block_index = int(call.get("block_index", merge_index))
        call["segment_id"] = f"{chrom}_block{block_index}_merged{merge_index}"
    return merged


def iter_seed_blocks(chrom_df, min_block_bins):
    eligible = chrom_df["calling_seed_eligible"].fillna(1).astype(int).eq(1).to_numpy(dtype=bool)
    blocks = []
    start = None
    for index, is_eligible in enumerate(eligible):
        if is_eligible and start is None:
            start = index
            continue
        if (not is_eligible) and start is not None:
            if (index - start) >= int(min_block_bins):
                blocks.append((start, index))
            start = None
    if start is not None and (len(eligible) - start) >= int(min_block_bins):
        blocks.append((start, len(eligible)))
    return blocks


def build_segment_level_events(args, chrom_df, logger):
    segment_rows = []
    seed_blocks = iter_seed_blocks(chrom_df, min_block_bins=max(args.min_bins, args.min_event_bins))
    for block_index, (block_left, block_right) in enumerate(seed_blocks, start=1):
        block = chrom_df.iloc[block_left:block_right].copy()
        values = block["robust_z"].to_numpy(dtype=np.float64)
        weights = block["calling_weight"].to_numpy(dtype=np.float64)
        block_segments = recursive_segment(values, weights, args.min_bins, args.split_threshold, args.max_segments_per_chrom)
        for segment_index, (left, right) in enumerate(block_segments, start=1):
            global_left = block_left + left
            global_right = block_left + right
            segment = chrom_df.iloc[global_left:global_right].copy()
            segment_id = f"{chrom_df['chrom'].iloc[0]}_block{block_index}_seg{segment_index}"
            mean_z = float(np.average(segment["robust_z"], weights=np.clip(segment["calling_weight"], 1e-6, None)))
            median_z = float(np.median(segment["robust_z"]))
            chrom_df.loc[segment.index, "segment_id"] = segment_id
            chrom_df.loc[segment.index, "segment_mean_robust_z"] = mean_z
            segment_rows.append(
                {
                    "chrom": str(chrom_df["chrom"].iloc[0]),
                    "block_index": int(block_index),
                    "segment_id": segment_id,
                    "left": global_left,
                    "right": global_right,
                    "mean_z": mean_z,
                    "median_z": median_z,
                }
            )
    if not segment_rows:
        logger.info("chrom=%s seed_blocks=0 segment_level_events=0", chrom_df["chrom"].iloc[0])
        return chrom_df, [], []
    sigma = max(float(np.std([row["mean_z"] for row in segment_rows])), 1.0)
    hmm_states = viterbi([row["mean_z"] for row in segment_rows], sigma, args.hmm_state_shift, args.hmm_stay_prob)

    segment_calls = []
    for hmm_state, meta in zip(hmm_states, segment_rows):
        segment = chrom_df.iloc[meta["left"]:meta["right"]].copy()
        state = reconcile_segment_state(hmm_state, meta["mean_z"], meta["median_z"], args.min_event_z)
        chrom_df.loc[segment.index, "hmm_state"] = state
        if state == "neutral" or len(segment) < args.min_event_bins:
            continue
        segment_calls.append(
            {
                "chrom": str(segment["chrom"].iloc[0]),
                "block_index": int(meta["block_index"]),
                "left": int(meta["left"]),
                "right": int(meta["right"]),
                "state": state,
                "mean_z": float(meta["mean_z"]),
                "median_z": float(meta["median_z"]),
                "source_segment_ids": [str(meta["segment_id"])],
            }
        )

    merged_segment_calls = merge_adjacent_segment_calls(segment_calls)
    segment_events = []
    for meta in merged_segment_calls:
        segment = chrom_df.iloc[int(meta["left"]):int(meta["right"])].copy()
        median_z = float(np.median(segment["robust_z"]))
        event = CandidateEvent(
            event_id=f"{args.sample_id}.{meta['segment_id']}.{meta['state']}",
            sample_id=args.sample_id,
            branch=args.branch,
            correction_model=args.correction_model,
            caller="segment_level_detector",
            caller_stage="segment_level",
            chrom=str(segment["chrom"].iloc[0]),
            start=int(segment["start"].iloc[0]),
            end=int(segment["end"].iloc[-1]),
            start_bin=int(segment["bin_index"].iloc[0]),
            end_bin=int(segment["bin_index"].iloc[-1]),
            n_bins=int(len(segment)),
            state=str(meta["state"]),
            svtype=STATE_TO_SVTYPE[str(meta["state"])],
            segment_id=str(meta["segment_id"]),
            segment_weight=float(segment["calling_weight"].sum()),
            segment_mean_signal=float(weighted_mean(segment["normalized_signal"], segment["calling_weight"])),
            segment_median_signal=float(np.median(segment["normalized_signal"])),
            segment_mean_robust_z=float(np.average(segment["robust_z"], weights=np.clip(segment["calling_weight"], 1e-6, None))),
            segment_median_robust_z=median_z,
            segment_abs_max_robust_z=float(np.max(np.abs(segment["robust_z"]))),
        )
        segment_events.append(event.to_dict())
    logger.info(
        "chrom=%s seed_blocks=%d segments=%d segment_level_events=%d",
        chrom_df["chrom"].iloc[0],
        len(seed_blocks),
        len(segment_rows),
        len(segment_events),
    )
    return chrom_df, segment_events, segment_rows


def build_candidate_events(args, chrom_df, logger):
    chrom_df, segment_events, segment_rows = build_segment_level_events(args, chrom_df, logger)
    chromosome_events = build_chromosome_dosage_event(args, chrom_df)
    masked_gap_rescue_events = build_masked_gap_rescue_events(args, chrom_df)
    events = fuse_candidate_events(
        segment_events + masked_gap_rescue_events,
        chromosome_events,
        iou_threshold=args.fusion_iou_threshold,
    )
    logger.info(
        "chrom=%s segments=%d segment_events=%d masked_gap_rescue_events=%d chromosome_events=%d candidate_events=%d",
        chrom_df["chrom"].iloc[0],
        len(segment_rows),
        len(segment_events),
        len(masked_gap_rescue_events),
        len(chromosome_events),
        len(events),
    )
    return chrom_df, events


def main():
    args = parse_args()
    logger = setup_logger("cnv_calling", args.log or None)
    bins_df, binsize, quality, signal_column = load_calling_bins(args)
    bins_df = annotate_region_risk(bins_df)
    bins_df["robust_z"] = 0.0
    bins_df["segment_id"] = ""
    bins_df["segment_mean_robust_z"] = 0.0
    bins_df["hmm_state"] = "neutral"
    bins_df["calling_weight"] = np.clip(
        bins_df["bin_weight"].to_numpy(dtype=np.float64) * bins_df["region_risk_weight"].to_numpy(dtype=np.float64),
        0.10,
        None,
    )
    bins_df, calling_background = annotate_calling_z_scores(bins_df, signal_column)

    chrom_frames = []
    all_events = []
    for _, chrom_df in bins_df.groupby("chrom", sort=False):
        chrom_df = chrom_df.copy()
        chrom_df, chrom_events = build_candidate_events(args, chrom_df, logger)
        chrom_frames.append(chrom_df)
        all_events.extend(chrom_events)

    output_bins = pd.concat(chrom_frames, ignore_index=True)
    output_events = pd.DataFrame(all_events)
    if output_events.empty:
        output_events = pd.DataFrame(columns=[field for field in CandidateEvent.__dataclass_fields__])

    summary = {
        "sample_id": args.sample_id,
        "branch": args.branch,
        "correction_model": args.correction_model,
        "binsize": binsize,
        "quality": None if not np.isfinite(quality) else float(quality),
        "bin_count": int(len(output_bins)),
        "candidate_event_count": int(len(output_events)),
        "median_effective_bin_count_per_chrom": float(
            np.median(
                [
                    effective_sample_size(frame["calling_weight"].to_numpy(dtype=np.float64))
                    for _, frame in output_bins.groupby("chrom", sort=False)
                ]
            )
        ) if not output_bins.empty else 0.0,
        "signal_column": signal_column,
        "calling_background_center": float(calling_background["background_center"]),
        "calling_background_scale": float(calling_background["background_scale"]),
        "chrom_shift_z": calling_background["chrom_shift_z"],
    }
    write_table(args.output_bins, output_bins)
    write_table(args.output_candidates, output_events)
    write_json(args.output_summary, summary)
    logger.info("wrote bins=%s candidates=%s summary=%s", args.output_bins, args.output_candidates, args.output_summary)


if __name__ == "__main__":
    main()
