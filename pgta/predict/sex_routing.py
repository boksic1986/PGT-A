#!/biosoftware/miniconda/envs/snakemake_env/bin/python
import argparse
import csv
import math
import re
import shlex
import subprocess
from pathlib import Path

from pgta.core.logging import setup_logger


def run_command(command, logger):
    logger.info("$ %s", " ".join(shlex.quote(part) for part in command))
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    if process.stdout is None:
        raise RuntimeError("Failed to capture command stdout.")

    output_lines = []
    for line in process.stdout:
        line = line.rstrip("\n")
        output_lines.append(line)
        if line:
            logger.info("[cmd] %s", line)

    return_code = process.wait()
    if return_code != 0:
        raise subprocess.CalledProcessError(return_code, command)
    return "\n".join(output_lines)


def parse_gender_output(raw_output):
    matches = re.findall(r"\b(Female|Male|F|M)\b", raw_output, flags=re.IGNORECASE)
    if not matches:
        raise ValueError(f"Unable to parse WisecondorX gender output: {raw_output!r}")

    token = matches[-1].upper()
    wise_gender = "F" if token in {"F", "FEMALE"} else "M"
    sex_call = "XX" if wise_gender == "F" else "XY"
    return sex_call, wise_gender


def sex_call_to_gender(sex_call):
    return "F" if str(sex_call).upper() == "XX" else "M"


def format_optional_float(value):
    if value is None or not math.isfinite(value):
        return ""
    return f"{value:.6f}"


def normalize_chrom_label(label):
    token = str(label).strip()
    if token.lower().startswith("chr"):
        token = token[3:]
    return token.upper()


def parse_idxstats_output(raw_output):
    autosome_depths = []
    x_depth = math.nan
    y_depth = math.nan
    for line in raw_output.splitlines():
        if not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 4:
            continue
        chrom, length_raw, mapped_raw, _ = fields
        if chrom == "*":
            continue
        try:
            length = float(length_raw)
            mapped = float(mapped_raw)
        except ValueError:
            continue
        if length <= 0.0:
            continue
        depth = mapped / length
        chrom_key = normalize_chrom_label(chrom)
        if chrom_key.isdigit() and 1 <= int(chrom_key) <= 22:
            autosome_depths.append(depth)
        elif chrom_key == "X":
            x_depth = depth
        elif chrom_key == "Y":
            y_depth = depth
    autosome_median_depth = float(math.nan)
    if autosome_depths:
        autosome_depths = sorted(float(value) for value in autosome_depths if math.isfinite(value))
        if autosome_depths:
            midpoint = len(autosome_depths) // 2
            if len(autosome_depths) % 2 == 1:
                autosome_median_depth = autosome_depths[midpoint]
            else:
                autosome_median_depth = (autosome_depths[midpoint - 1] + autosome_depths[midpoint]) / 2.0
    x_relative_depth = math.nan
    y_relative_depth = math.nan
    y_to_x_ratio = math.nan
    if math.isfinite(autosome_median_depth) and autosome_median_depth > 0.0:
        if math.isfinite(x_depth):
            x_relative_depth = x_depth / autosome_median_depth
        if math.isfinite(y_depth):
            y_relative_depth = y_depth / autosome_median_depth
    if math.isfinite(x_depth) and x_depth > 0.0 and math.isfinite(y_depth):
        y_to_x_ratio = y_depth / x_depth
    return {
        "autosome_median_depth": autosome_median_depth,
        "chrX_depth": x_depth,
        "chrY_depth": y_depth,
        "x_relative_depth": x_relative_depth,
        "y_relative_depth": y_relative_depth,
        "y_to_x_ratio": y_to_x_ratio,
    }


def infer_bam_sex(metrics, xx_min_x_relative, xy_max_x_relative, xy_min_y_relative, xx_max_y_relative):
    x_relative = metrics.get("x_relative_depth", math.nan)
    y_relative = metrics.get("y_relative_depth", math.nan)
    if not math.isfinite(x_relative) or not math.isfinite(y_relative):
        return "", "bam_depth_unavailable"
    if x_relative >= xx_min_x_relative:
        return "XX", "bam_x_depth_high"
    if x_relative <= xy_max_x_relative and y_relative >= xy_min_y_relative:
        return "XY", "bam_x_low_y_supported"
    if y_relative <= xx_max_y_relative:
        return "XX", "bam_y_depth_low"
    return "", "bam_depth_ambiguous"


def run_samtools_idxstats(samtools, bam_path, logger):
    return run_command([samtools, "idxstats", str(bam_path)], logger)


def resolve_final_sex_call(wise_sex_call, bam_inferred_sex, method):
    normalized_method = str(method).strip().lower()
    if normalized_method == "wisecondorx_only":
        return wise_sex_call, sex_call_to_gender(wise_sex_call), "wisecondorx_only"
    if normalized_method == "bam_depth_only":
        if bam_inferred_sex in {"XX", "XY"}:
            return bam_inferred_sex, sex_call_to_gender(bam_inferred_sex), "bam_depth_only"
        return wise_sex_call, sex_call_to_gender(wise_sex_call), "bam_depth_only_fallback_wisecondorx"
    if bam_inferred_sex in {"XX", "XY"} and bam_inferred_sex != wise_sex_call:
        return bam_inferred_sex, sex_call_to_gender(bam_inferred_sex), "bam_depth_override"
    if bam_inferred_sex in {"XX", "XY"}:
        return wise_sex_call, sex_call_to_gender(wise_sex_call), "wisecondorx_bam_consensus"
    return wise_sex_call, sex_call_to_gender(wise_sex_call), "wisecondorx_fallback"


def write_gender_tsv(
    output_tsv,
    sample_id,
    sex_call,
    wise_gender,
    predict_gender,
    sex_call_source,
    bam_inferred_sex,
    bam_metrics,
    bam_reason,
    raw_output,
):
    output_path = Path(output_tsv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    compact_output = re.sub(r"\s+", " ", raw_output).strip()
    with open(output_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "sample_id",
                "sex_call",
                "wise_gender",
                "predict_gender",
                "sex_call_source",
                "bam_inferred_sex",
                "bam_inference_reason",
                "bam_autosome_median_depth",
                "bam_chrX_depth",
                "bam_chrY_depth",
                "bam_x_relative_depth",
                "bam_y_relative_depth",
                "bam_y_to_x_ratio",
                "raw_output",
            ]
        )
        writer.writerow(
            [
                sample_id,
                sex_call,
                wise_gender,
                predict_gender,
                sex_call_source,
                bam_inferred_sex,
                bam_reason,
                format_optional_float(bam_metrics.get("autosome_median_depth", math.nan)),
                format_optional_float(bam_metrics.get("chrX_depth", math.nan)),
                format_optional_float(bam_metrics.get("chrY_depth", math.nan)),
                format_optional_float(bam_metrics.get("x_relative_depth", math.nan)),
                format_optional_float(bam_metrics.get("y_relative_depth", math.nan)),
                format_optional_float(bam_metrics.get("y_to_x_ratio", math.nan)),
                compact_output,
            ]
        )


def run_wisecondorx_gender(
    wisecondorx,
    sample_npz,
    gender_reference,
    output_tsv,
    sample_id,
    logger,
    bam_path="",
    samtools="",
    method="wisecondorx_plus_bam_depth",
    xx_min_x_relative=0.95,
    xy_max_x_relative=0.80,
    xy_min_y_relative=0.20,
    xx_max_y_relative=0.15,
):
    raw_output = run_command(
        [
            wisecondorx,
            "gender",
            str(sample_npz),
            str(gender_reference),
        ],
        logger,
    )
    wise_sex_call, wise_gender = parse_gender_output(raw_output)
    bam_metrics = {
        "autosome_median_depth": math.nan,
        "chrX_depth": math.nan,
        "chrY_depth": math.nan,
        "x_relative_depth": math.nan,
        "y_relative_depth": math.nan,
        "y_to_x_ratio": math.nan,
    }
    bam_inferred_sex = ""
    bam_reason = "bam_depth_skipped"
    if bam_path and samtools:
        idxstats_output = run_samtools_idxstats(samtools, bam_path, logger)
        bam_metrics = parse_idxstats_output(idxstats_output)
        bam_inferred_sex, bam_reason = infer_bam_sex(
            bam_metrics,
            xx_min_x_relative=xx_min_x_relative,
            xy_max_x_relative=xy_max_x_relative,
            xy_min_y_relative=xy_min_y_relative,
            xx_max_y_relative=xx_max_y_relative,
        )
    sex_call, predict_gender, sex_call_source = resolve_final_sex_call(
        wise_sex_call=wise_sex_call,
        bam_inferred_sex=bam_inferred_sex,
        method=method,
    )
    write_gender_tsv(
        output_tsv=output_tsv,
        sample_id=sample_id,
        sex_call=sex_call,
        wise_gender=wise_gender,
        predict_gender=predict_gender,
        sex_call_source=sex_call_source,
        bam_inferred_sex=bam_inferred_sex,
        bam_metrics=bam_metrics,
        bam_reason=bam_reason,
        raw_output=raw_output,
    )
    logger.info(
        "gender call completed: sample=%s sex_call=%s wise_gender=%s bam_inferred_sex=%s source=%s x_rel=%s y_rel=%s",
        sample_id,
        sex_call,
        wise_gender,
        bam_inferred_sex or "NA",
        sex_call_source,
        format_optional_float(bam_metrics["x_relative_depth"]) or "NA",
        format_optional_float(bam_metrics["y_relative_depth"]) or "NA",
    )


def main():
    parser = argparse.ArgumentParser(description="Run WisecondorX gender and write a structured TSV result.")
    parser.add_argument("--wisecondorx", required=True)
    parser.add_argument("--sample-npz", required=True)
    parser.add_argument("--gender-reference", required=True)
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--bam", default="")
    parser.add_argument("--samtools", default="")
    parser.add_argument("--method", default="wisecondorx_plus_bam_depth")
    parser.add_argument("--xx-min-x-relative", type=float, default=0.95)
    parser.add_argument("--xy-max-x-relative", type=float, default=0.80)
    parser.add_argument("--xy-min-y-relative", type=float, default=0.20)
    parser.add_argument("--xx-max-y-relative", type=float, default=0.15)
    parser.add_argument("--log", default="", help="Optional log file path")
    args = parser.parse_args()

    logger = setup_logger("wisecondorx_gender", args.log or None)
    run_wisecondorx_gender(
        wisecondorx=args.wisecondorx,
        sample_npz=args.sample_npz,
        gender_reference=args.gender_reference,
        output_tsv=args.output_tsv,
        sample_id=args.sample_id,
        logger=logger,
        bam_path=args.bam,
        samtools=args.samtools,
        method=args.method,
        xx_min_x_relative=args.xx_min_x_relative,
        xy_max_x_relative=args.xy_max_x_relative,
        xy_min_y_relative=args.xy_min_y_relative,
        xx_max_y_relative=args.xx_max_y_relative,
    )


if __name__ == "__main__":
    main()
