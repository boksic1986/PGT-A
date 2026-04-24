#!/biosoftware/miniconda/envs/snakemake_env/bin/python
import argparse
import csv
import re
from pathlib import Path

from pgta.core.logging import setup_logger


def natural_key(value):
    parts = re.split(r"(\d+)", str(value))
    key = []
    for part in parts:
        if not part:
            continue
        if part.isdigit():
            key.append((0, int(part)))
        else:
            key.append((1, part.lower()))
    return tuple(key)


def parse_decisions(raw):
    values = [item.strip().upper() for item in str(raw).split(",") if item.strip()]
    if not values:
        raise ValueError("Empty decision filter.")
    return values


def load_selected_samples(summary_tsv, decisions):
    decisions = {item.upper() for item in decisions}
    selected = []
    with open(summary_tsv, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = {"sample_id", "qc_decision"} - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"Missing required columns in {summary_tsv}: {','.join(sorted(missing))}"
            )
        for row in reader:
            sample_id = str(row.get("sample_id", "")).strip()
            decision = str(row.get("qc_decision", "")).strip().upper()
            if sample_id and decision in decisions:
                selected.append(sample_id)
    selected = sorted(dict.fromkeys(selected), key=natural_key)
    if not selected:
        raise ValueError(
            f"No samples matched decisions={','.join(sorted(decisions))} in {summary_tsv}"
        )
    return selected


def main():
    parser = argparse.ArgumentParser(description="Select reference cohort sample IDs from baseline QC summary.")
    parser.add_argument("--summary-tsv", required=True)
    parser.add_argument("--decisions", required=True, help="Comma-separated QC decisions, e.g. PASS,WARN")
    parser.add_argument("--output", required=True)
    parser.add_argument("--log", default="")
    args = parser.parse_args()

    logger = setup_logger("select_reference_cohorts", args.log or None)
    decisions = parse_decisions(args.decisions)
    selected = load_selected_samples(args.summary_tsv, decisions)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("".join(f"{sample_id}\n" for sample_id in selected), encoding="utf-8")
    logger.info(
        "reference cohort selected: decisions=%s samples=%d output=%s",
        ",".join(decisions),
        len(selected),
        output_path,
    )


if __name__ == "__main__":
    main()
