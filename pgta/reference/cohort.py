#!/biosoftware/miniconda/envs/snakemake_env/bin/python
import argparse
import csv
import re
from pathlib import Path

from pgta.core.logging import setup_logger

RESEQUENCING_REQUIRED_COLUMNS = {
    "sample_id",
    "source_sample_id",
    "library_id",
    "fastq_r1",
    "fastq_r2",
    "sex_group",
    "batch_group",
    "reseq_reason",
    "replacement_policy",
    "status",
    "decision_reason",
}
RESEQUENCING_VALID_POLICIES = {"replace_original", "compare_then_replace"}
RESEQUENCING_DEFAULT_REFERENCE_STATUSES = ("promoted",)


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


def load_resequencing_manifest(path_value):
    path = Path(path_value)
    rows = []
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = RESEQUENCING_REQUIRED_COLUMNS - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"Missing required columns in {path}: {','.join(sorted(missing))}"
            )
        for row_number, row in enumerate(reader, start=2):
            normalized = {key: str(row.get(key, "")).strip() for key in RESEQUENCING_REQUIRED_COLUMNS}
            sample_id = normalized["sample_id"]
            source_sample_id = normalized["source_sample_id"]
            status = normalized["status"].lower()
            policy = normalized["replacement_policy"].lower()
            sex_group = normalized["sex_group"].upper()
            if not sample_id:
                raise ValueError(f"Empty sample_id in {path}:{row_number}")
            if not source_sample_id:
                raise ValueError(f"Empty source_sample_id in {path}:{row_number}")
            if sex_group not in {"XX", "XY"}:
                raise ValueError(f"Invalid sex_group for {sample_id} in {path}:{row_number}: {sex_group}")
            if policy not in RESEQUENCING_VALID_POLICIES:
                raise ValueError(
                    f"Invalid replacement_policy for {sample_id} in {path}:{row_number}: {policy}"
                )
            normalized["status"] = status
            normalized["replacement_policy"] = policy
            normalized["sex_group"] = sex_group
            rows.append(normalized)
    return rows


def apply_resequencing_manifest(selected, manifest_rows, allowed_statuses=None):
    allowed = {
        str(status).strip().lower()
        for status in (allowed_statuses or RESEQUENCING_DEFAULT_REFERENCE_STATUSES)
        if str(status).strip()
    }
    if not allowed:
        raise ValueError("resequencing allowed statuses must not be empty.")

    selected_lookup = set(selected)
    rows_by_sample = {}
    replacement_sources = set()
    for row in manifest_rows:
        sample_id = row["sample_id"]
        if sample_id in rows_by_sample:
            raise ValueError(f"Duplicate re-sequencing sample_id in manifest: {sample_id}")
        rows_by_sample[sample_id] = row
        if row["status"] in allowed and row["replacement_policy"] == "replace_original":
            replacement_sources.add(row["source_sample_id"])

    filtered = []
    for sample_id in selected:
        row = rows_by_sample.get(sample_id)
        if row is not None and row["status"] not in allowed:
            continue
        if sample_id in replacement_sources and sample_id not in rows_by_sample:
            continue
        if sample_id not in selected_lookup:
            continue
        filtered.append(sample_id)
    return sorted(dict.fromkeys(filtered), key=natural_key)


def load_selected_samples(summary_tsv, decisions, resequencing_manifest=None, resequencing_allowed_statuses=None):
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
    if resequencing_manifest:
        selected = apply_resequencing_manifest(
            selected,
            load_resequencing_manifest(resequencing_manifest),
            allowed_statuses=resequencing_allowed_statuses,
        )
    if not selected:
        raise ValueError(
            f"No samples matched decisions={','.join(sorted(decisions))} in {summary_tsv}"
        )
    return selected


def main():
    parser = argparse.ArgumentParser(description="Select reference cohort sample IDs from baseline QC summary.")
    parser.add_argument("--summary-tsv", required=True)
    parser.add_argument("--decisions", required=True, help="Comma-separated QC decisions, e.g. PASS,WARN")
    parser.add_argument("--resequencing-manifest", default="")
    parser.add_argument("--resequencing-allowed-status", action="append", default=[])
    parser.add_argument("--output", required=True)
    parser.add_argument("--log", default="")
    args = parser.parse_args()

    logger = setup_logger("select_reference_cohorts", args.log or None)
    decisions = parse_decisions(args.decisions)
    selected = load_selected_samples(
        args.summary_tsv,
        decisions,
        resequencing_manifest=args.resequencing_manifest or None,
        resequencing_allowed_statuses=args.resequencing_allowed_status or None,
    )

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
