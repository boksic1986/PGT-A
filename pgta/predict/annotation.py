from __future__ import annotations

import argparse
import json
import logging
import sqlite3
from pathlib import Path

import pandas as pd

from pgta.core.logging import setup_logger


ANNOTATION_COLUMNS = [
    "event_id",
    "candidate_id",
    "sample_id",
    "chrom",
    "start",
    "end",
    "state",
    "event_layer",
    "cytoband",
    "genes",
    "gene_number",
    "gene_location",
    "omim_genes",
    "omim_phenotypes",
    "hpo_terms",
    "region_context",
    "annotation_backend",
    "annotation_status",
    "annotation_bundle_id",
    "genome_build",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Annotate PGTA CNV report/review events with a PGTA-owned SQLite bundle."
    )
    parser.add_argument("--event-tsv", action="append", default=[])
    parser.add_argument("--classification-tsv", action="append", default=[])
    parser.add_argument("--branch-s-evidence", action="append", default=[])
    parser.add_argument("--bundle-db", required=True)
    parser.add_argument("--genome-build", default="hg19")
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--output-json", required=True)
    parser.add_argument("--log", default="")
    return parser.parse_args()


def ensure_parent(path_value):
    path = Path(path_value)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def read_table(path_value):
    if not path_value:
        return pd.DataFrame()
    path = Path(path_value)
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def normalize_chromosome(value) -> str:
    text = str(value).strip()
    if not text:
        return ""
    if text.lower().startswith("chr"):
        suffix = text[3:]
    else:
        suffix = text
    suffix_upper = suffix.upper()
    if suffix_upper == "23":
        suffix_upper = "X"
    elif suffix_upper == "24":
        suffix_upper = "Y"
    elif suffix_upper in {"M", "MT", "25"}:
        suffix_upper = "M"
    elif suffix_upper.isdigit():
        suffix_upper = str(int(suffix_upper))
    return f"chr{suffix_upper}"


def text_or_empty(value) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except TypeError:
        pass
    return str(value).strip()


def numeric_or_none(value):
    number = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if pd.isna(number):
        return None
    return int(number)


def normalize_event_frame(frame: pd.DataFrame, event_layer: str) -> pd.DataFrame:
    if frame is None or frame.empty:
        return pd.DataFrame(columns=ANNOTATION_COLUMNS)
    rows = []
    for idx, row in frame.iterrows():
        chrom = normalize_chromosome(row.get("chrom", row.get("chr", "")))
        start = numeric_or_none(row.get("start", row.get("pos_start", "")))
        end = numeric_or_none(row.get("end", row.get("pos_end", "")))
        if not chrom or start is None or end is None or end <= start:
            continue
        sample_id = text_or_empty(row.get("sample_id", row.get("sample", "")))
        state = text_or_empty(row.get("state", row.get("type", "")))
        candidate_id = text_or_empty(row.get("candidate_id", ""))
        event_id = text_or_empty(row.get("event_id", "")) or candidate_id
        if not event_id:
            event_id = f"{sample_id}:{chrom}:{start}-{end}:{state}:{event_layer}:{idx}"
        row_event_layer = text_or_empty(row.get("event_layer", "")) or str(event_layer)
        rows.append(
            {
                "event_id": event_id,
                "candidate_id": candidate_id,
                "sample_id": sample_id,
                "chrom": chrom,
                "start": start,
                "end": end,
                "state": state,
                "event_layer": row_event_layer,
            }
        )
    return pd.DataFrame(rows)


def collect_events(event_paths=None, classification_paths=None, branch_s_evidence_paths=None) -> pd.DataFrame:
    frames = []
    for path_value in event_paths or []:
        frame = read_table(path_value)
        if not frame.empty:
            frames.append(normalize_event_frame(frame, "autosomal_report"))
    for path_value in classification_paths or []:
        frame = read_table(path_value)
        if frame.empty:
            continue
        layer = frame.get("v2_report_layer_class", pd.Series("", index=frame.index)).fillna("").astype(str)
        visibility = frame.get("v2_report_visibility", pd.Series("", index=frame.index)).fillna("").astype(str)
        include = (
            layer.isin({"report_event", "internal_review_event", "technical_risk_review", "background_unknown_review", "branch_s_event"})
            | visibility.isin({"report_strong_event", "report_weak_event", "internal_review_event", "branch_s_event"})
        ) & ~layer.eq("filtered_event") & ~visibility.eq("filtered_event")
        if include.any():
            prepared = frame.loc[include].copy()
            prepared["event_layer"] = layer.loc[include].replace("", "internal_review_event")
            frames.append(
                pd.concat(
                    [
                        normalize_event_frame(subframe, str(layer_name))
                        for layer_name, subframe in prepared.groupby("event_layer", dropna=False)
                    ],
                    ignore_index=True,
                )
            )
    for path_value in branch_s_evidence_paths or []:
        frame = read_table(path_value)
        if frame.empty:
            continue
        counts = pd.to_numeric(frame.get("a_candidate_count", pd.Series(0, index=frame.index)), errors="coerce").fillna(0)
        prepared = frame.loc[counts.gt(0)].copy()
        if prepared.empty:
            continue
        prepared["event_id"] = prepared.apply(
            lambda row: (
                f"{text_or_empty(row.get('sample_id'))}:branch_s:"
                f"{text_or_empty(row.get('region_class')) or normalize_chromosome(row.get('chrom'))}:"
                f"{int(row.get('start'))}-{int(row.get('end'))}"
            ),
            axis=1,
        )
        prepared["state"] = prepared.get("region_class", pd.Series("", index=prepared.index)).astype(str)
        frames.append(normalize_event_frame(prepared, "branch_s_review"))
    if not frames:
        return pd.DataFrame(columns=ANNOTATION_COLUMNS[:8])
    combined = pd.concat(frames, ignore_index=True)
    combined["_priority"] = combined["event_layer"].map(
        {
            "autosomal_report": 0,
            "report_event": 0,
            "internal_review_event": 1,
            "technical_risk_review": 1,
            "background_unknown_review": 1,
            "branch_s_event": 2,
            "branch_s_review": 2,
        }
    ).fillna(3)
    combined = combined.sort_values("_priority")
    key_cols = ["sample_id", "chrom", "start", "end", "state"]
    combined = combined.drop_duplicates(subset=key_cols, keep="first").drop(columns=["_priority"])
    return combined.reset_index(drop=True)


def unique_join(values) -> str:
    result = []
    for value in values:
        text = text_or_empty(value)
        if text and text not in result:
            result.append(text)
    return ",".join(result)


def table_exists(conn: sqlite3.Connection, table_name: str) -> bool:
    row = conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name=?",
        (table_name,),
    ).fetchone()
    return row is not None


def query_interval_table(conn, table_name, chrom, start, end):
    if not table_exists(conn, table_name):
        return []
    return conn.execute(
        f"SELECT * FROM {table_name} WHERE chrom = ? AND start < ? AND end > ?",
        (chrom, int(end), int(start)),
    ).fetchall()


def fetch_metadata(conn: sqlite3.Connection) -> dict[str, str]:
    if not table_exists(conn, "bundle_metadata"):
        return {}
    rows = conn.execute("SELECT key, value FROM bundle_metadata").fetchall()
    return {str(key): str(value) for key, value in rows}


def fetch_interval_records(conn: sqlite3.Connection, table_name: str, chrom: str, start: int, end: int) -> list[dict]:
    if not table_exists(conn, table_name):
        return []
    cursor = conn.execute(
        f"SELECT * FROM {table_name} WHERE chrom = ? AND start < ? AND end > ?",
        (chrom, int(end), int(start)),
    )
    columns = [item[0] for item in cursor.description or []]
    return [dict(zip(columns, row)) for row in cursor.fetchall()]


def fetch_gene_link_records(conn: sqlite3.Connection, table_name: str, genes: list[str]) -> list[dict]:
    if not genes or not table_exists(conn, table_name):
        return []
    placeholders = ",".join(["?"] * len(genes))
    cursor = conn.execute(
        f"SELECT * FROM {table_name} WHERE gene_symbol IN ({placeholders})",
        genes,
    )
    columns = [item[0] for item in cursor.description or []]
    return [dict(zip(columns, row)) for row in cursor.fetchall()]


def missing_annotation_rows(events: pd.DataFrame, genome_build: str, status: str) -> pd.DataFrame:
    output = events.copy()
    for column in ANNOTATION_COLUMNS:
        if column not in output.columns:
            output[column] = ""
    output["annotation_backend"] = "pgta_sqlite"
    output["annotation_status"] = status
    output["annotation_bundle_id"] = ""
    output["genome_build"] = str(genome_build)
    output["gene_number"] = 0
    return output[ANNOTATION_COLUMNS]


def annotate_events_dataframe(events: pd.DataFrame, bundle_db, genome_build="hg19") -> pd.DataFrame:
    if events is None or events.empty:
        return pd.DataFrame(columns=ANNOTATION_COLUMNS)
    normalized = normalize_event_frame(events, "autosomal_report")
    bundle_path = Path(bundle_db)
    if not bundle_path.exists() or bundle_path.stat().st_size == 0:
        return missing_annotation_rows(normalized, genome_build, "missing_backend")

    rows = []
    conn = sqlite3.connect(bundle_path)
    try:
        metadata = fetch_metadata(conn)
        bundle_id = metadata.get("bundle_id", bundle_path.name)
        for _, event in normalized.iterrows():
            chrom = normalize_chromosome(event["chrom"])
            start = int(event["start"])
            end = int(event["end"])
            cytobands = fetch_interval_records(conn, "cytoband", chrom, start, end)
            genes = fetch_interval_records(conn, "gene", chrom, start, end)
            region_context = fetch_interval_records(conn, "region_context", chrom, start, end)
            gene_symbols = [text_or_empty(item.get("gene_symbol")) for item in genes if text_or_empty(item.get("gene_symbol"))]
            omim = fetch_gene_link_records(conn, "omim_gene", gene_symbols)
            hpo = fetch_gene_link_records(conn, "hpo_term", gene_symbols)
            payload = event.to_dict()
            payload.update(
                {
                    "chrom": chrom,
                    "cytoband": unique_join(item.get("name", "") for item in cytobands),
                    "genes": unique_join(gene_symbols),
                    "gene_number": int(len(set(gene_symbols))),
                    "gene_location": unique_join(item.get("gene_location", "") for item in genes),
                    "omim_genes": unique_join(
                        f"{item.get('gene_symbol')}:{item.get('omim_id')}" for item in omim
                    ),
                    "omim_phenotypes": unique_join(item.get("phenotype", "") for item in omim),
                    "hpo_terms": unique_join(
                        f"{item.get('hpo_id')}:{item.get('hpo_name')}" for item in hpo
                    ),
                    "region_context": unique_join(item.get("context_label", "") for item in region_context),
                    "annotation_backend": "pgta_sqlite",
                    "annotation_status": "annotated",
                    "annotation_bundle_id": bundle_id,
                    "genome_build": str(genome_build),
                }
            )
            rows.append(payload)
    finally:
        conn.close()
    output = pd.DataFrame(rows)
    for column in ANNOTATION_COLUMNS:
        if column not in output.columns:
            output[column] = ""
    return output[ANNOTATION_COLUMNS]


def summarize_annotations(annotated: pd.DataFrame, bundle_db: str, genome_build: str) -> dict:
    if annotated is None or annotated.empty:
        return {
            "status": "empty",
            "row_count": 0,
            "bundle_db": str(bundle_db),
            "genome_build": str(genome_build),
        }
    status_counts = annotated["annotation_status"].fillna("").astype(str).value_counts().to_dict()
    return {
        "status": "completed",
        "row_count": int(len(annotated)),
        "sample_count": int(annotated["sample_id"].astype(str).nunique()) if "sample_id" in annotated.columns else 0,
        "status_counts": {str(key): int(value) for key, value in status_counts.items()},
        "bundle_db": str(bundle_db),
        "genome_build": str(genome_build),
        "annotation_backend": "pgta_sqlite",
        "annotation_changes_event_classification": False,
        "annotation_changes_filtering": False,
        "cnvpro_filter_logic_imported": False,
    }


def main():
    args = parse_args()
    logger = setup_logger("cnv_event_annotation", args.log or None)
    events = collect_events(
        event_paths=args.event_tsv,
        classification_paths=args.classification_tsv,
        branch_s_evidence_paths=args.branch_s_evidence,
    )
    annotated = annotate_events_dataframe(events, bundle_db=args.bundle_db, genome_build=args.genome_build)
    ensure_parent(args.output_tsv)
    annotated.to_csv(args.output_tsv, sep="\t", index=False)
    summary = summarize_annotations(annotated, bundle_db=args.bundle_db, genome_build=args.genome_build)
    ensure_parent(args.output_json).write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    logger.info("wrote CNV event annotations rows=%d to %s", len(annotated), args.output_tsv)


if __name__ == "__main__":
    main()
