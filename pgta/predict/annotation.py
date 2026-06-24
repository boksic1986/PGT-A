from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import logging
import sqlite3
import tempfile
from datetime import datetime, timezone
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
    "gene_source_status",
    "omim_source_status",
    "hpo_source_status",
    "annotation_backend",
    "annotation_status",
    "annotation_bundle_id",
    "genome_build",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Annotate PGTA CNV report/review events with a PGTA-owned SQLite bundle."
    )
    parser.add_argument("--build-bundle", action="store_true")
    parser.add_argument("--event-tsv", action="append", default=[])
    parser.add_argument("--classification-tsv", action="append", default=[])
    parser.add_argument("--branch-s-evidence", action="append", default=[])
    parser.add_argument("--bundle-db", required=True)
    parser.add_argument("--genome-build", default="hg19")
    parser.add_argument("--output-tsv", default="")
    parser.add_argument("--output-json", default="")
    parser.add_argument("--log", default="")
    parser.add_argument("--bundle-id", default="")
    parser.add_argument("--cytoband-source", default="")
    parser.add_argument("--refgene-bed", default="")
    parser.add_argument("--omim-source", action="append", default=[])
    parser.add_argument("--hpo-source", action="append", default=[])
    parser.add_argument("--source-manifest", default="")
    parser.add_argument("--source-root", default="")
    parser.add_argument("--pgta-source-root", default="")
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


def open_text(path_value):
    path = Path(path_value)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def md5sum(path_value) -> str:
    digest = hashlib.md5()
    with Path(path_value).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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


def load_cytoband_rows(path_value) -> list[tuple]:
    rows = []
    if not path_value:
        return rows
    with open_text(path_value) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            chrom = normalize_chromosome(parts[0])
            start = numeric_or_none(parts[1])
            end = numeric_or_none(parts[2])
            if not chrom or start is None or end is None or end <= start:
                continue
            rows.append((chrom, start, end, text_or_empty(parts[3]), text_or_empty(parts[4])))
    return rows


def load_refgene_rows(path_value) -> list[tuple]:
    rows = []
    seen = set()
    if not path_value:
        return rows
    with open_text(path_value) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            chrom = normalize_chromosome(parts[0])
            start = numeric_or_none(parts[1])
            end = numeric_or_none(parts[2])
            gene_symbol = text_or_empty(parts[4])
            transcript_id = text_or_empty(parts[5])
            strand = text_or_empty(parts[3])
            if not chrom or start is None or end is None or end <= start or not gene_symbol:
                continue
            key = (chrom, start, end, gene_symbol, transcript_id, strand)
            if key in seen:
                continue
            seen.add(key)
            rows.append((chrom, start, end, gene_symbol, "refGene", transcript_id, strand))
    return rows


def load_omim_rows(paths) -> list[tuple]:
    rows = []
    seen = set()
    for path_value in paths or []:
        path = Path(path_value)
        if not path.exists():
            continue
        with open_text(path) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames:
                continue
            field_lookup = {name.lower(): name for name in reader.fieldnames}
            gene_col = field_lookup.get("genes") or field_lookup.get("gene") or field_lookup.get("gene_symbol")
            mim_col = field_lookup.get("mim number") or field_lookup.get("mim_number") or field_lookup.get("omim_id")
            phenotype_col = field_lookup.get("phenotypes") or field_lookup.get("phenotype")
            inheritance_col = field_lookup.get("inheritance")
            if not gene_col:
                continue
            for row in reader:
                gene_symbol = text_or_empty(row.get(gene_col, ""))
                if not gene_symbol:
                    continue
                omim_id = text_or_empty(row.get(mim_col, "")) if mim_col else ""
                phenotype = text_or_empty(row.get(phenotype_col, "")) if phenotype_col else ""
                inheritance = text_or_empty(row.get(inheritance_col, "")) if inheritance_col else ""
                if not omim_id and not phenotype:
                    continue
                key = (gene_symbol, omim_id, phenotype, inheritance)
                if key in seen:
                    continue
                seen.add(key)
                rows.append((gene_symbol, omim_id, phenotype, inheritance, path.name))
    return rows


def load_hpo_rows(paths) -> list[tuple]:
    rows = []
    seen = set()
    for path_value in paths or []:
        path = Path(path_value)
        if not path.exists():
            continue
        with open_text(path) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames:
                continue
            field_lookup = {name.lower(): name for name in reader.fieldnames}
            gene_col = field_lookup.get("gene_symbol") or field_lookup.get("genes") or field_lookup.get("gene")
            hpo_col = field_lookup.get("hpo_id") or field_lookup.get("hpo")
            name_col = field_lookup.get("hpo_name") or field_lookup.get("hpo term") or field_lookup.get("hpo_term")
            if not gene_col or not hpo_col:
                continue
            for row in reader:
                gene_symbol = text_or_empty(row.get(gene_col, ""))
                hpo_id = text_or_empty(row.get(hpo_col, ""))
                hpo_name = text_or_empty(row.get(name_col, "")) if name_col else ""
                if not gene_symbol or not hpo_id:
                    continue
                key = (gene_symbol, hpo_id, hpo_name)
                if key in seen:
                    continue
                seen.add(key)
                rows.append((gene_symbol, hpo_id, hpo_name, path.name))
    return rows


def source_manifest_row(path_value, resource_type, genome_build, source_root="", pgta_source_root="") -> dict:
    path = Path(path_value)
    source_path = ""
    pgta_path = str(path)
    path_md5 = md5sum(path) if path.exists() else ""
    if source_root and pgta_source_root:
        try:
            relative = path.resolve().relative_to(Path(pgta_source_root).resolve())
            candidate = Path(source_root) / relative
            if candidate.exists():
                source_path = str(candidate)
        except ValueError:
            source_path = ""
        if not source_path and path.exists():
            matches = []
            for candidate in Path(source_root).rglob(path.name):
                try:
                    if candidate.is_file() and md5sum(candidate) == path_md5:
                        matches.append(str(candidate))
                except OSError:
                    continue
            if len(matches) == 1:
                source_path = matches[0]
    return {
        "resource_type": str(resource_type),
        "genome_build": str(genome_build),
        "source_path": source_path,
        "pgta_path": pgta_path,
        "md5": path_md5,
        "size_bytes": int(path.stat().st_size) if path.exists() else 0,
        "copy_date": datetime.now(timezone.utc).date().isoformat(),
    }


def write_source_manifest(rows: list[dict], path_value) -> None:
    if not path_value:
        return
    path = ensure_parent(path_value)
    columns = [
        "resource_type",
        "genome_build",
        "source_path",
        "pgta_path",
        "md5",
        "size_bytes",
        "copy_date",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def build_annotation_bundle(
    bundle_db,
    cytoband_source="",
    refgene_bed="",
    omim_sources=None,
    hpo_sources=None,
    source_manifest="",
    source_root="",
    pgta_source_root="",
    bundle_id="",
    genome_build="hg19",
) -> dict:
    bundle_path = ensure_parent(bundle_db)
    cytoband_rows = load_cytoband_rows(cytoband_source)
    gene_rows = load_refgene_rows(refgene_bed)
    omim_rows = load_omim_rows(omim_sources or [])
    hpo_rows = load_hpo_rows(hpo_sources or [])
    if not cytoband_rows:
        raise ValueError("cytoband_source did not yield any cytoband rows")
    if refgene_bed and not gene_rows:
        raise ValueError("refgene_bed did not yield any gene rows")
    if omim_sources and not omim_rows:
        raise ValueError("omim_sources did not yield any OMIM rows")

    metadata = {
        "bundle_id": bundle_id or bundle_path.stem,
        "genome_build": str(genome_build),
        "cytoband_source_status": "ready" if cytoband_rows else "missing_cytoband_source",
        "gene_source_status": "ready" if gene_rows else "missing_gene_source",
        "omim_source_status": "ready" if omim_rows else "missing_omim_source",
        "hpo_source_status": "ready" if hpo_rows else "missing_hpo_source",
        "cytoband_row_count": str(len(cytoband_rows)),
        "gene_row_count": str(len(gene_rows)),
        "omim_row_count": str(len(omim_rows)),
        "hpo_row_count": str(len(hpo_rows)),
        "bundle_build_time_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }

    with tempfile.NamedTemporaryFile(
        prefix=f"{bundle_path.name}.",
        suffix=".tmp",
        dir=str(bundle_path.parent),
        delete=False,
    ) as handle:
        tmp_path = Path(handle.name)
    conn = sqlite3.connect(tmp_path)
    try:
        conn.execute("CREATE TABLE cytoband (chrom TEXT, start INTEGER, end INTEGER, name TEXT, stain TEXT)")
        conn.execute(
            "CREATE TABLE gene (chrom TEXT, start INTEGER, end INTEGER, gene_symbol TEXT, gene_location TEXT, transcript_id TEXT, strand TEXT)"
        )
        conn.execute(
            "CREATE TABLE omim_gene (gene_symbol TEXT, omim_id TEXT, phenotype TEXT, inheritance TEXT, source TEXT)"
        )
        conn.execute("CREATE TABLE hpo_term (gene_symbol TEXT, hpo_id TEXT, hpo_name TEXT, source TEXT)")
        conn.execute(
            "CREATE TABLE region_context (chrom TEXT, start INTEGER, end INTEGER, context_label TEXT, context_reason TEXT)"
        )
        conn.execute("CREATE TABLE bundle_metadata (key TEXT PRIMARY KEY, value TEXT)")
        conn.executemany("INSERT INTO cytoband VALUES (?, ?, ?, ?, ?)", cytoband_rows)
        conn.executemany("INSERT INTO gene VALUES (?, ?, ?, ?, ?, ?, ?)", gene_rows)
        conn.executemany("INSERT INTO omim_gene VALUES (?, ?, ?, ?, ?)", omim_rows)
        conn.executemany("INSERT INTO hpo_term VALUES (?, ?, ?, ?)", hpo_rows)
        conn.executemany(
            "INSERT INTO bundle_metadata VALUES (?, ?)",
            [(str(key), str(value)) for key, value in sorted(metadata.items())],
        )
        conn.execute("CREATE INDEX idx_cytoband_interval ON cytoband(chrom, start, end)")
        conn.execute("CREATE INDEX idx_gene_interval ON gene(chrom, start, end)")
        conn.execute("CREATE INDEX idx_omim_gene ON omim_gene(gene_symbol)")
        conn.execute("CREATE INDEX idx_hpo_gene ON hpo_term(gene_symbol)")
        conn.commit()
    finally:
        conn.close()
    tmp_path.replace(bundle_path)

    manifest_rows = []
    if cytoband_source:
        manifest_rows.append(source_manifest_row(cytoband_source, "cytoband", genome_build, source_root, pgta_source_root))
    if refgene_bed:
        manifest_rows.append(source_manifest_row(refgene_bed, "refgene", genome_build, source_root, pgta_source_root))
    for path_value in omim_sources or []:
        manifest_rows.append(source_manifest_row(path_value, "omim", genome_build, source_root, pgta_source_root))
    for path_value in hpo_sources or []:
        manifest_rows.append(source_manifest_row(path_value, "hpo", genome_build, source_root, pgta_source_root))
    write_source_manifest(manifest_rows, source_manifest)

    return {
        "status": "completed",
        "bundle_db": str(bundle_path),
        "bundle_id": metadata["bundle_id"],
        "genome_build": str(genome_build),
        "cytoband_row_count": len(cytoband_rows),
        "gene_row_count": len(gene_rows),
        "omim_row_count": len(omim_rows),
        "hpo_row_count": len(hpo_rows),
        "gene_source_status": metadata["gene_source_status"],
        "omim_source_status": metadata["omim_source_status"],
        "hpo_source_status": metadata["hpo_source_status"],
        "source_manifest": str(source_manifest) if source_manifest else "",
    }


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
    output["gene_source_status"] = status
    output["omim_source_status"] = status
    output["hpo_source_status"] = status
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
        gene_source_status = metadata.get("gene_source_status", "unknown")
        omim_source_status = metadata.get("omim_source_status", "unknown")
        hpo_source_status = metadata.get("hpo_source_status", "unknown")
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
                    "gene_source_status": gene_source_status,
                    "omim_source_status": omim_source_status,
                    "hpo_source_status": hpo_source_status,
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
    gene_source_status = unique_join(annotated.get("gene_source_status", pd.Series("", index=annotated.index)).fillna("").astype(str))
    omim_source_status = unique_join(annotated.get("omim_source_status", pd.Series("", index=annotated.index)).fillna("").astype(str))
    hpo_source_status = unique_join(annotated.get("hpo_source_status", pd.Series("", index=annotated.index)).fillna("").astype(str))
    bundle_ids = unique_join(annotated.get("annotation_bundle_id", pd.Series("", index=annotated.index)).fillna("").astype(str))
    return {
        "status": "completed",
        "row_count": int(len(annotated)),
        "sample_count": int(annotated["sample_id"].astype(str).nunique()) if "sample_id" in annotated.columns else 0,
        "status_counts": {str(key): int(value) for key, value in status_counts.items()},
        "bundle_db": str(bundle_db),
        "annotation_bundle_id": bundle_ids,
        "genome_build": str(genome_build),
        "annotation_backend": "pgta_sqlite",
        "gene_source_status": gene_source_status,
        "omim_source_status": omim_source_status,
        "hpo_source_status": hpo_source_status,
        "gene_nonempty_row_count": int(annotated.get("genes", pd.Series("", index=annotated.index)).fillna("").astype(str).ne("").sum()),
        "omim_nonempty_row_count": int(annotated.get("omim_genes", pd.Series("", index=annotated.index)).fillna("").astype(str).ne("").sum()),
        "hpo_nonempty_row_count": int(annotated.get("hpo_terms", pd.Series("", index=annotated.index)).fillna("").astype(str).ne("").sum()),
        "annotation_changes_event_classification": False,
        "annotation_changes_filtering": False,
        "cnvpro_filter_logic_imported": False,
    }


def main():
    args = parse_args()
    logger = setup_logger("cnv_event_annotation", args.log or None)
    if args.build_bundle:
        summary = build_annotation_bundle(
            bundle_db=args.bundle_db,
            cytoband_source=args.cytoband_source,
            refgene_bed=args.refgene_bed,
            omim_sources=args.omim_source,
            hpo_sources=args.hpo_source,
            source_manifest=args.source_manifest,
            source_root=args.source_root,
            pgta_source_root=args.pgta_source_root,
            bundle_id=args.bundle_id,
            genome_build=args.genome_build,
        )
        if args.output_json:
            ensure_parent(args.output_json).write_text(
                json.dumps(summary, indent=2, sort_keys=True),
                encoding="utf-8",
            )
        logger.info(
            "built annotation bundle %s gene_rows=%s omim_rows=%s hpo_status=%s",
            args.bundle_db,
            summary["gene_row_count"],
            summary["omim_row_count"],
            summary["hpo_source_status"],
        )
        return 0
    if not args.output_tsv or not args.output_json:
        raise ValueError("--output-tsv and --output-json are required unless --build-bundle is used")
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
