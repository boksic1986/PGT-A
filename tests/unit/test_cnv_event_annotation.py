import gzip
import sqlite3
from pathlib import Path

import pandas as pd

from pgta.predict.annotation import (
    annotate_events_dataframe,
    build_annotation_bundle,
    collect_events,
    fetch_metadata,
    normalize_chromosome,
)


def _write_test_bundle(path: Path) -> None:
    conn = sqlite3.connect(path)
    try:
        conn.execute(
            "CREATE TABLE cytoband (chrom TEXT, start INTEGER, end INTEGER, name TEXT, stain TEXT)"
        )
        conn.execute(
            "CREATE TABLE gene (chrom TEXT, start INTEGER, end INTEGER, gene_symbol TEXT, gene_location TEXT)"
        )
        conn.execute(
            "CREATE TABLE omim_gene (gene_symbol TEXT, omim_id TEXT, phenotype TEXT)"
        )
        conn.execute(
            "CREATE TABLE hpo_term (gene_symbol TEXT, hpo_id TEXT, hpo_name TEXT)"
        )
        conn.execute(
            "CREATE TABLE region_context (chrom TEXT, start INTEGER, end INTEGER, context_label TEXT, context_reason TEXT)"
        )
        conn.execute(
            "CREATE TABLE bundle_metadata (key TEXT PRIMARY KEY, value TEXT)"
        )
        conn.execute("INSERT INTO cytoband VALUES ('chr1', 100, 300, 'p36.33', 'gneg')")
        conn.execute("INSERT INTO gene VALUES ('chr1', 120, 180, 'GENE1', 'overlap')")
        conn.execute("INSERT INTO omim_gene VALUES ('GENE1', '123456', 'example phenotype')")
        conn.execute("INSERT INTO hpo_term VALUES ('GENE1', 'HP:0000001', 'example hpo')")
        conn.execute(
            "INSERT INTO region_context VALUES ('chr1', 90, 310, 'test_context', 'unit test')"
        )
        conn.execute("INSERT INTO bundle_metadata VALUES ('bundle_id', 'unit-test-bundle')")
        conn.execute("INSERT INTO bundle_metadata VALUES ('genome_build', 'hg19')")
        conn.commit()
    finally:
        conn.close()


def test_normalize_chromosome_accepts_common_aliases():
    assert normalize_chromosome("1") == "chr1"
    assert normalize_chromosome("chr1") == "chr1"
    assert normalize_chromosome("23") == "chrX"
    assert normalize_chromosome("24") == "chrY"
    assert normalize_chromosome("X") == "chrX"
    assert normalize_chromosome("Y") == "chrY"


def test_annotation_bundle_overlaps_cytoband_gene_omim_hpo(tmp_path):
    bundle = tmp_path / "pgta_cnv_annotation.hg19.sqlite"
    _write_test_bundle(bundle)
    events = pd.DataFrame(
        [
            {
                "event_id": "E1",
                "sample_id": "S1",
                "chrom": "1",
                "start": 150,
                "end": 250,
                "state": "gain",
                "event_layer": "autosomal_report",
            }
        ]
    )

    annotated = annotate_events_dataframe(events, bundle_db=bundle, genome_build="hg19")

    assert len(annotated) == 1
    row = annotated.iloc[0]
    assert row["chrom"] == "chr1"
    assert row["cytoband"] == "p36.33"
    assert row["genes"] == "GENE1"
    assert row["gene_number"] == 1
    assert row["omim_genes"] == "GENE1:123456"
    assert "example phenotype" in row["omim_phenotypes"]
    assert "HP:0000001" in row["hpo_terms"]
    assert row["region_context"] == "test_context"
    assert row["annotation_backend"] == "pgta_sqlite"
    assert row["annotation_status"] == "annotated"
    assert row["annotation_bundle_id"] == "unit-test-bundle"


def test_missing_annotation_bundle_preserves_rows_with_missing_status(tmp_path):
    missing_bundle = tmp_path / "missing.sqlite"
    events = pd.DataFrame(
        [
            {
                "event_id": "E1",
                "sample_id": "S1",
                "chrom": "chr2",
                "start": 1000,
                "end": 2000,
                "state": "loss",
                "event_layer": "autosomal_report",
            }
        ]
    )

    annotated = annotate_events_dataframe(events, bundle_db=missing_bundle, genome_build="hg19")

    assert len(annotated) == 1
    row = annotated.iloc[0]
    assert row["annotation_status"] == "missing_backend"
    assert row["annotation_backend"] == "pgta_sqlite"
    assert row["genes"] == ""
    assert row["cytoband"] == ""


def test_collect_events_deduplicates_report_events_against_classifier(tmp_path):
    report_path = tmp_path / "report_events.tsv"
    classifier_path = tmp_path / "sample.candidate_classification.tsv"
    base_row = {
        "sample_id": "S1",
        "candidate_id": "C1",
        "event_id": "E1",
        "chrom": "chr1",
        "start": 100,
        "end": 200,
        "state": "gain",
    }
    pd.DataFrame([{**base_row, "v2_report_layer_class": "report_event"}]).to_csv(
        report_path, sep="\t", index=False
    )
    pd.DataFrame([{**base_row, "v2_report_layer_class": "report_event"}]).to_csv(
        classifier_path, sep="\t", index=False
    )

    events = collect_events(event_paths=[report_path], classification_paths=[classifier_path])

    assert len(events) == 1
    assert events.iloc[0]["event_layer"] == "autosomal_report"


def test_build_annotation_bundle_from_refgene_and_omim_sources_records_missing_hpo(tmp_path):
    cytoband = tmp_path / "cytoBand.txt.gz"
    with gzip.open(cytoband, "wt", encoding="utf-8") as handle:
        handle.write("chr1\t0\t300\tp36.33\tgneg\n")

    refgene = tmp_path / "refGene.sorted.bed"
    refgene.write_text(
        "\n".join(
            [
                "1\t100\t260\t+\tGENE1\tNM_000001\t120\t250\t100,\t260,",
                "1\t1000\t1300\t-\tGENE2\tNM_000002\t1020\t1280\t1000,\t1300,",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    omim = tmp_path / "20181210_OMIMannotations.tsv.gz"
    with gzip.open(omim, "wt", encoding="utf-8") as handle:
        handle.write("genes\tMim Number\tPhenotypes\tInheritance\n")
        handle.write("GENE1\t123456\texample phenotype\tAD\n")

    bundle = tmp_path / "pgta_cnv_annotation.hg19.sqlite"
    manifest = tmp_path / "source_manifest.tsv"

    summary = build_annotation_bundle(
        bundle_db=bundle,
        cytoband_source=cytoband,
        refgene_bed=refgene,
        omim_sources=[omim],
        hpo_sources=[],
        source_manifest=manifest,
        source_root=tmp_path,
        pgta_source_root=tmp_path,
        bundle_id="unit-test-gene-omim",
        genome_build="hg19",
    )

    assert summary["gene_row_count"] == 2
    assert summary["omim_row_count"] == 1
    assert summary["hpo_row_count"] == 0
    assert summary["hpo_source_status"] == "missing_hpo_source"
    assert manifest.exists()
    manifest_text = manifest.read_text(encoding="utf-8")
    assert "refgene" in manifest_text
    assert "omim" in manifest_text

    conn = sqlite3.connect(bundle)
    try:
        metadata = fetch_metadata(conn)
    finally:
        conn.close()
    assert metadata["bundle_id"] == "unit-test-gene-omim"
    assert metadata["gene_source_status"] == "ready"
    assert metadata["omim_source_status"] == "ready"
    assert metadata["hpo_source_status"] == "missing_hpo_source"

    annotated = annotate_events_dataframe(
        pd.DataFrame(
            [
                {
                    "event_id": "E1",
                    "sample_id": "S1",
                    "chrom": "chr1",
                    "start": 150,
                    "end": 240,
                    "state": "gain",
                    "event_layer": "autosomal_report",
                }
            ]
        ),
        bundle_db=bundle,
        genome_build="hg19",
    )
    row = annotated.iloc[0]
    assert row["genes"] == "GENE1"
    assert row["gene_number"] == 1
    assert row["gene_location"] == "refGene"
    assert row["omim_genes"] == "GENE1:123456"
    assert row["omim_phenotypes"] == "example phenotype"
    assert row["gene_source_status"] == "ready"
    assert row["omim_source_status"] == "ready"
    assert row["hpo_source_status"] == "missing_hpo_source"
