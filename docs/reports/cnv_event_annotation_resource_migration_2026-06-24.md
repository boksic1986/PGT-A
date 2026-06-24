# PGTA CNV Event Annotation Resource Migration

Date: 2026-06-24

## Scope

This loop implements a PGTA-owned event-level annotation sidecar. It does not
import CNVpro calling, cghFLasso, CN thresholds, `filterCNV.py` filtering, or
AnnotSV ranking decisions.

## Resource Audit

Source checked on `fengxian`:

- `/home/jiucheng/project/cnvpro/cnvpro_minimal_20260614`

Annotation-like files found in the minimal CNVpro package:

- `code/cnvseqpipe/cytoBand.txt.gz`
- `code/cnvseqpipe/cytoBand.hg38.txt.gz`
- AnnotSV 2.2 Tcl/program files under `resources/utility/AnnotSV_2.2`

The minimal CNVpro package did not expose standalone hg19 gene, OMIM, or HPO
tables suitable for direct PGTA migration. PGTA therefore migrated only the
confirmed hg19 cytoband resource in this first pass.

## PGTA-Owned Bundle

PGTA resource directory:

- `/data/project/CNV/PGT-A/resources/annotation/hg19/`

Migrated source:

- `/data/project/CNV/PGT-A/resources/annotation/hg19/source_from_cnvpro/cytoBand.txt.gz`
- MD5: `161003a509f74247805d8ff346a8149e`

Generated files:

- `/data/project/CNV/PGT-A/resources/annotation/hg19/source_manifest.tsv`
- `/data/project/CNV/PGT-A/resources/annotation/hg19/pgta_cnv_annotation.hg19.sqlite`

SQLite bundle metadata:

- `bundle_id`: `pgta_cnv_annotation.hg19.v20260624`
- `genome_build`: `hg19`
- cytoband rows: `862`
- gene/OMIM/HPO status:
  `not_present_in_cnvpro_minimal_bundle_pending_independent_pgta_resource`

Tables created:

- `cytoband`
- `gene`
- `omim_gene`
- `hpo_term`
- `region_context`
- `bundle_metadata`

The gene, OMIM, HPO, and region context tables are present but currently empty.
They are placeholders for PGTA-owned resources; they are not populated from
Arraytools and are not silently shared with CNVpro.

## Workflow Integration

New formal sidecar target:

- `cnv_event_annotation`

New production rule:

- `cnv_event_annotation_build`

Inputs:

- Branch B V2 report events
- Branch B V2 per-sample classifier rows for internal/review context
- Branch S evidence rows for sex-chromosome review context
- PGTA-owned SQLite bundle

Outputs:

- `wisecondorx/cnv/postprocess_gap2m/event_annotation/cnv_event_annotation.tsv`
- `wisecondorx/cnv/postprocess_gap2m/event_annotation/cnv_event_annotation.json`

`cnv_report` now consumes the annotation TSV and displays cytoband/gene fields
where available. Annotation status and bundle identity are written into
`report_contract`.

## Runtime Contract

Annotation is display-only:

- does not change event classification
- does not change Branch B V2 filtering
- does not change Branch S state
- does not change PASS/FAIL or TP/FN/FP
- does not pull `branch_b_v2_report_ablation` into production `cnv_report`

If the PGTA bundle is absent, annotation outputs `annotation_status=missing_backend`
and preserves input event rows; report generation is not blocked.

## Validation

Remote pytest on `fengxian`:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_cnv_event_annotation.py \
  tests/unit/test_cnv_report.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py -q

43 passed in 1.87s
```

Remote dry-run on Y/H/G/0615 active configs:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile <active_config> --cores 1 -n cnv_event_annotation cnv_report
```

All four active configs passed dry-run. The planned DAG adds:

```text
cnv_event_annotation_build -> cnv_event_annotation -> cnv_report_summary -> cnv_report
```

0615 smoke materialization:

- annotation rows: `147`
- samples: `5`
- cytoband nonempty rows: `147`
- gene nonempty rows: `0`
- report contract annotation backend: `pgta_sqlite`
- report contract bundle ID: `pgta_cnv_annotation.hg19.v20260624`

## Remaining Work

1. Add PGTA-owned gene/OMIM/HPO resources from an approved independent source.
2. Populate `gene`, `omim_gene`, `hpo_term`, and `region_context` tables.
3. Add coverage tests requiring non-empty gene annotation only after those
   resources are present.
4. If AnnotSV rank is evaluated later, keep it as a shadow comparison with a
   separate validation gate.
