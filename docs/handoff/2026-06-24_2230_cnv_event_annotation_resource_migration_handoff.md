# Handoff: CNV Event Annotation Resource Migration

Date: 2026-06-24 22:30

## Context

Active task: migrate useful CNVpro annotation resources into a PGTA-owned event
annotation sidecar without sharing Arraytools databases or importing CNVpro
calling/filter logic.

Adopted previous handoff:

- `docs/handoff/2026-06-24_0031_sample_master_manifest_handoff.md`

## Completed

- Added `pgta.predict.annotation` as the PGTA event annotation sidecar.
- Added dispatcher action `cnv_event_annotation`.
- Added Snakemake target/rule `cnv_event_annotation` /
  `cnv_event_annotation_build`.
- Added event annotation outputs under
  `wisecondorx/cnv/postprocess_gap2m/event_annotation/`.
- Added `cnv_report --event-annotation-tsv` integration.
- `cnv_report` now writes annotation status, backend, row count, and bundle ID
  into `report_contract`.
- `cnv_report` top event display can include cytoband/gene fields when present.
- Added unit tests for annotation overlap, missing backend fallback, workflow
  contract, and report merge behavior.

## Remote Resource State

PGTA-owned resource directory:

- `/data/project/CNV/PGT-A/resources/annotation/hg19/`

Migrated CNVpro minimal resource:

- `source_from_cnvpro/cytoBand.txt.gz`
- MD5: `161003a509f74247805d8ff346a8149e`

Built bundle:

- `/data/project/CNV/PGT-A/resources/annotation/hg19/pgta_cnv_annotation.hg19.sqlite`
- bundle ID: `pgta_cnv_annotation.hg19.v20260624`
- cytoband rows: `862`

The minimal CNVpro bundle did not expose standalone gene/OMIM/HPO annotation
tables. PGTA-created `gene`, `omim_gene`, `hpo_term`, and `region_context`
tables exist in SQLite but are empty placeholders.

## Validation

Remote unit tests on `fengxian`:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_cnv_event_annotation.py \
  tests/unit/test_cnv_report.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py -q

43 passed in 1.87s
```

Remote dry-run:

- Y/H/G/0615 active configs all passed:
  `cnv_event_annotation cnv_report -n`

0615 smoke run:

- command target: `cnv_event_annotation`
- annotation rows: `147`
- cytoband nonempty rows: `147`
- gene nonempty rows: `0`
- report contract bundle ID:
  `pgta_cnv_annotation.hg19.v20260624`

## Important Boundaries

- Annotation does not alter Branch A, Branch B V2, Branch S, reference build,
  sex calling, report-event classification, filtering, PASS/FAIL, or TP/FN/FP.
- PGTA runtime reads only PGTA-owned bundle paths.
- CNVpro `filterCNV.py`, cghFLasso, CN thresholds, and AnnotSV rank were not
  migrated into decision logic.
- If the SQLite bundle is missing, annotation outputs `missing_backend` and
  report generation continues.

## Files Modified

- `Snakefile`
- `rules/pipeline_modes.smk`
- `rules/predict_layout.smk`
- `rules/predict_workflow.smk`
- `rules/script_entrypoints.smk`
- `rules/target_assembly.smk`
- `scripts/_compat_entry.py`
- `pgta/predict/annotation.py`
- `pgta/predict/report.py`
- `tests/unit/test_cnv_event_annotation.py`
- `tests/unit/test_cnv_report.py`
- `tests/unit/test_branch_ab_phase12_workflow_contract.py`
- `docs/reports/cnv_event_annotation_resource_migration_2026-06-24.md`
- `docs/CURRENT_CONTEXT_INDEX.md`
- `CURRENT_STATE.md`
- `PLANS.md`

## Next Step

Add PGTA-owned gene/OMIM/HPO data from an approved independent source, populate
the SQLite tables, and then add coverage tests requiring non-empty gene
annotation. Do not use Arraytools database paths and do not promote AnnotSV rank
without a separate validation gate.
