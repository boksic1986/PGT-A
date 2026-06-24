# CNV Event Annotation Gene/OMIM/HPO Handoff

Date: 2026-06-25

## Current Context

Active branch:

```text
codex/pgta-annotation-gene-omim-hpo-audit
```

Current report:

```text
docs/reports/cnv_event_annotation_gene_omim_hpo_bundle_2026-06-25.md
```

This loop upgrades the PGTA-owned hg19 annotation bundle from cytoband-only to
gene/OMIM-ready. HPO remains explicitly marked as missing.

## Files Changed

Code and tests:

- `pgta/predict/annotation.py`
- `pgta/predict/report.py`
- `tests/unit/test_cnv_event_annotation.py`
- `tests/unit/test_cnv_report.py`

Docs:

- `docs/reports/cnv_event_annotation_gene_omim_hpo_bundle_2026-06-25.md`
- `docs/handoff/2026-06-25_0005_cnv_event_annotation_gene_omim_hpo_handoff.md`
- `docs/CURRENT_CONTEXT_INDEX.md`
- `CURRENT_STATE.md`
- `PLANS.md`

Remote PGTA resources, not git-tracked:

- `/data/project/CNV/PGT-A/resources/annotation/hg19/source_from_cnvpro/`
- `/data/project/CNV/PGT-A/resources/annotation/hg19/source_manifest.tsv`
- `/data/project/CNV/PGT-A/resources/annotation/hg19/pgta_cnv_annotation.hg19.sqlite`
- `/data/project/CNV/PGT-A/resources/annotation/hg19/pgta_cnv_annotation.hg19.build_summary.json`

## Resource State

Runtime bundle:

```text
/data/project/CNV/PGT-A/resources/annotation/hg19/pgta_cnv_annotation.hg19.sqlite
```

Bundle ID:

```text
pgta_cnv_annotation.hg19.v20260625_gene_omim
```

Bundle row counts:

- cytoband: 862
- gene: 70384
- OMIM: 13736
- HPO: 0

Source status:

- gene: `ready`
- OMIM: `ready`
- HPO: `missing_hpo_source`

## Runtime Boundary

PGTA runtime reads only the PGTA-owned SQLite bundle. It does not read CNVpro
source paths or Arraytools databases at runtime.

Do not import the following into PGTA decision logic:

- CNVpro `filterCNV.py`
- CNVpro CN thresholds
- cghFLasso / `CNVcalling.R`
- AnnotSV rank

Annotation remains display-only and must not change Branch A/B/S, filtering,
report-event counts, PASS/FAIL, or TP/FN/FP.

## Validation

Remote unit tests on `fengxian`:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_cnv_event_annotation.py \
  tests/unit/test_cnv_report.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py \
  tests/unit/test_current_context_index.py -q

48 passed in 1.42s
```

Remote dry-run:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile <Y/H/G/0615 active config> --cores 1 -n \
  cnv_event_annotation cnv_report
```

All four active configs passed dry-run.

0615 forced smoke run:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml \
  --cores 4 --forcerun cnv_event_annotation_build cnv_event_annotation cnv_report
```

Smoke output:

- annotation rows: 147
- cytoband nonempty rows: 147
- gene nonempty rows: 147
- OMIM nonempty rows: 146
- HPO nonempty rows: 0
- report contract source status:
  `gene=ready`, `omim=ready`, `hpo=missing_hpo_source`

## Next Recommended Step

Commit, create PR, review, merge to `main`, and sync tracked `main` content to
the non-git remote mirror. Keep SQLite/resources as remote PGTA resources and
out of git.

After merge, rerun the targeted annotation/report pytest and at least one
active `cnv_event_annotation cnv_report -n` dry-run on the synced mirror.

