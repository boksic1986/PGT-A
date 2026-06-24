# PGTA CNV Event Annotation Gene/OMIM/HPO Bundle

Date: 2026-06-25

## Scope

This loop upgrades the PGTA-owned hg19 event annotation bundle from
`cytoband-ready` to `gene/OMIM-ready`. HPO remains explicitly unavailable
because no independent HPO source was present in the audited CNVpro minimal
AnnotSV GRCh37 resources.

This is annotation/report display work only. It does not modify Branch A,
Branch B V2, Branch S, reference build, sex calling, report-event
classification, filtering, PASS/FAIL, or TP/FN/FP.

## Resource Audit

Audited CNVpro minimal source path on `fengxian`:

```text
/home/jiucheng/project/cnvpro/cnvpro_minimal_20260614/resources/utility/AnnotSV_2.2/share/doc/AnnotSV/Annotations_Human/
```

Resources copied into the PGTA-owned resource area:

```text
/data/project/CNV/PGT-A/resources/annotation/hg19/source_from_cnvpro/
```

Included resources:

| resource_type | source | md5 | size_bytes |
|---|---|---|---:|
| cytoband | `code/cnvseqpipe/cytoBand.txt.gz` | `161003a509f74247805d8ff346a8149e` | 6261 |
| refgene | `RefGene/GRCh37/refGene.sorted.bed` | `dfbc8fcd3cb4f9b393409096dd9d075f` | 17241605 |
| omim | `Genes-based/OMIM/20181126_morbidGenes.tsv.gz` | `53413557c88221c329c029a03f01c52c` | 29199 |
| omim | `Genes-based/OMIM/20181126_morbidGenesCandidates.tsv.gz` | `e5a1146cc82087e7691956c66721aa82` | 9138 |
| omim | `Genes-based/OMIM/20181210_OMIMannotations.tsv.gz` | `c2c82a489301286b6a396a4e5a3dbce9` | 383517 |

No standalone HPO source was found in the audited resource subtree. The bundle
therefore records `hpo_source_status=missing_hpo_source` and keeps
`hpo_term` empty.

## PGTA-Owned Bundle

Runtime bundle:

```text
/data/project/CNV/PGT-A/resources/annotation/hg19/pgta_cnv_annotation.hg19.sqlite
```

Bundle metadata after rebuild:

| field | value |
|---|---|
| `bundle_id` | `pgta_cnv_annotation.hg19.v20260625_gene_omim` |
| `genome_build` | `hg19` |
| `cytoband_row_count` | 862 |
| `gene_row_count` | 70384 |
| `omim_row_count` | 13736 |
| `hpo_row_count` | 0 |
| `gene_source_status` | `ready` |
| `omim_source_status` | `ready` |
| `hpo_source_status` | `missing_hpo_source` |

The runtime path is PGTA-owned. The workflow does not read CNVpro source paths
or Arraytools databases at runtime.

## Code Changes

Updated modules:

- `pgta/predict/annotation.py`
- `pgta/predict/report.py`
- `tests/unit/test_cnv_event_annotation.py`
- `tests/unit/test_cnv_report.py`

`pgta.predict.annotation` now supports a bundle-build action in the same module
as the annotation runtime. This avoids adding a parallel one-off script.

The builder populates:

- `cytoband`
- `gene`
- `omim_gene`
- `hpo_term`
- `region_context`
- `bundle_metadata`

`cnv_event_annotation` now reports source status columns per event:

- `gene_source_status`
- `omim_source_status`
- `hpo_source_status`

`cnv_report` writes the same source statuses into `report_contract`.

## Validation

TDD red check on remote mirror:

```text
ImportError: cannot import name 'build_annotation_bundle' from 'pgta.predict.annotation'
```

Remote unit tests after implementation:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_cnv_event_annotation.py \
  tests/unit/test_cnv_report.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py \
  tests/unit/test_current_context_index.py -q

48 passed in 1.42s
```

Remote dry-run for active Y/H/G/0615 configs:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile <active_config> --cores 1 -n cnv_event_annotation cnv_report
```

All four active configs passed dry-run.

0615 forced smoke materialization:

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
- report contract bundle ID:
  `pgta_cnv_annotation.hg19.v20260625_gene_omim`
- report contract source status:
  `gene=ready`, `omim=ready`, `hpo=missing_hpo_source`

## Boundary

Explicitly not imported:

- CNVpro `filterCNV.py`
- CNVpro CN thresholds
- cghFLasso / `CNVcalling.R`
- AnnotSV rank as a report decision rule
- Arraytools database paths

Annotation remains display-only and does not change event counts or report
classification.

## Known Limitation

The SQLite bundle is a remote PGTA resource and is not tracked in git. Because
the bundle path is passed as a workflow parameter to preserve missing-backend
fallback behavior, resource rebuilds should use:

```text
--forcerun cnv_event_annotation_build
```

when validating a newly built bundle against already materialized annotation
outputs.

