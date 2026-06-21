# Handoff: Report Main Convergence And CNV Plot

Date: 2026-06-22 04:37

## Context Source

Active previous handoff:
`docs/handoff/2026-06-22_0930_lowres_branch_bs_integration_handoff.md`

New report:
`docs/reports/report_main_convergence_cnv_plot_2026-06-22.md`

## Scope Completed

This loop changed only report-layer behavior, Branch B V2 report visibility, and
CNV plot output.

Unchanged:

- Branch A WisecondorX/CBS caller
- active 1Mb reference `h_r0_shadow_ref_20260619`
- reference build
- BAM/mapping
- Branch S production status
- default `merge_gap_bp=0`

## Code Changes

Files changed:

- `pgta/predict/branch_b/plot.py`
- `pgta/predict/branch_b/v2_classifier.py`
- `pgta/predict/branch_b/v2_benchmark.py`
- `pgta/predict/report.py`
- `rules/predict_layout.smk`
- `rules/predict_workflow.smk`
- `tests/unit/test_branch_b_plot.py`
- `tests/unit/test_branch_b_v2_classifier.py`
- `tests/unit/test_branch_b_v2_benchmark.py`
- `tests/unit/test_cnv_report.py`
- `tests/unit/test_branch_ab_phase12_workflow_contract.py`

Main behavior:

- Branch B V2 emits report visibility:
  `report_strong_event`, `report_weak_event`, `internal_review_event`,
  `filtered_event`, `branch_s_event`.
- V2 benchmark writes report events as final autosomal main-table event input.
- `cnv_report` consumes V2 report events and V2 sample summary when V2 benchmark
  is available.
- `cnv_report` keeps zero-report-event samples visible in `cnv_summary.tsv`.
- CNV plot uses only `calibrated_z`.
- Per-sample plot bins are written as `{sample}.plot_bins.tsv`.
- SVG highlights only final autosomal report events as `dup` or `del`.

## Remote Commands Run

Remote path:
`/data/project/CNV/PGT-A/refactor_validation_20260419`

Unit tests:

```bash
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_cnv_report.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py \
  tests/unit/test_branch_b_plot.py \
  tests/unit/test_branch_b_v2_classifier.py \
  tests/unit/test_branch_b_v2_benchmark.py -q
```

Result:

```text
79 passed
```

Dry-run:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake \
  -s Snakefile \
  --configfile <active_lowres_gap2m_config.yaml> \
  --cores 1 -n \
  branch_b_v2_benchmark branch_s_review cnv_report
```

All four active configs returned `STATUS 0`.

Report refresh:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake \
  -s Snakefile \
  --configfile <active_lowres_gap2m_config.yaml> \
  --cores 4 \
  --forcerun cnv_report_summary \
  cnv_report
```

All four active configs completed.

## Materialized Outputs

Primary report outputs:

- `results_y_h_r0_shadow_ref_20260619/wisecondorx/cnv/report/cnv_summary.*`
- `results_h_20260608_h_r0_shadow_ref_20260619/wisecondorx/cnv/report/cnv_summary.*`
- `results_g1_g8_h_r0_shadow_ref_20260621/wisecondorx/cnv/report/cnv_summary.*`
- `results_20260615_h_r0_shadow_ref_20260619/wisecondorx/cnv/report/cnv_summary.*`

Plot outputs:

- `wisecondorx/cnv/plots/{sample}.final_cnv.svg`
- `wisecondorx/cnv/plots/{sample}.plot_bins.tsv`

## Acceptance Summary

| cohort | samples | truth | preserved | FN | truth filtered | report | strong | weak | internal | filtered | Branch S | plots |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Y1-Y8 | 8 | 10 | 10 | 0 | 0 | 40 | 21 | 19 | 31 | 13 | 13 | 8/8 |
| H1-H16 | 16 | 10 | 10 | 0 | 0 | 23 | 6 | 17 | 32 | 8 | 42 | 16/16 |
| G1-G8 | 8 | 10 | 10 | 0 | 0 | 26 | 15 | 11 | 29 | 7 | 13 | 8/8 |
| 2026-06-15 | 5 | 0 | 0 | 0 | 0 | 71 | 52 | 19 | 57 | 23 | 14 | 5/5 |

Specific checks:

- H6 chr21 gain remains visible as `report_weak_event`.
- H6 chrX remains Branch S.
- G2 chr8 truth remains visible as `internal_review_event`, not filtered.
- G2 chr11 truth remains `report_strong_event`.
- 2026-06-15 remains burden/context only.

Plot contract:

- required columns present in every plot TSV
- `report_state` limited to `dup`, `del`, `neutral`
- SVG legend limited to `dup`, `del`, `neutral bin`, `smooth z trend`
- no Branch A/B, internal, filtered, mask, or rejected plot legend terms

## Important Implementation Note

`cnv_report` is an aggregate target. Forcing only `cnv_report` does not rerun
`report.py`. To refresh report files after code changes, force
`cnv_report_summary`:

```bash
snakemake --forcerun cnv_report_summary cnv_report
```

## Current Limitations

- This remains development-only.
- Branch B V2 and Branch S are not production-final.
- The shadow reference is not production-promoted.
- 2026-06-15 has no locked truth labels, so no TP/FN/FP conclusion is allowed.
- Report event burden remains high in 2026-06-15.

## Next Recommended Gate

Inspect remaining `report_strong_event` and `report_weak_event` rows by
candidate-level evidence, especially in 2026-06-15:

- lowres same-direction evidence
- ref-MAD stability
- clean support
- high-risk region context
- B-side signal support
- SCA/sex-chromosome separation

Do not use sample-level event counts or 2026-06-15 burden counts to reverse
engineer new filters. Any new demotion/filter rule must preserve Y/H/G truth
with FN=0 and keep H6 chr21 and G2 truth visible.
