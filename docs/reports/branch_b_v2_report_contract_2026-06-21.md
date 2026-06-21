# Branch B V2 Report Contract Integration

date: 2026-06-21
status: materialized_development_only
active_reference_id: h_r0_shadow_ref_20260619
branch_a_input: explicit_merge_gap_bp_2000000_overlay

## Scope

This loop connects the already materialized Branch B V2 burden stratification
outputs to `cnv_report` as development-only report evidence.

It does not change Branch A, does not add hard filters, does not promote
Branch B V2, and does not promote Branch S or SCA final calling.

## Contract

`cnv_report` now has an explicit Branch B V2 burden-display contract. The
report may read:

- `wisecondorx/cnv/postprocess_gap2m/v2_benchmark/summary.json`
- `wisecondorx/cnv/postprocess_gap2m/v2_benchmark/sample_summary.tsv`

Report-level fields:

- `branch_b_v2_burden_status`
- `branch_b_v2_background_unknown_review_count`
- `branch_b_v2_branch_s_review_count`
- `branch_b_v2_technical_risk_review_count`
- `branch_b_v2_report_candidate_count`
- `branch_b_v2_legacy_fields_used=false`
- `branch_b_v2_final_impact=development_review_only`

Sample-level TSV/JSON/MD/HTML output displays the same V2 burden context per
sample when the benchmark outputs are configured as `cnv_report` inputs.

## Interpretation

Branch B V2 burden display is audit/report context only:

- it is not FP-reduction proof;
- it is not final promotion;
- it does not make `UNKNOWN_BACKGROUND` benign;
- it does not turn length, GC/RC, clean support, or B-side signal context into
  hard filters;
- it does not use legacy/current-code Branch B decision fields.

The report package remains `development_only_not_final_release` until the
upstream Branch A/B/S gates are explicitly promoted.

## Validation Status

Remote validation was run on:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419
```

Unit tests:

```bash
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_cnv_report.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py \
  tests/unit/test_current_context_index.py -q
```

Result:

```text
28 passed in 0.93s
```

Snakemake dry-run:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake \
  -s Snakefile --configfile <active gap2m config> --cores 1 -n \
  branch_b_v2_benchmark branch_s_review cnv_report
```

Passed for:

- `config_predict_y_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`
- `config_predict_h_20260608_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`
- `config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`

Materialized `cnv_report` for all three active gap2m configs.

| cohort | report contract | V2 burden status | V2 truth preserved | FN | hard-suppressed truth | key display check |
|---|---|---|---:|---:|---:|---|
| Y1-Y8 | `development_only_not_final_release` | `ready` | 10/10 | 0 | 0 | V2 burden counts displayed |
| H1-H16 | `development_only_not_final_release` | `ready` | 10/10 | 0 | 0 | H6 chr21 visible |
| 2026-06-15 | `development_only_not_final_release` | `skipped_no_truth` | 0/0 | 0 | 0 | burden/context only |

Report JSON contract checks:

- `reference_id=h_r0_shadow_ref_20260619`
- `branch_b_v2_final_impact=development_review_only`
- `branch_b_v2_legacy_fields_used=false`
- Y1-Y8 burden counts: background unknown review `84`, Branch S review `13`,
  technical risk review `0`, report candidate `0`
- H1-H16 burden counts: background unknown review `63`, Branch S review `42`,
  technical risk review `0`, report candidate `0`
- 2026-06-15 burden counts: background unknown review `151`, Branch S review
  `14`, technical risk review `0`, report candidate `0`
