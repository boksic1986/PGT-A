# Branch B V2 Sex-Route Refinement

Date: 2026-06-21

Status: `materialized_shadow_refinement`

## Scope

This report records the first Branch B V2 refinement after the gap2m benchmark:
sex-chromosome candidates are routed out of the autosomal Branch B positive
support class and into Branch S review.

The active contract is unchanged:

```text
reference_id=h_r0_shadow_ref_20260619
branch_a.merge_gap_bp=2_000_000
branch_a.strong_z=10.0
postprocess.output_dir=wisecondorx/cnv/postprocess_gap2m
final_report_impact=none_shadow_only
```

This is not final SCA promotion. It does not change WisecondorX predict,
Branch A candidate generation, mosaic logic, sex calling, or final report
release status.

## Classifier Change

Before applying autosomal Branch B V2 classes, the classifier now checks the
candidate chromosome:

```text
chrX / X / chrY / Y
  -> v2_candidate_class=V2_SEX_CHROMOSOME_REVIEW
  -> v2_classifier_action=V2_ROUTE_BRANCH_S_REVIEW
```

The original evidence tier, evidence gate, review priority, and
`none_shadow_only` final-report impact are preserved. This keeps chrX/chrY
signals visible for Branch S/SCA review while preventing them from inflating the
autosomal Branch B positive-support burden.

## Remote Materialized Results

Remote mirror:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419
```

Per-sample classifier outputs:

```text
results_y_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m/v2_classifier/*.candidate_classification.tsv
results_h_20260608_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m/v2_classifier/*.candidate_classification.tsv
results_20260615_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m/v2_classifier/*.candidate_classification.tsv
```

Benchmark summaries:

```text
results_y_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m/v2_benchmark/summary.json
results_h_20260608_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m/v2_benchmark/summary.json
results_20260615_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m/v2_benchmark/summary.json
```

## Class Burden After Sex Routing

| cohort | candidates | V2 positive-support class | V2 sex-chromosome review | V2 no-call contract risk | truth preserved | FN | hard-suppressed truth |
|---|---:|---:|---:|---:|---:|---:|---:|
| Y1-Y8 | 97 | 78 | 13 | 6 | 10/10 | 0 | 0 |
| H1-H16 | 105 | 57 | 42 | 6 | 10/10 | 0 | 0 |
| H7-H16 subset | 57 | 24 | 31 | 2 | no truth | NA | NA |
| 2026-06-15 | 165 | 127 | 14 | 24 | no truth | NA | NA |

Important distinction:

- `v2_positive_support_candidate_count` in existing benchmark summary JSON is
  evidence-tier based and can still include sex-chromosome positive-support
  evidence.
- The table above uses classifier classes from
  `v2_classifier/*.candidate_classification.tsv`, so it separates autosomal
  Branch B review from Branch S/SCA review.

## Truth Preservation Notes

Y1-Y8:

- Locked truth remains preserved 10/10 with FN=0.
- Y3 chrX loss is now `V2_SEX_CHROMOSOME_REVIEW` /
  `V2_ROUTE_BRANCH_S_REVIEW`.

H1-H6 within H1-H16:

- Locked truth remains preserved 10/10 with FN=0.
- H5/H6 chrX loss truth events are now `V2_SEX_CHROMOSOME_REVIEW` /
  `V2_ROUTE_BRANCH_S_REVIEW`.
- H6 chr21 gain remains autosomal `V2_POSITIVE_SUPPORT_REVIEW` with
  `top_a_abs_zscore=7.1135`.

## Interpretation

This refinement supports:

- Branch B V2 can separate autosomal CNV review evidence from sex-chromosome
  review routing without losing locked truth.
- The H7-H16 context burden is clearer: 31/57 candidate rows are sex-chromosome
  review, while 24/57 remain autosomal positive-support review.
- The 2026-06-15 five-sample context remains no-truth burden/context only.

This refinement does not support:

- final SCA promotion;
- 2026-06-15 TP/FN/FP conclusions;
- promotion of `merge_gap_bp=2_000_000` to default;
- promotion of `h_r0_shadow_ref_20260619` to production;
- Branch B V2 final-report decisions.

## Remote Commands

Unit tests:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 /biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest tests/unit/test_branch_b_v2_classifier.py tests/unit/test_branch_b_v2_benchmark.py tests/unit/test_branch_ab_phase12_workflow_contract.py -q'
```

Result:

```text
28 passed in 0.95s
```

Forced materialization:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile <branch_b_v2_gap2m_config> --cores 24 --forcerun cnv_branch_b_v2_classifier cnv_branch_b_v2_benchmark branch_b_v2_benchmark
```

All three configs completed:

- `config_predict_y_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`
- `config_predict_h_20260608_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`
- `config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`

## Next Gate

The next gate is still Branch B V2 review-burden refinement plus Branch S
review-reportable development. The immediate next useful step is to decide
which Branch S summary fields should be carried into the next workflow-generated
P6/report package as review labels, still with explicit uncertainty and no
final SCA promotion.
