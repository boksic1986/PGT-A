# Branch A Validation Under H R0 Shadow Reference

Date: 2026-06-20

## Scope

This report records the P2 Branch A no-FN validation materialized under:

```text
reference_id: h_r0_shadow_ref_20260619
reference status: shadow-only
```

It validates Branch A sensitivity and candidate burden only. It does not promote
Branch B, Branch S/SCA, CNVpro-inspired evidence, or 2026-06-15 reports.

## Workflow Contract

New P2 target:

```text
branch_a_validation
```

The target depends on:

```text
WisecondorX predict output -> Branch A candidates
CNV QC TSV
gender TSV when sex-specific prediction is enabled
truth TSV when configured
```

It does not consume:

```text
Branch B final events
Branch B artifact summaries
matched-negative outputs
Branch S outputs
cnv_report outputs
```

## Remote Outputs

Y1-Y8:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419/results_y_h_r0_shadow_ref_20260619/wisecondorx/cnv/branch_a_validation/
```

H1-H16:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419/results_h_20260608_h_r0_shadow_ref_20260619/wisecondorx/cnv/branch_a_validation/
```

2026-06-15 exploratory samples:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419/results_20260615_h_r0_shadow_ref_20260619/wisecondorx/cnv/branch_a_validation/
```

Each output directory contains:

```text
sample_summary.tsv
truth_metrics.tsv
summary.json
```

## Summary

| cohort | sample_count | truth_events | truth_detected | FN | Branch A candidates | key note |
|---|---:|---:|---:|---:|---:|---|
| Y1-Y8 | 8 | 10 | 10 | 0 | 131 | known positives detected |
| H1-H16 | 16 | 10 | 10 | 0 | 221 | H6 chr21 gain detected |
| 2026-06-15 | 5 | 0 | 0 | 0 | 201 | exploratory burden only; no truth table |

## Per-Sample Candidate Burden

### Y1-Y8

| sample | truth | detected | FN | A candidates | A strong | A sensitive | top abs z |
|---|---:|---:|---:|---:|---:|---:|---:|
| Y1 | 1 | 1 | 0 | 20 | 14 | 6 | 132.59 |
| Y2 | 2 | 2 | 0 | 37 | 13 | 24 | 94.78 |
| Y3 | 1 | 1 | 0 | 6 | 5 | 1 | 59.59 |
| Y4 | 1 | 1 | 0 | 15 | 3 | 12 | 123.90 |
| Y5 | 1 | 1 | 0 | 10 | 3 | 7 | 79.32 |
| Y6 | 2 | 2 | 0 | 18 | 16 | 2 | 131.34 |
| Y7 | 1 | 1 | 0 | 3 | 2 | 1 | 40.90 |
| Y8 | 1 | 1 | 0 | 22 | 16 | 6 | 156.22 |

### H1-H16

| sample | truth | detected | FN | A candidates | A strong | A sensitive | top abs z | H6 chr21 |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| H1 | 1 | 1 | 0 | 26 | 15 | 11 | 132.80 | not_applicable |
| H2 | 1 | 1 | 0 | 26 | 17 | 9 | 130.31 | not_applicable |
| H3 | 1 | 1 | 0 | 3 | 2 | 1 | 36.76 | not_applicable |
| H4 | 2 | 2 | 0 | 8 | 3 | 5 | 111.55 | not_applicable |
| H5 | 3 | 3 | 0 | 5 | 5 | 0 | 62.16 | not_applicable |
| H6 | 2 | 2 | 0 | 6 | 5 | 1 | 62.53 | detected |
| H7 | 0 | 0 | 0 | 14 | 12 | 2 | 131.80 | not_applicable |
| H8 | 0 | 0 | 0 | 15 | 12 | 3 | 134.63 | not_applicable |
| H9 | 0 | 0 | 0 | 12 | 12 | 0 | 133.10 | not_applicable |
| H10 | 0 | 0 | 0 | 12 | 12 | 0 | 135.09 | not_applicable |
| H11 | 0 | 0 | 0 | 12 | 12 | 0 | 132.69 | not_applicable |
| H12 | 0 | 0 | 0 | 13 | 13 | 0 | 136.15 | not_applicable |
| H13 | 0 | 0 | 0 | 23 | 14 | 9 | 124.50 | not_applicable |
| H14 | 0 | 0 | 0 | 22 | 12 | 10 | 148.03 | not_applicable |
| H15 | 0 | 0 | 0 | 12 | 12 | 0 | 131.31 | not_applicable |
| H16 | 0 | 0 | 0 | 12 | 12 | 0 | 133.42 | not_applicable |

### 2026-06-15 Exploratory Samples

| sample | A candidates | A strong | A sensitive | top abs z |
|---|---:|---:|---:|---:|
| JZ26125843-56-56 | 28 | 4 | 24 | 13.90 |
| JZ26125844-59-59 | 38 | 20 | 18 | 30.11 |
| JZ26125845-60-60 | 57 | 49 | 8 | 122.76 |
| JZ26125846-61-61 | 44 | 30 | 14 | 121.88 |
| JZ26125847-62-62 | 34 | 17 | 17 | 125.74 |

## Interpretation

P2 currently supports the statement that Branch A recall does not regress on the
available Y1-Y8 and H1-H6 truth sets under the `h_r0_shadow_ref_20260619`
shadow reference.

The result does not support report release. Candidate burden remains high in
H7-H16 and especially in the 2026-06-15 exploratory samples. Those candidates
must be handled by later Branch B evidence refinement and Branch S/SCA gates,
without treating Branch B legacy kept counts or N1 shadow background as final
filters.

## Remote Validation

Remote tests:

```text
16 passed in 0.73s
```

Remote Snakemake dry-runs passed for:

```text
config_predict_y_h_r0_shadow_branch_a_validation_20260620.yaml
config_predict_h_20260608_h_r0_shadow_branch_a_validation_20260620.yaml
config_predict_20260615_h_r0_shadow_branch_a_validation_20260620.yaml
```

Remote materialization completed for the same three configs.

## Next Gate

Proceed to P3/P5 planning:

```text
P3: Branch B candidate-level evidence refinement
P5: Branch S/SCA report-boundary validation
```

Do not generate a formal 2026-06-15 report until G3 and G5 pass under the same
reference/config contract.
