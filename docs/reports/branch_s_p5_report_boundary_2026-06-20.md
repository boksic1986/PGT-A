# Branch S P5 Report Boundary

Date: 2026-06-20

## Scope

This report records the P5 Branch S/SCA report-boundary materialization under:

```text
reference_id: h_r0_shadow_ref_20260619
reference status: shadow-only
```

P5 defines how SCA evidence can be represented in workflow outputs. It does not
promote SCA evidence to final reportable calls.

## Workflow Contract

New P5 target:

```text
branch_s_review
```

The target outputs only:

```text
wisecondorx/cnv/postprocess/branch_s/{sample}.sex_chrom_evidence.tsv
wisecondorx/cnv/postprocess/branch_s/{sample}.sca_state_scores.tsv
wisecondorx/cnv/postprocess/branch_s/{sample}.summary.json
```

The target does not request:

```text
Branch B V2 classifier
matched-negative percentile
cnv_report
```

## Summary Schema

Each summary JSON now includes the P5 report-boundary fields:

```text
sample
sex_call
expected_x_ploidy
expected_y_ploidy
x_nonpar_direction
x_par_context
y_nonpar_direction
y_par_or_homology_context
branch_a_x_support
branch_a_y_support
sca_candidate_state
sca_confidence_tier
sca_output_mode
sca_uncertainty_reason
report_text_status
```

Current output mode remains:

```text
sca_output_mode=review_development_only
report_text_status=development_only_not_final_reportable
```

## Remote Outputs

Y1-Y8:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419/results_y_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess/branch_s/
```

H1-H16:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419/results_h_20260608_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess/branch_s/
```

2026-06-15 exploratory samples:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419/results_20260615_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess/branch_s/
```

## Materialized Output Summary

| cohort | evidence tables | score tables | summary JSON | output mode |
|---|---:|---:|---:|---|
| Y1-Y8 | 8 | 8 | 8 | review_development_only |
| H1-H16 | 16 | 16 | 16 | review_development_only |
| 2026-06-15 | 5 | 5 | 5 | review_development_only |

SCA confidence tier counts:

| cohort | SCA_REVIEW_STRONG | SCA_REVIEW_WEAK | SCA_NO_CALL |
|---|---:|---:|---:|
| Y1-Y8 | 4 | 2 | 2 |
| H1-H16 | 14 | 1 | 1 |
| 2026-06-15 | 3 | 2 | 0 |

## Remote Validation

Remote unit tests:

```text
40 passed in 0.85s
```

Remote dry-runs passed for:

```text
config_predict_y_h_r0_shadow_branch_s_review_20260620.yaml
config_predict_h_20260608_h_r0_shadow_branch_s_review_20260620.yaml
config_predict_20260615_h_r0_shadow_branch_s_review_20260620.yaml
```

Remote forced materialization completed for the same three configs by rerunning
only `cnv_branch_s_shadow` under the `branch_s_review` target.

Post-materialization schema check:

```text
Y: evidence=8, scores=8, summaries=8, missing_schema=none
H: evidence=16, scores=16, summaries=16, missing_schema=none
0615: evidence=5, scores=5, summaries=5, missing_schema=none
```

## Interpretation

P5-review boundary is represented in workflow outputs. Every summary explicitly
states expected sex ploidy, X/Y non-PAR direction, PAR context, Branch A support,
candidate SCA state, review confidence tier, output mode, and uncertainty.

P5-final is not passed. Locked SCA truth coverage is still incomplete for X gain,
XXY/XXX/XYY-like samples, Y loss, mosaic SCA fraction series, clean XX/XY
negatives across batches, and PAR/XY-homology edge cases.

This means Branch S can be carried into a report package only as visible
review/development evidence, not as final SCA confirmation or suppression.
