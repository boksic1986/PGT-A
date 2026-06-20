# Branch B P3 Evidence Contract

Date: 2026-06-20

## Scope

This report records the P3 Branch B candidate-level evidence refinement
materialized under:

```text
reference_id: h_r0_shadow_ref_20260619
reference status: shadow-only
```

P3 refines Branch A candidates only. It does not create B-only final events and
does not promote, suppress, or release report events.

## Workflow Contract

New P3 target:

```text
branch_b_evidence
```

The target depends on existing current-code Branch B evidence inputs:

```text
Branch A candidates
Branch B calibrated bins
Branch B current-code final/artifact-rule events
run metadata
```

The target outputs only:

```text
wisecondorx/cnv/postprocess/evidence_ledger/{sample}.candidate_evidence.tsv
wisecondorx/cnv/postprocess/evidence_ledger/{sample}.summary.json
```

The target does not request:

```text
Branch B V2 classifier
matched-negative percentile
Branch S
cnv_report
```

## Schema Contract

Each Branch A candidate has one Branch B evidence row. The P3-specific fields
are:

```text
branch_b_direction_support
copy_number_like_amplitude
mosaic_proxy
loh_evidence
upd_evidence
background_source
background_status
region_risk_context
sample_noise_context
cnvpro_consistency_status
disposition
disposition_reason
report_impact
```

Current missing or unavailable evidence remains explicit:

```text
background_source=UNKNOWN_BACKGROUND
background_status=UNKNOWN_BACKGROUND
loh_evidence=not_available
upd_evidence=not_available
report_impact=none_shadow_only
```

Legacy Branch B final dispositions are retained as context in
`final_disposition` and `disposition_reason`, but do not become P3 hard
confirmation or suppression.

## Remote Outputs

Y1-Y8:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419/results_y_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess/evidence_ledger/
```

H1-H16:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419/results_h_20260608_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess/evidence_ledger/
```

2026-06-15 exploratory samples:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419/results_20260615_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess/evidence_ledger/
```

## Summary

| cohort | ledger files | summary files | Branch A candidate rows | P3 disposition |
|---|---:|---:|---:|---|
| Y1-Y8 | 8 | 8 | 131 | REVIEW_REQUIRED |
| H1-H16 | 16 | 16 | 221 | REVIEW_REQUIRED |
| 2026-06-15 | 5 | 5 | 201 | REVIEW_REQUIRED |

## Remote Validation

Remote unit tests:

```text
35 passed in 0.89s
```

Remote dry-runs passed for:

```text
config_predict_y_h_r0_shadow_branch_b_evidence_20260620.yaml
config_predict_h_20260608_h_r0_shadow_branch_b_evidence_20260620.yaml
config_predict_20260615_h_r0_shadow_branch_b_evidence_20260620.yaml
```

Remote forced materialization completed for the same three configs by rerunning
only `cnv_branch_b_evidence_ledger` under the `branch_b_evidence` target.

Post-materialization schema check:

```text
Y: ledger_files=8, summary_files=8, rows=131, missing_schema=none
H: ledger_files=16, summary_files=16, rows=221, missing_schema=none
0615: ledger_files=5, summary_files=5, rows=201, missing_schema=none
```

## Interpretation

P3 currently satisfies the engineering contract that every fixed Branch A
candidate has a Branch B evidence row. The output is review-safe: missing
background remains `UNKNOWN_BACKGROUND`, unavailable LOH/UPD remains
`not_available`, and no Branch B row is promoted into a final report call.

This does not satisfy P5/SCA or P6 report release. The next active gate is P5
Branch S/SCA report-boundary validation.
