---
status: validation_complete
decision_use: current_branch_b_v2_method_correction
date: 2026-06-21
active_reference_id: h_r0_shadow_ref_20260619
branch_a_overlay: merge_gap_bp_2000000
final_impact: development_review_only
---

# Branch B V2 Pass2 Correction

## Scope

This correction retracts the previous Branch B V2 report-burden pass2
candidate demotion rule. The retracted rule used same-sample report-event
burden as part of the evidence for changing a single candidate from
`report_event` to `internal_review_event`.

The corrected boundary is:

```text
sample report_event count >= 3
```

is a sample-level burden audit signal only. It must not change a candidate's
`v2_report_layer_class`, `v2_report_visibility`, `v2_filter_action`, or
`v2_burden_reduction_action`.

## Superseded Result

The previous document
`docs/reports/branch_b_v2_report_event_audit_2026-06-21.md` is retained only as
`superseded_methodology_audit_only`.

It can be used to understand output burden exploration, but it cannot be used
as the current Branch B V2 decision rule or performance evidence.

## Corrected Contract

- Branch A remains unchanged: WisecondorX/CBS candidates with explicit
  `merge_gap_bp=2_000_000` overlay for this R&D branch.
- Branch B V2 remains Branch-A-anchored and cannot create B-only events.
- Legacy/current-code Branch B fields remain excluded from V2 decisions:
  `final_disposition`, `branch_b_keep_event`, and legacy artifact labels.
- `report_events.tsv/json` remain valid benchmark/audit outputs.
- `sample_report_burden_flag` and `sample_report_burden_reason` are emitted
  only in benchmark sample summaries and aggregate benchmark JSON.
- Future demotion/filtering must use candidate-level combination evidence, not
  truth labels, sample report counts, or 2026-06-15 burden counts.

## Validation Plan

Remote validation must rerun on `ssh fengxian` in:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419
```

Required unit tests:

```text
tests/unit/test_branch_b_v2_classifier.py
tests/unit/test_branch_b_v2_benchmark.py
tests/unit/test_cnv_report.py
tests/unit/test_branch_ab_phase12_workflow_contract.py
tests/unit/test_current_context_index.py
```

Required dry-runs:

```text
branch_b_v2_benchmark branch_s_review cnv_report
```

for the Y/H/G/2026-06-15 active gap2m configs.

Required materialized checks:

| cohort | required gate |
|---|---|
| Y1-Y8 | truth preserved 10/10, FN=0, truth filtered=0 |
| H1-H16 | H1-H6 truth preserved 10/10, FN=0, H6 chr21 visible |
| G1-G8 | truth preserved 10/10, FN=0, G2 truth not filtered |
| 2026-06-15 | burden/context only, no TP/FN/FP claims |

## Materialized Result

Remote full validation and materialization completed on `fengxian`.

| cohort | candidates | truth | truth preserved | FN | truth filtered | report | internal review | filtered audit-only | Branch S | sample burden flags |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Y1-Y8 | 97 | 10 | 10 | 0 | 0 | 21 | 50 | 13 | 13 | 2 |
| H1-H16 | 105 | 10 | 10 | 0 | 0 | 6 | 49 | 8 | 42 | 1 |
| G1-G8 | 75 | 10 | 10 | 0 | 0 | 15 | 40 | 7 | 13 | 1 |
| 2026-06-15 | 165 | 0 | 0 | 0 | 0 | 52 | 76 | 23 | 14 | 5 |

Interpretation:

- The corrected rule restores the previous valid report-layer filter counts.
- `sample_report_burden_flag` is now audit-only and does not demote candidate
  events.
- Y/H/G locked truth remains preserved with FN=0 and no report-layer filtered
  truth.
- H6 chr21 remains visible as `internal_review_event`.
- G2 locked truth remains visible: chr8 gain is `internal_review_event`; chr11
  loss is `report_event`.
- 2026-06-15 remains burden/context only and must not be used for TP/FN/FP.

## Remote Validation

```text
tests/unit/test_branch_b_v2_classifier.py
tests/unit/test_branch_b_v2_benchmark.py
tests/unit/test_cnv_report.py
tests/unit/test_branch_ab_phase12_workflow_contract.py
tests/unit/test_current_context_index.py
73 passed in 1.24s
```

Snakemake dry-runs passed for `branch_b_v2_benchmark branch_s_review
cnv_report` under all four active gap2m configs:

```text
config_predict_y_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml
config_predict_h_20260608_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml
config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml
config_predict_g1_g8_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml
```

Forced materialization completed for Branch B V2 classifier, benchmark, Branch
S, and report outputs. After adding report-level sample burden display, report
summary outputs were refreshed for all four configs.

`cnv_summary.json` now carries:

```text
branch_b_v2_sample_report_burden_flag_count
branch_b_v2_sample_report_burden_threshold
branch_b_v2_final_impact=development_review_only
branch_b_v2_legacy_fields_used=false
```

## Next Rule Gate

The next Branch B V2 burden reduction loop must start from candidate-level
evidence audit:

- A z strength;
- length tier;
- clean support;
- high-risk region / acrocentric / telomere / centromere context;
- B-side support;
- GC/RC context;
- mosaic-sensitive weak-positive protection.

Single indicators remain insufficient for demotion or filtering.
