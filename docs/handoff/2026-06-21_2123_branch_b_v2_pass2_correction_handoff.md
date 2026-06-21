---
status: active_handoff
decision_use: current_context
created_at: 2026-06-21 21:23 Asia/Shanghai
active_reference_id: h_r0_shadow_ref_20260619
branch_a_overlay: merge_gap_bp_2000000
final_impact: development_review_only
---

# Branch B V2 Pass2 Correction Handoff

## Required Context

Read these first in the next session:

1. `docs/CURRENT_CONTEXT_INDEX.md`
2. `docs/handoff/2026-06-21_2123_branch_b_v2_pass2_correction_handoff.md`
3. `AGENTS.md`
4. `skills/conversation_handoff/SKILL.md`
5. `skills/pgta_reference_modeling_analysis/SKILL.md`
6. Current result files and configs for the active task

## Completed This Loop

- Retracted the previous Branch B V2 pass2 candidate-demotion rule that used
  same-sample report-event burden.
- Removed sample report-event count from candidate-level demotion/filtering.
- Kept report-event burden as sample-level audit output only:
  `sample_report_burden_flag` and `sample_report_burden_reason`.
- Marked the previous pass2 report and handoff as superseded methodology audit.
- Added the current correction report:
  `docs/reports/branch_b_v2_pass2_correction_2026-06-21.md`.

## Current Method Boundary

`sample report_event count >= 3` is not candidate evidence.

It may appear in benchmark sample summary and aggregate summary as:

```text
sample_report_burden_flag
sample_report_burden_reason
sample_report_burden_threshold
sample_report_burden_flag_count
```

It must not change:

```text
v2_report_layer_class
v2_report_visibility
v2_filter_action
v2_burden_reduction_action
```

## Preserved Constraints

- Do not modify Branch A.
- Do not rebuild reference.
- Do not promote `merge_gap_bp=2_000_000` to default.
- Do not promote Branch B V2 or Branch S to production-final.
- Do not use legacy/current-code Branch B fields as V2 decision evidence.
- Do not use 2026-06-15 as truth.

## Remote Validation Status

Remote validation and materialization completed.

```text
tests/unit/test_branch_b_v2_classifier.py
tests/unit/test_branch_b_v2_benchmark.py
tests/unit/test_cnv_report.py
tests/unit/test_branch_ab_phase12_workflow_contract.py
tests/unit/test_current_context_index.py
73 passed in 1.24s
```

Dry-run target:

```text
branch_b_v2_benchmark branch_s_review cnv_report
```

Result: RC=0 for Y, H, G, and 2026-06-15 active gap2m configs.

Forced materialization completed for Branch B V2 classifier, benchmark, Branch
S, and report outputs. Report summaries were refreshed after adding sample
burden display fields.

## Materialized Counts

| cohort | candidates | truth | preserved | FN | truth filtered | report | internal review | filtered | Branch S | burden flags |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Y1-Y8 | 97 | 10 | 10 | 0 | 0 | 21 | 50 | 13 | 13 | 2 |
| H1-H16 | 105 | 10 | 10 | 0 | 0 | 6 | 49 | 8 | 42 | 1 |
| G1-G8 | 75 | 10 | 10 | 0 | 0 | 15 | 40 | 7 | 13 | 1 |
| 2026-06-15 | 165 | 0 | 0 | 0 | 0 | 52 | 76 | 23 | 14 | 5 |

H6 chr21 remains visible as `internal_review_event` with
`top_a_abs_zscore=7.1135`.

G2 locked truth remains visible:

- chr8 gain: `internal_review_event`;
- chr11 loss: `report_event`.

## Next Gate

The next burden-reduction design must begin with candidate-level evidence
audit, not sample-level report count.
