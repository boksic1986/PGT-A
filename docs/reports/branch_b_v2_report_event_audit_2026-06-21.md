---
status: superseded_methodology_audit_only
decision_use: historical_audit_not_current_branch_b_v2_rule
date: 2026-06-21
active_reference_id: h_r0_shadow_ref_20260619
branch_a_overlay: merge_gap_bp_2000000
final_impact: development_review_only
---

# Branch B V2 Report Event Audit

> Superseded on 2026-06-21 by
> `docs/reports/branch_b_v2_pass2_correction_2026-06-21.md`.
> The `same sample has >=3 report events` condition is retained only as a
> sample-level burden audit flag. It must not drive per-candidate demotion,
> filtering, or report-layer class changes.

## Scope

This audit is the second Branch B V2 report-layer burden reduction pass. It does
not modify Branch A, does not rebuild the reference, and does not promote the
active shadow reference or Branch B V2 to production-final status.

The pass runs on the fixed contract:

- reference: `h_r0_shadow_ref_20260619`;
- Branch A input: WisecondorX/CBS candidates with explicit
  `merge_gap_bp=2_000_000` overlay;
- truth cohorts: Y1-Y8, H1-H6, and G1-G8;
- context-only cohorts: H7-H16 and 2026-06-15.

## Pass2 Rule

The new pass only demotes `report_event` to `internal_review_event`. It does not
hard-filter or delete events unless a workflow/reference contract risk exists.

The demotion requires combined evidence:

```text
same sample has >=3 report events
+ background/no-null context is unknown or missing null support
+ GC/RC context is unstable
+ Branch A absolute z-score is <50
```

The emitted reason is:

```text
multi_report_unknown_background_gc_rc_unstable_internal_review
```

Single indicators are still insufficient for demotion or filtering. Length,
background unknown, B-side signal context, GC/RC attenuation/amplification,
blacklist/soft-mask burden, or `a_abs_zscore < 10` alone must not suppress a
Branch A candidate.

## Burden Change

| cohort | candidates | truth | truth preserved | FN | truth filtered | report before | report after | internal review | filtered audit-only | Branch S |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Y1-Y8 | 97 | 10 | 10 | 0 | 0 | 21 | 9 | 62 | 13 | 13 |
| H1-H16 | 105 | 10 | 10 | 0 | 0 | 6 | 4 | 51 | 8 | 42 |
| G1-G8 | 75 | 10 | 10 | 0 | 0 | 15 | 12 | 43 | 7 | 13 |
| 2026-06-15 | 165 | 0 | 0 | 0 | 0 | 52 | 15 | 113 | 23 | 14 |

Interpretation:

- The main report-event burden dropped in all materialized cohorts.
- The pass is conservative for truth cohorts: Y/H/G truth remains 10/10
  preserved with FN=0.
- 2026-06-15 has no locked truth; it is burden/context only and must not be
  used to claim FP reduction or tune truth thresholds.
- Demoted events remain visible in internal review output and are not erased
  from audit evidence.

## Sample-Level Counts

### Y1-Y8

| sample | candidates | truth | preserved | FN | report | internal review | filtered | Branch S |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Y1 | 11 | 1 | 1 | 0 | 2 | 6 | 0 | 3 |
| Y2 | 34 | 2 | 2 | 0 | 1 | 24 | 8 | 1 |
| Y3 | 3 | 1 | 1 | 0 | 0 | 1 | 0 | 2 |
| Y4 | 15 | 1 | 1 | 0 | 1 | 12 | 1 | 1 |
| Y5 | 10 | 1 | 1 | 0 | 1 | 7 | 2 | 0 |
| Y6 | 8 | 2 | 2 | 0 | 2 | 2 | 1 | 3 |
| Y7 | 3 | 1 | 1 | 0 | 1 | 2 | 0 | 0 |
| Y8 | 13 | 1 | 1 | 0 | 1 | 8 | 1 | 3 |

### H1-H16

| sample | candidates | truth | preserved | FN | report | internal review | filtered | Branch S |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| H1 | 17 | 1 | 1 | 0 | 1 | 10 | 3 | 3 |
| H2 | 16 | 1 | 1 | 0 | 1 | 12 | 0 | 3 |
| H3 | 3 | 1 | 1 | 0 | 0 | 3 | 0 | 0 |
| H4 | 7 | 2 | 2 | 0 | 1 | 2 | 3 | 1 |
| H5 | 2 | 3 | 3 | 0 | 0 | 0 | 0 | 2 |
| H6 | 3 | 2 | 2 | 0 | 0 | 1 | 0 | 2 |
| H7 | 5 | 0 | 0 | 0 | 0 | 1 | 1 | 3 |
| H8 | 6 | 0 | 0 | 0 | 0 | 2 | 1 | 3 |
| H9 | 3 | 0 | 0 | 0 | 0 | 0 | 0 | 3 |
| H10 | 3 | 0 | 0 | 0 | 0 | 0 | 0 | 3 |
| H11 | 3 | 0 | 0 | 0 | 0 | 0 | 0 | 3 |
| H12 | 4 | 0 | 0 | 0 | 0 | 0 | 0 | 4 |
| H13 | 14 | 0 | 0 | 0 | 1 | 10 | 0 | 3 |
| H14 | 13 | 0 | 0 | 0 | 0 | 10 | 0 | 3 |
| H15 | 3 | 0 | 0 | 0 | 0 | 0 | 0 | 3 |
| H16 | 3 | 0 | 0 | 0 | 0 | 0 | 0 | 3 |

H6 chr21 gain remains preserved and visible as `internal_review_event`, with
`top_a_abs_zscore=7.1135`, `GC_RC_ATTENUATED`, and
`B_SIGNAL_SUPPORTED_A_DIRECTION`.

### G1-G8

| sample | candidates | truth | preserved | FN | report | internal review | filtered | Branch S |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| G1 | 4 | 1 | 1 | 0 | 1 | 3 | 0 | 0 |
| G2 | 28 | 2 | 2 | 0 | 3 | 19 | 3 | 3 |
| G3 | 13 | 2 | 2 | 0 | 1 | 10 | 0 | 2 |
| G4 | 4 | 1 | 1 | 0 | 2 | 2 | 0 | 0 |
| G5 | 8 | 1 | 1 | 0 | 1 | 3 | 2 | 2 |
| G6 | 8 | 1 | 1 | 0 | 2 | 3 | 0 | 3 |
| G7 | 4 | 1 | 1 | 0 | 1 | 1 | 2 | 0 |
| G8 | 6 | 1 | 1 | 0 | 1 | 2 | 0 | 3 |

G2 remains the main G-batch burden outlier. This pass reduced G2 report events
from 6 to 3 without filtering either G2 locked truth event. Further filtering
must first inspect G2 event-level evidence rather than add broad rules.

### 2026-06-15

| sample | candidates | report | internal review | filtered | Branch S |
|---|---:|---:|---:|---:|---:|
| JZ26125843-56-56 | 26 | 2 | 14 | 8 | 2 |
| JZ26125844-59-59 | 35 | 1 | 27 | 6 | 1 |
| JZ26125845-60-60 | 47 | 4 | 38 | 0 | 5 |
| JZ26125846-61-61 | 33 | 6 | 23 | 1 | 3 |
| JZ26125847-62-62 | 24 | 2 | 11 | 8 | 3 |

These five samples still have no locked truth table. Their counts are report
burden/context only.

## Output Files

New benchmark audit outputs are generated under each cohort's active gap2m
benchmark directory:

```text
wisecondorx/cnv/postprocess_gap2m/v2_benchmark/report_events.tsv
wisecondorx/cnv/postprocess_gap2m/v2_benchmark/report_events.json
```

Existing filtered-event audit files remain:

```text
wisecondorx/cnv/postprocess_gap2m/v2_benchmark/filtered_events.tsv
wisecondorx/cnv/postprocess_gap2m/v2_benchmark/filtered_events.json
```

## Remote Validation

Remote unit tests on `fengxian`:

```text
tests/unit/test_branch_b_v2_classifier.py
tests/unit/test_branch_b_v2_benchmark.py
tests/unit/test_cnv_report.py
tests/unit/test_branch_ab_phase12_workflow_contract.py
tests/unit/test_current_context_index.py
71 passed in 1.10s
```

Snakemake dry-runs passed for
`branch_b_v2_benchmark branch_s_review cnv_report` under all four active gap2m
configs:

```text
config_predict_y_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml
config_predict_h_20260608_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml
config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml
config_predict_g1_g8_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml
```

Forced materialization refreshed Branch B V2 classifier, benchmark, Branch S,
and report outputs for Y, H, G, and 2026-06-15.

## Interpretation

This pass is useful because it makes Branch B V2 perform a real report-layer
burden reduction while preserving locked truth. It is still not production
promotion:

- the reference remains a fixed shadow baseline;
- `merge_gap_bp=2Mb` remains an explicit overlay, not the default;
- Branch S remains `review_reportable_with_limitations`;
- 2026-06-15 remains context-only;
- further demotion/filter rules must be ablated on Y/H/G truth first.
