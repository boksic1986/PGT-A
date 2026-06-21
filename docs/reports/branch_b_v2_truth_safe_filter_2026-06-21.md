# Branch B V2 Truth-Safe Filter Contract

Date: 2026-06-21

Status: `development_shadow_contract`

## Purpose

This document records the next Branch B V2 contract after the method reset:
Branch B V2 is a Branch-A-anchored evidence and filter-action layer. It still
does not create B-only events, does not replace WisecondorX/CBS, and does not
promote any candidate into a final report by itself.

The goal of this step is to make filter behavior explicit and auditable while
keeping locked truth safe. The filter layer can downgrade or route candidates
for review, but it must not hard-suppress a Branch A positive unless the reason
is a workflow/reference contract risk.

## Inputs

Input remains one row per Branch A candidate:

- `sample_id`, `chrom`, `start`, `end`, `state`
- Branch A support fields such as `a_abs_zscore` and `a_support_level`
- Branch B V2 auxiliary evidence such as clean-bin support, GC/RC context,
  background/null context, B-side signal context, and refmap status

Legacy/current-code Branch B fields remain ignored for V2 decision evidence:

- `final_disposition`
- `branch_b_keep_event`
- `branch_b_report_class`
- `branch_b_artifact_status`

## New Output Fields

The classifier now emits these filter-specific fields:

| field | meaning |
|---|---|
| `v2_filter_version` | filter contract version |
| `v2_filter_action` | explicit candidate-level filter action |
| `v2_filter_reason` | machine-readable reason for the action |
| `v2_filter_scope` | scope of the action, such as review, Branch S, or workflow contract |
| `v2_filter_hard_suppression_allowed` | `1` only when hard suppression is allowed by contract |

The benchmark also records:

- `top_v2_filter_action` in truth-overlap metrics;
- `v2_filter_suppressed_count` per sample;
- `v2_filter_suppressed_candidate_count` in the benchmark summary.

## Action Semantics

Current action set:

| action | behavior |
|---|---|
| `keep_report_candidate` | keep candidate in report-tier review context |
| `keep_background_unknown_review` | keep candidate because background is unknown or not informative |
| `keep_review_candidate` | keep candidate as default review context |
| `keep_background_compatible_review` | keep candidate as background-compatible review context, not benign by itself |
| `downgrade_to_technical_risk_review` | downgrade to technical review, without hard suppression |
| `route_to_branch_s_review` | route sex-chromosome candidate to Branch S review |
| `suppress_workflow_contract_risk` | hard suppression allowed only for workflow/reference contract risk |

## Hard-Suppression Policy

Only workflow/reference contract risk can use a hard-suppression action in V2.
The current implemented case is:

```text
refmap_status = SAME_CHROM_REF_LEAKAGE
-> v2_filter_action = suppress_workflow_contract_risk
```

The following are not hard-suppression reasons in this contract:

- GC/RC attenuation;
- B-side signal discordance with Branch A direction;
- unknown background or no calibration null;
- length tier;
- low clean support or high region-risk burden.

Low clean support and high-risk-region burden can downgrade a candidate to
`technical_risk_review`, but this is a review disposition, not a biological
negative call.

## Benchmark Semantics

The benchmark treats `v2_filter_action=suppress_workflow_contract_risk` as hard
suppression. If this action overlaps locked truth, it counts as FN and as
`truth_hard_suppressed_count`.

This is intentional: the workflow must not hide truth loss behind a contract
risk label. Any future hard action must prove truth safety through Y1-Y8 and
H1-H6 locked truth before it can be considered.

## Current Validation

Remote unit validation on `ssh fengxian`:

```bash
cd /data/project/CNV/PGT-A/refactor_validation_20260419
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_v2_classifier.py \
  tests/unit/test_branch_b_v2_benchmark.py -q
```

Initial result:

```text
26 passed in 0.90s
```

Broader remote unit validation after synchronizing this contract:

```text
55 passed in 1.37s
```

Snakemake dry-runs succeeded for
`branch_b_v2_benchmark branch_s_review cnv_report` under all three active gap2m
configs:

- `config_predict_y_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`
- `config_predict_h_20260608_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`
- `config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`

Forced materialization then refreshed:

```text
cnv_branch_b_v2_classifier
cnv_branch_b_v2_benchmark
branch_b_v2_benchmark
```

Remote materialized benchmark result:

| cohort | candidates | truth events | truth preserved | FN | hard-suppressed truth | filter-suppressed candidates | status |
|---|---:|---:|---:|---:|---:|---:|---|
| Y1-Y8 | 97 | 10 | 10 | 0 | 0 | 0 | ready |
| H1-H16 | 105 | 10 | 10 | 0 | 0 | 0 | ready |
| 2026-06-15 | 165 | 0 | 0 | 0 | 0 | 0 | skipped_no_truth |

Remote result-file check confirmed:

- classifier TSV contains `v2_filter_action`;
- truth metrics TSV contains `top_v2_filter_action`;
- sample summary TSV contains `v2_filter_suppressed_count`.

This is now a materialized contract validation. It still is not a final
performance-promotion claim because no FP/review-burden reduction has been
demonstrated.

## Current Conclusion

The Branch B V2 filter contract now separates:

- evidence labels;
- review/report disposition;
- explicit filter action;
- hard-suppression permission.

This is the correct direction for reducing FP/review burden later, but this
step itself does not claim FP reduction, Branch S finalization, report release,
or production reference promotion.
