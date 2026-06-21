---
status: superseded_handoff
decision_use: historical_audit_not_current_context
created_at: 2026-06-21 20:49 Asia/Shanghai
active_reference_id: h_r0_shadow_ref_20260619
branch_a_overlay: merge_gap_bp_2000000
final_impact: development_review_only
---

# Branch B V2 Report Burden Pass2 Handoff

> Superseded on 2026-06-21 by
> `docs/handoff/2026-06-21_2123_branch_b_v2_pass2_correction_handoff.md`.
> The sample-level `report_event_count >= 3` condition was removed from
> candidate-level demotion logic. The materialized pass2 result below is kept
> only as historical output-burden exploration.

## Required Context

Read these first in the next session:

1. `docs/CURRENT_CONTEXT_INDEX.md`
2. `docs/handoff/2026-06-21_2049_branch_b_v2_report_burden_pass2_handoff.md`
3. `AGENTS.md`
4. `skills/conversation_handoff/SKILL.md`
5. `skills/pgta_reference_modeling_analysis/SKILL.md`
6. Current result files and configs for the active task

## Completed This Loop

- Added a second Branch B V2 report-layer burden pass.
- The pass demotes selected `report_event` rows to `internal_review_event`.
- It does not hard-filter truth-overlap candidates and does not use legacy
  Branch B decision fields.
- Added benchmark outputs:
  - `report_events.tsv`
  - `report_events.json`
- Reran remote unit tests, dry-runs, and materialization for Y, H, G, and
  2026-06-15.
- Added the report-event audit:
  `docs/reports/branch_b_v2_report_event_audit_2026-06-21.md`.

## Changed Files

Core workflow/code:

- `pgta/predict/branch_b/v2_classifier.py`
- `pgta/predict/branch_b/v2_benchmark.py`
- `rules/predict_layout.smk`
- `rules/predict_workflow.smk`
- `rules/target_assembly.smk`

Tests:

- `tests/unit/test_branch_b_v2_classifier.py`
- `tests/unit/test_branch_b_v2_benchmark.py`
- `tests/unit/test_branch_ab_phase12_workflow_contract.py`
- `tests/unit/test_current_context_index.py`

Docs:

- `docs/CURRENT_CONTEXT_INDEX.md`
- `CURRENT_STATE.md`
- `PLANS.md`
- `docs/reports/branch_b_v2_report_event_audit_2026-06-21.md`
- `docs/handoff/2026-06-21_2049_branch_b_v2_report_burden_pass2_handoff.md`

## Pass2 Rule

The new demotion rule applies only when combined evidence is present:

```text
sample report_event count >= 3
+ unknown/no-null background context
+ unstable GC/RC context
+ a_abs_zscore < 50
```

Action:

```text
report_event -> internal_review_event
```

Reason:

```text
multi_report_unknown_background_gc_rc_unstable_internal_review
```

No single weak indicator is sufficient to demote or filter. The rule does not
create B-only events.

## Materialized Result

| cohort | candidates | truth | preserved | FN | truth filtered | report before | report after | internal review | filtered | Branch S |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Y1-Y8 | 97 | 10 | 10 | 0 | 0 | 21 | 9 | 62 | 13 | 13 |
| H1-H16 | 105 | 10 | 10 | 0 | 0 | 6 | 4 | 51 | 8 | 42 |
| G1-G8 | 75 | 10 | 10 | 0 | 0 | 15 | 12 | 43 | 7 | 13 |
| 2026-06-15 | 165 | 0 | 0 | 0 | 0 | 52 | 15 | 113 | 23 | 14 |

H6 chr21 remains preserved and visible as `internal_review_event` with
`top_a_abs_zscore=7.1135`.

G2 remains the G-batch burden outlier, now reduced from 6 report events to 3
without filtering its truth events.

## Remote Commands

Unit tests:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 /biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest tests/unit/test_branch_b_v2_classifier.py tests/unit/test_branch_b_v2_benchmark.py tests/unit/test_cnv_report.py tests/unit/test_branch_ab_phase12_workflow_contract.py tests/unit/test_current_context_index.py -q'
```

Result:

```text
71 passed in 1.10s
```

Dry-run target for each active gap2m config:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile --configfile <config> --cores 1 -n branch_b_v2_benchmark branch_s_review cnv_report
```

Result:

```text
RC=0 for Y, H, G, and 2026-06-15 configs.
```

Forced materialization target for each active gap2m config:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile --configfile <config> --cores 16 --forcerun cnv_branch_b_v2_classifier cnv_branch_b_v2_benchmark cnv_report_summary branch_b_v2_benchmark branch_s_review cnv_report
```

Result:

```text
RC=0 for Y, H, G, and 2026-06-15 configs.
```

## Current Boundaries

- Reference remains `fixed_shadow_baseline_not_production`.
- `merge_gap_bp=2_000_000` remains an explicit overlay, not default.
- Branch B V2 remains development/report-layer logic, not production-final.
- Branch S remains `review_reportable_with_limitations`.
- 2026-06-15 remains context only; do not calculate TP/FN/FP for it.
- CNVpro/CNVseq concepts are evidence tags and tiering inspiration only; they
  do not replace WisecondorX/CBS, `CNVcalling.R`, or cghFLasso.

## Next Step

Before adding any further demotion or filtering rule:

1. Inspect remaining `report_events.tsv` rows, especially G2 and 0615 sample
   `JZ26125846-61-61`.
2. Confirm whether remaining report events are strong A signals, truth-like
   weak positives, region-risk artifacts, or unresolved background cases.
3. Run any proposed rule against Y1-Y8, H1-H6, and G1-G8 first.
4. Keep H6 chr21 and all G truth events protected.
