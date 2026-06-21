---
status: active_handoff
decision_use: current_context
created_at: 2026-06-22 09:30 Asia/Shanghai
active_reference_id: h_r0_shadow_ref_20260619
branch_a_overlay: merge_gap_bp_2000000
final_impact: development_review_only
lowres_predict_materialization: complete
---

# Lowres Branch B/S Integration Handoff

## Context

Use this handoff after `docs/CURRENT_CONTEXT_INDEX.md`.

This loop connects the already completed 2Mb/3Mb shadow reference builds to
Branch B and Branch S as auxiliary evidence only.

The active 1Mb chain remains:

```text
h_r0_shadow_ref_20260619
WisecondorX predict/CBS Branch A
explicit merge_gap_bp=2_000_000 overlay
Branch B V2 report-layer filtering
Branch S review/development-only
```

## Completed In This Loop

Code:

- `rules/predict_layout.smk`
  - supports `lowres_evidence.reference_npz_glob`;
  - keeps explicit `lowres_evidence.reference_npz` support;
  - fails only when lowres is enabled and neither source is available.
- `rules/predict_workflow.smk`
  - passes optional lowres 2Mb/3Mb event files into Branch S when enabled.
- `pgta/predict/branch_s.py`
  - adds segment-level X/Y non-PAR evidence for local sex-chromosome CNV review;
  - keeps whole X/Y non-PAR median for global SCA trend;
  - treats PAR evidence as secondary context only;
  - emits lowres sex-chrom context fields with
    `sex_chrom_lowres_final_impact=context_only_not_filter`.

Tests:

- `tests/unit/test_branch_s_shadow.py`
  - local X/Y segment evidence can preserve sex-chrom CNV review even when
    global non-PAR median is neutral;
  - PAR-only signal remains context and does not create SCA call;
  - lowres absence for short sex-chrom segments is context only.
- `tests/unit/test_branch_ab_phase12_workflow_contract.py`
  - lowres glob and Branch S lowres inputs are covered by static contract.

Configs added to the branch:

- 8 lowres predict configs for Y/H/G/2026-06-15 at 2Mb and 3Mb.
- 4 lowres-enabled main-chain configs for Y/H/G/2026-06-15.

Report:

- `docs/reports/branch_b_s_lowres_integration_2026-06-22.md`.

## Remote Validation

Remote unit tests on `fengxian` passed:

```text
68 passed in 1.53s
```

Command:

```text
cd /data/project/CNV/PGT-A/refactor_validation_20260419
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_s_shadow.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py \
  tests/unit/test_lowres_evidence.py \
  tests/unit/test_ref_stability.py \
  tests/unit/test_branch_b_v2_classifier.py -q
```

2Mb/3Mb predict dry-runs parsed successfully for Y/H/G/2026-06-15 and planned
existing-BAM WisecondorX convert/predict jobs. No BWA/mapping rules were in the
job stats.

## Long Task Status

Lowres predict materialization was launched as a background remote task.

```text
PID: 4690
log: /data/project/CNV/PGT-A/refactor_validation_20260419/logs/driver/lowres_predict_2mb_3mb_20260622.log
command: /data/project/CNV/PGT-A/refactor_validation_20260419/logs/driver/lowres_predict_2mb_3mb_20260622.command.sh
target: 2Mb/3Mb predict outputs for Y, H, G, and 2026-06-15
runtime history: monitor/runtime.db absent; no historical estimate available
```

Final observed progress:

- Y 2Mb/3Mb complete: 8/8 predict done, 8/8 A-branch events, 8/8 V2 classifier,
  8/8 Branch S summary at each binsize.
- H 2Mb/3Mb complete: 16/16 at each output level and binsize.
- G 2Mb/3Mb complete: 8/8 at each output level and binsize.
- 2026-06-15 2Mb/3Mb complete: 5/5 at each output level and binsize.

Lowres-enabled main-chain dry-runs and materialization also passed for all four
active configs.

## Materialized Acceptance

| cohort | samples | candidates | truth | preserved | filtered truth | report | internal review | filtered audit | Branch S | status |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| Y1-Y8 | 8 | 97 | 10 | 10 | 0 | 21 | 50 | 13 | 13 | ready |
| H1-H16 | 16 | 105 | 10 | 10 | 0 | 6 | 49 | 8 | 42 | ready |
| G1-G8 | 8 | 75 | 10 | 10 | 0 | 15 | 40 | 7 | 13 | ready |
| 2026-06-15 | 5 | 165 | 0 | 0 | 0 | 52 | 76 | 23 | 14 | skipped_no_truth |

H6 chr21 remains visible as `internal_review_event` with
`LOWRES_2MB_3MB_SAME_DIRECTION_SUPPORT` and `REF_STABILITY_STABLE`.

## Current Design Rules

- Lowres same-direction support can improve confidence/explanation.
- Lowres absence cannot independently demote or filter.
- High ref-MAD weakens negative lowres interpretation.
- Small or boundary-diluted sex-chromosome candidates must not be hidden by
  neutral whole-chromosome non-PAR median or absent 2Mb/3Mb support.
- Branch S remains review/development-only and does not replace sex calling.

## Next Actions

1. Inspect lowres/ref-MAD evidence by truth event and by current report event.
2. Use same-direction lowres support to improve confidence/explanation only.
3. Keep lowres absence out of standalone demotion/filtering.
4. Review high ref-MAD regions before interpreting no-support context.
5. Continue Branch S/SCA validation separately with a broader locked truth panel.
