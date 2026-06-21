# Handoff: Branch B V2 Background Context Review Label

Date: 2026-06-21 11:03 Asia/Shanghai

Status: `active_current_handoff`

## 1. Goal

Continue Branch B V2 evidence/disposition refinement by making the unresolved
background state explicit:

```text
UNKNOWN_BACKGROUND + NO_NULL_SUPPORT
```

This loop adds review-context labels only. It does not change final report
logic, hard suppression, Branch A, Branch S, or WisecondorX predict.

## 2. Context Restored

Read before execution:

- `C:\Users\11217\.codex\attachments\b29e8565-0516-429a-91ff-991e2ad43c59\goal-objective.md`
- `docs/CURRENT_CONTEXT_INDEX.md`
- `docs/handoff/2026-06-21_1042_branch_b_v2_direction_support_label_handoff.md`
- `AGENTS.md`
- `skills/conversation_handoff/SKILL.md`
- `skills/pgta_reference_modeling_analysis/SKILL.md`
- `REPO_MAP.md`
- `CURRENT_STATE.md`
- `PLANS.md`

Active reference:

```text
h_r0_shadow_ref_20260619
```

Active Branch A overlay:

```text
merge_gap_bp=2_000_000
strong_z=10.0
postprocess.output_dir=wisecondorx/cnv/postprocess_gap2m
```

## 3. Completed Work

Implemented new Branch B V2 classifier review fields:

```text
v2_background_context_label
v2_background_context_reason
```

The sample summary now includes:

```text
background_context_label_counts
```

Primary materialized labels:

- `UNKNOWN_BACKGROUND_NO_NULL_SUPPORT`
- `SHADOW_BACKGROUND_NO_NULL_SUPPORT`

These fields are context labels only. They intentionally do not alter:

- WisecondorX predict;
- Branch A candidate discovery;
- `v2_candidate_class`;
- `v2_classifier_action`;
- `v2_evidence_tier`;
- `v2_direction_support_label`;
- `v2_final_report_impact`;
- Branch S status.

New report:

- `docs/reports/branch_b_v2_background_context_label_2026-06-21.md`

Code/tests changed:

- `pgta/predict/branch_b/v2_classifier.py`
- `tests/unit/test_branch_b_v2_classifier.py`

Core file updates:

- `CURRENT_STATE.md`
- `PLANS.md`
- `docs/CURRENT_CONTEXT_INDEX.md`

## 4. TDD / Verification

TDD red check on remote:

```text
3 failed, 15 passed
```

The failures were expected missing-field failures for:

- `v2_background_context_label`;
- `background_context_label_counts`.

Remote unit test after implementation:

```bash
cd /data/project/CNV/PGT-A/refactor_validation_20260419
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_v2_classifier.py \
  tests/unit/test_branch_b_v2_benchmark.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py -q
```

Result:

```text
35 passed in 0.91s
```

Remote Snakemake dry-runs after materialization:

```text
Nothing to be done (all requested files are present and up to date).
```

Configs:

- `config_predict_y_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`
- `config_predict_h_20260608_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`
- `config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`

Forced materialization used:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake \
  -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile \
  --configfile <config> \
  --cores 8 \
  --forcerun cnv_branch_b_v2_classifier cnv_branch_b_v2_benchmark \
  branch_b_v2_benchmark
```

All three forced materializations completed.

## 5. Current Evidence

Background context counts:

| cohort | rows | UNKNOWN background no null | SHADOW background no null |
|---|---:|---:|---:|
| Y1-Y8 truth cohort | 97 | 86 | 11 |
| H1-H16 | 105 | 69 | 36 |
| 2026-06-15 context | 165 | 155 | 10 |

Truth preservation:

| cohort | truth events | truth preserved | FN | hard-suppressed truth |
|---|---:|---:|---:|---:|
| Y1-Y8 | 10 | 10 | 0 | 0 |
| H1-H16 | 10 | 10 | 0 | 0 |
| 2026-06-15 | 0 | 0 | 0 | 0 |

H6 chr21 remains:

```text
candidate_id=H6.A0003
v2_candidate_class=V2_POSITIVE_SUPPORT_REVIEW
v2_background_context_label=UNKNOWN_BACKGROUND_NO_NULL_SUPPORT
v2_direction_support_label=B_DIRECTION_SUPPORTED
v2_final_report_impact=none_shadow_only
a_abs_zscore=7.113507302991461
```

## 6. Current Conclusion

Branch B V2 now exposes missing matched-negative and calibration-null support as
explicit review context. This prevents unknown background from being
misinterpreted as benign or background-compatible.

This is not FP reduction and not Branch B promotion. The next loop should use
the explicit context labels to design truth-safe burden stratification without
hard suppression.

## 7. Suggested Next Step

Continue Branch B V2 refinement from this state:

1. Keep `UNKNOWN_BACKGROUND_NO_NULL_SUPPORT` as review context only.
2. Preserve Y1-Y8/H1-H6 truth 10/10, H6 chr21, and hard-suppressed truth=0.
3. Continue excluding legacy/current-code Branch B decision fields.
4. Design the next truth-safe burden stratification using the new context
   label plus existing direction/region/sample-noise evidence.

## 8. Key Files

- `pgta/predict/branch_b/v2_classifier.py`
- `tests/unit/test_branch_b_v2_classifier.py`
- `docs/reports/branch_b_v2_background_context_label_2026-06-21.md`
- `CURRENT_STATE.md`
- `PLANS.md`
- `docs/CURRENT_CONTEXT_INDEX.md`

Remote result roots:

- `/data/project/CNV/PGT-A/refactor_validation_20260419/results_y_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m`
- `/data/project/CNV/PGT-A/refactor_validation_20260419/results_h_20260608_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m`
- `/data/project/CNV/PGT-A/refactor_validation_20260419/results_20260615_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m`

## 9. Core File Sync

- `CURRENT_STATE.md`: updated with background context labels and materialized
  verification.
- `PLANS.md`: updated so the next Branch B step uses explicit background
  context labels but does not hard-filter them.
- `docs/CURRENT_CONTEXT_INDEX.md`: updated to point to this handoff and report.
- `AGENTS.md`: not updated; hard constraints did not change.
- `REPO_MAP.md`: not updated; stable repository structure did not change.

## 10. Do Not Misread

- Do not treat `UNKNOWN_BACKGROUND_NO_NULL_SUPPORT` as benign.
- Do not treat `SHADOW_BACKGROUND_NO_NULL_SUPPORT` as a formal N0/matched
  negative gate.
- Do not convert these labels into hard filters without locked truth-safe
  ablation.
- Do not treat 2026-06-15 as TP/FN/FP evidence.
- Do not promote Branch B V2, Branch S, the gap2m overlay, or the shadow
  reference from this loop.
