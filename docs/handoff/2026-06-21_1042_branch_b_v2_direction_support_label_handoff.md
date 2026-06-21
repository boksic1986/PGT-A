# Handoff: Branch B V2 Direction-Support Review Label

Date: 2026-06-21 10:42 Asia/Shanghai

Status: `active_current_handoff`

## 1. Goal

Continue Branch B V2 evidence/disposition refinement by implementing the safe
output from the autosomal burden audit: Branch B-side direction support is now a
review label only, not a hard filter or final-disposition rule.

## 2. Context Restored

Read before execution:

- `C:\Users\11217\.codex\attachments\b29e8565-0516-429a-91ff-991e2ad43c59\goal-objective.md`
- `docs/CURRENT_CONTEXT_INDEX.md`
- `docs/handoff/2026-06-21_1021_branch_b_v2_autosomal_burden_audit_handoff.md`
- `AGENTS.md`
- `skills/conversation_handoff/SKILL.md`
- `skills/pgta_reference_modeling_analysis/SKILL.md`
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
v2_direction_support_label
v2_direction_support_reason
```

These fields are emitted in `v2_classifier/*.candidate_classification.tsv` and
summarized as `direction_support_label_counts` in each V2 classifier summary.

The implementation intentionally does not alter:

- WisecondorX predict;
- Branch A candidate discovery;
- `v2_candidate_class`;
- `v2_classifier_action`;
- `v2_final_report_impact`;
- hard-suppression behavior;
- Branch S status.

New report:

- `docs/reports/branch_b_v2_direction_support_label_2026-06-21.md`

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
3 failed, 11 passed
```

The failures were expected missing-field failures for:

- `v2_direction_support_label`;
- `direction_support_label_counts`.

Remote unit test after implementation:

```text
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_v2_classifier.py \
  tests/unit/test_branch_b_v2_benchmark.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py -q
```

Result:

```text
32 passed in 0.71s
```

Additional TDD red check before the final conflict-priority fix:

```text
1 failed
```

The failure confirmed that amplitude direction conflict must take priority over
fraction-based support in the review label.

Remote forced materialization used:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake
```

Targets/rules rerun:

```text
cnv_branch_b_v2_classifier
cnv_branch_b_v2_benchmark
branch_b_v2_benchmark
```

Configs:

- `config_predict_y_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`
- `config_predict_h_20260608_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`
- `config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`

All three forced materializations completed.

## 5. Current Evidence

Direction label counts:

| cohort | rows | B direction supported | A-only weak B direction | B direction conflict |
|---|---:|---:|---:|---:|
| Y1-Y8 truth cohort | 97 | 66 | 20 | 11 |
| H1-H16 | 105 | 68 | 26 | 11 |
| 2026-06-15 context | 165 | 97 | 40 | 28 |

Autosomal `V2_POSITIVE_SUPPORT_REVIEW` rows:

| cohort | autosomal positive rows | B direction supported | A-only weak B direction |
|---|---:|---:|---:|
| Y1-Y8 truth cohort | 78 | 58 | 20 |
| H1-H16 | 57 | 35 | 22 |
| 2026-06-15 context | 127 | 89 | 38 |

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
v2_direction_support_label=B_DIRECTION_SUPPORTED
v2_final_report_impact=none_shadow_only
top_a_abs_zscore=7.113507302991461
```

## 6. Current Conclusion

Direction support is now represented in workflow outputs as review evidence.

It is still not a Branch B V2 promotion gate, not a hard filter, and not final
report evidence. Weak Branch B-side direction support remains unsafe as a
universal downgrade because known autosomal truth candidates can be
`A_ONLY_WEAK_B_DIRECTION`.

The main unresolved Branch B V2 condition remains:

```text
UNKNOWN_BACKGROUND + NO_NULL_SUPPORT
```

## 7. Suggested Next Step

Move to the next Branch B V2 refinement only after preserving the current
direction-label contract:

1. Keep direction support as review evidence only.
2. Continue excluding legacy/current-code Branch B decision fields.
3. Target `UNKNOWN_BACKGROUND + NO_NULL_SUPPORT` with a truth-safe
   background/evidence contract.
4. Preserve Y1-Y8/H1-H6 truth 10/10, H6 chr21, and hard-suppressed truth=0.

## 8. Key Files

- `pgta/predict/branch_b/v2_classifier.py`
- `tests/unit/test_branch_b_v2_classifier.py`
- `docs/reports/branch_b_v2_direction_support_label_2026-06-21.md`
- `CURRENT_STATE.md`
- `PLANS.md`
- `docs/CURRENT_CONTEXT_INDEX.md`

Remote result roots:

- `/data/project/CNV/PGT-A/refactor_validation_20260419/results_y_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m`
- `/data/project/CNV/PGT-A/refactor_validation_20260419/results_h_20260608_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m`
- `/data/project/CNV/PGT-A/refactor_validation_20260419/results_20260615_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m`

## 9. Core File Sync

- `CURRENT_STATE.md`: updated with implemented direction-support review label
  and materialized verification.
- `PLANS.md`: updated so the next Branch B step targets background/evidence
  contracts, not direction hard filtering.
- `docs/CURRENT_CONTEXT_INDEX.md`: updated to point to this handoff and report.
- `AGENTS.md`: not updated; hard constraints did not change.
- `REPO_MAP.md`: not updated; stable repository structure did not change.

## 10. Do Not Misread

- Do not treat `B_DIRECTION_SUPPORTED` as a PASS call.
- Do not treat `A_ONLY_WEAK_B_DIRECTION` as benign/artifact.
- Do not hard-filter direction conflicts solely from this label.
- Do not treat 2026-06-15 as TP/FN/FP evidence.
- Do not promote Branch B V2, Branch S, the gap2m overlay, or the shadow
  reference from this loop.
