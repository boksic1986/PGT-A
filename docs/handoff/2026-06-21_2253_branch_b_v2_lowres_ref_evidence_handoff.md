---
status: active_handoff
decision_use: current_context
created_at: 2026-06-21 22:53 Asia/Shanghai
active_reference_id: h_r0_shadow_ref_20260619
branch_a_overlay: merge_gap_bp_2000000
final_impact: development_review_only
---

# Branch B V2 Low-Resolution Ref Evidence Handoff

## 1. 本轮目标

Implement the first low-resolution evidence contract for Branch B V2:

- 2Mb/3Mb same-direction evidence as auxiliary context;
- reference-bin MAD stability context;
- no Branch A change;
- no low-res standalone filtering;
- no Branch B V2 or Branch S promotion.

This loop implements the interface and workflow wiring only. It does not
materialize 2Mb/3Mb shadow references.

## 2. 已确认的项目状态

- Current context entrypoint:
  `docs/CURRENT_CONTEXT_INDEX.md`.
- Previous handoff:
  `docs/handoff/2026-06-21_2123_branch_b_v2_pass2_correction_handoff.md`.
- Active reference remains `h_r0_shadow_ref_20260619`, still
  `fixed_shadow_baseline_not_production`.
- Active Branch A benchmark input remains the explicit
  `merge_gap_bp=2_000_000` overlay; default Branch A behavior is unchanged.
- Branch B V2 remains Branch-A-anchored candidate evidence/refinement only.
- Branch S remains `review_reportable_with_limitations`.

## 3. 已完成事项

Implemented:

- `pgta/predict/branch_b/lowres_evidence.py`
- `pgta/predict/branch_b/ref_stability.py`
- low-res/ref-MAD context fields in
  `pgta/predict/branch_b/v2_classifier.py`
- dispatcher actions in `scripts/_compat_entry.py`
- Snakemake entrypoint constants in `rules/script_entrypoints.smk`
- config-gated low-res paths and fail-fast checks in
  `rules/predict_layout.smk`
- optional low-res/ref-MAD workflow rules in `rules/predict_workflow.smk`
- low-res targets in `rules/target_assembly.smk`
- unit tests for low-res evidence, ref-MAD stability, classifier safety, and
  workflow contract.

New classifier fields:

```text
v2_lowres_context_label
v2_lowres_context_reason
v2_ref_stability_context
```

New optional low-res ledger fields include:

```text
lowres_2mb_support_label
lowres_2mb_same_direction
lowres_2mb_overlap_fraction
lowres_2mb_z_or_score
lowres_3mb_support_label
lowres_3mb_same_direction
lowres_3mb_overlap_fraction
lowres_3mb_z_or_score
lowres_consensus_label
event_ref_mad_median
event_ref_mad_p90
high_ref_mad_bin_fraction
ref_stability_context
```

## 4. 当前结论

Low-resolution evidence is now an implemented optional Branch B V2 context
interface.

It is not yet a validated biological filtering layer because the actual 2Mb and
3Mb shadow references and low-res predict outputs have not been built in this
loop.

Current safety contract:

- low-res same-direction support can increase review confidence;
- low-res absence is never a single filter;
- H6 chr21-like short/weak positives remain at least internal review;
- high ref-MAD weakens negative low-res interpretation;
- chrX/chrY/SCA candidates remain Branch S context, not autosomal low-res
  filter inputs.

## 5. 待验证问题

- Build `h_r0_shadow_ref_20260619_2mb`.
- Build `h_r0_shadow_ref_20260619_3mb`.
- Rerun low-res WisecondorX predict for Y1-Y8, H1-H16, G1-G8, and
  2026-06-15.
- Materialize low-res evidence ledgers.
- Measure low-res same-direction support among locked Y/H/G truth events,
  especially events `>=4Mb`.
- Confirm H6 chr21 remains visible and is not filtered by low-res absence.
- Produce Branch S burden table separately for chrX/chrY/SCA.

## 6. 建议下一步

Start a remote long-task gate to build and predict the 2Mb/3Mb shadow
references. Use explicit configs and return PID/log/target paths immediately
instead of blocking until completion.

Suggested gate output:

- low-res support table for all locked truth candidates;
- low-res support table for 2026-06-15 report/internal-review candidates;
- ref-MAD stability table for candidate regions;
- Branch S burden table.

## 7. 关键文件与路径

Report:

- `docs/reports/branch_b_v2_lowres_ref_evidence_2026-06-21.md`

Code:

- `pgta/predict/branch_b/lowres_evidence.py`
- `pgta/predict/branch_b/ref_stability.py`
- `pgta/predict/branch_b/v2_classifier.py`

Workflow:

- `rules/predict_layout.smk`
- `rules/predict_workflow.smk`
- `rules/script_entrypoints.smk`
- `rules/target_assembly.smk`
- `scripts/_compat_entry.py`

Tests:

- `tests/unit/test_lowres_evidence.py`
- `tests/unit/test_ref_stability.py`
- `tests/unit/test_branch_b_v2_classifier.py`
- `tests/unit/test_branch_ab_phase12_workflow_contract.py`

## 8. 环境约束与注意事项

- All real validation must run on `ssh fengxian`.
- Remote mirror:
  `/data/project/CNV/PGT-A/refactor_validation_20260419`.
- Python:
  `/biosoftware/miniconda/envs/snakemake_env/bin/python`.
- Snakemake:
  `/biosoftware/miniconda/envs/snakemake_env/bin/snakemake`.
- Do not use local test results as pass evidence.
- Do not stage local temp scripts, Excel files, local `reports/`, images, or
  `validate_npz.py`.

## 9. 命令执行记录

Remote validation executed after syncing the changed files to the mirror:

```text
cd /data/project/CNV/PGT-A/refactor_validation_20260419
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_lowres_evidence.py \
  tests/unit/test_ref_stability.py \
  tests/unit/test_branch_b_v2_classifier.py \
  tests/unit/test_branch_b_v2_benchmark.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py \
  tests/unit/test_cnv_report.py \
  tests/unit/test_current_context_index.py -q
```

Result: pending final rerun after this handoff is synced.
Result:

```text
81 passed in 1.69s
```

Default-path dry-runs to run:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile <active gap2m config> --cores 1 -n \
  branch_b_v2_benchmark branch_s_review cnv_report
```

Default-path dry-run result: RC=0 for all four active gap2m configs:

```text
config_predict_y_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml
config_predict_h_20260608_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml
config_predict_g1_g8_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml
config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml
```

All four returned `Nothing to be done`.

A temporary remote lowres-enabled dry-run config was created at:

```text
/tmp/pgta_lowres_evidence_dryrun_20260621.yaml
```

It used two existing reference NPZ inputs only to validate rule wiring. Result:
RC=0. Planned jobs:

```text
cnv_branch_b_ref_stability_bins: 1
cnv_branch_b_ref_stability: 8
cnv_branch_b_lowres_evidence: 8
cnv_branch_b_v2_classifier: 8
cnv_branch_b_v2_benchmark: 1
branch_b_v2_benchmark: 1
total: 27
```

No lowres outputs were materialized; this was dry-run only.

## 10. Core file sync

- `docs/CURRENT_CONTEXT_INDEX.md`: updated to this handoff.
- `CURRENT_STATE.md`: updated with low-res/ref-MAD interface status.
- `PLANS.md`: updated with low-res/ref-MAD next gate.
- `REPO_MAP.md`: not updated; no stable directory-level ownership change.
- `AGENTS.md`: not updated; no hard-rule change.
