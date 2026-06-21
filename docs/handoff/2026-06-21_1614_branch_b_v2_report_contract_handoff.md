# Handoff: Branch B V2 Report Contract Integration

Date: 2026-06-21 16:14 Asia/Shanghai

Status: `active_current_handoff`

## 1. Goal

Connect the already materialized Branch B V2 burden stratification outputs to
`cnv_report` as `development_review_only` report-display evidence.

This does not change Branch A, does not add hard filters, and does not promote
Branch B V2 or Branch S.

## 2. Restored Project State

Context was restored from:

- `docs/CURRENT_CONTEXT_INDEX.md`
- `docs/handoff/2026-06-21_1527_branch_b_v2_burden_stratification_handoff.md`
- `AGENTS.md`
- `skills/conversation_handoff/SKILL.md`
- `skills/pgta_reference_modeling_analysis/SKILL.md`
- `CURRENT_STATE.md`
- `PLANS.md`

Active reference:

```text
h_r0_shadow_ref_20260619
```

Active Branch A overlay for this report-contract loop:

```text
merge_gap_bp=2_000_000
strong_z=10.0
postprocess.output_dir=wisecondorx/cnv/postprocess_gap2m
```

## 3. Completed In This Loop

Completed code changes:

- `pgta/predict/report.py` loads Branch B V2 benchmark summary and
  sample-summary inputs.
- `rules/predict_workflow.smk` wires V2 benchmark outputs into
  `cnv_report_summary` when `branch_b_v2_benchmark` is requested.
- `cnv_summary.tsv/json/md/html` exposes Branch B V2 burden/status fields as
  development-only evidence.

Report contract fields:

- `branch_b_v2_burden_status`
- `branch_b_v2_background_unknown_review_count`
- `branch_b_v2_branch_s_review_count`
- `branch_b_v2_technical_risk_review_count`
- `branch_b_v2_report_candidate_count`
- `branch_b_v2_legacy_fields_used=false`
- `branch_b_v2_final_impact=development_review_only`

## 4. Current Interpretation

This loop is report-contract integration only. It does not prove FP/review
burden reduction and does not make Branch B V2 final-report logic.

Branch S remains `review_reportable_with_limitations`.

## 5. Verification

Remote validation on `ssh fengxian`:

```bash
cd /data/project/CNV/PGT-A/refactor_validation_20260419
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_cnv_report.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py \
  tests/unit/test_current_context_index.py -q
```

Result:

```text
28 passed in 0.93s
```

Snakemake dry-run passed for
`branch_b_v2_benchmark branch_s_review cnv_report` under all three active
gap2m configs:

- `config_predict_y_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`
- `config_predict_h_20260608_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`
- `config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`

Materialized `cnv_report` for all three active gap2m configs.

Materialized contract checks:

- Y1-Y8: V2 truth preserved 10/10, FN=0, hard-suppressed truth=0.
- H1-H16: V2 truth preserved 10/10, FN=0, hard-suppressed truth=0; H6 chr21
  remains visible in the report.
- 2026-06-15: V2 status is `skipped_no_truth`; report is burden/context only.
- All three report JSON contracts keep
  `branch_b_v2_final_impact=development_review_only` and
  `branch_b_v2_legacy_fields_used=false`.

## 6. Core File Sync

- `CURRENT_STATE.md`: updated with report-contract integration status.
- `PLANS.md`: updated with report-contract validation/materialization gate.
- `docs/CURRENT_CONTEXT_INDEX.md`: updated to point to this handoff and report.
- `AGENTS.md`: not updated; no repository-level hard rule changed.
