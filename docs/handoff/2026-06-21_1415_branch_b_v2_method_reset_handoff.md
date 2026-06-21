# Handoff: Branch B V2 Method Reset And Threshold Inventory

Date: 2026-06-21 14:15 Asia/Shanghai

Status: `active_current_handoff`

## 1. Goal

Reset Branch B V2 away from legacy/current-code Branch B threshold filtering and
toward Branch-A-anchored evidence stratification.

Branch A remains the WisecondorX/CBS-derived primary discovery layer. Branch B
V2 does not generate B-only final events and does not use GC/RC corrected signal
direction to contradict Branch A.

## 2. Context Restored

Read before execution:

- `C:\Users\11217\.codex\attachments\b29e8565-0516-429a-91ff-991e2ad43c59\goal-objective.md`
- `docs/CURRENT_CONTEXT_INDEX.md`
- `docs/handoff/2026-06-21_1103_branch_b_v2_background_context_label_handoff.md`
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

Updated Branch B V2 classifier outputs:

```text
v2_signal_strength_tier
v2_length_tier
v2_clean_support_label
v2_gc_rc_context_label
v2_b_signal_context_label
v2_b_signal_context_reason
v2_disposition
```

Updated semantics:

- `B_DIRECTION_CONFLICT` is replaced by
  `B_SIGNAL_DISCORDANT_WITH_A_DIRECTION`.
- `B_SIGNAL_DISCORDANT_WITH_A_DIRECTION` is review context only. It does not
  mean Branch B called an opposite CNV event.
- GC/RC context is auxiliary evidence only and cannot hard-suppress Branch A
  positive candidates.
- First-round V2 dispositions are conservative:
  `report_candidate`, `review_candidate`, `technical_risk_review`,
  `background_unknown_review`, and `sca_branch_s_review`.

Updated benchmark output:

- `truth_metrics.tsv` now records top-candidate V2 disposition, length tier,
  clean-support label, GC/RC context label, B-signal context label, and
  attenuation ratio.

New report:

- `docs/reports/branch_b_v2_method_reset_threshold_inventory_2026-06-21.md`

Core file updates:

- `CURRENT_STATE.md`
- `PLANS.md`
- `docs/CURRENT_CONTEXT_INDEX.md`

## 4. Code And Tests Changed

- `pgta/predict/branch_b/v2_classifier.py`
- `pgta/predict/branch_b/v2_benchmark.py`
- `tests/unit/test_branch_b_v2_classifier.py`
- `tests/unit/test_branch_b_v2_benchmark.py`

## 5. Verification

Remote TDD red check before implementation:

```text
5 failed, 17 passed
```

The failures were expected missing-field / old-semantics failures.

Remote unit test after implementation:

```bash
cd /data/project/CNV/PGT-A/refactor_validation_20260419
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_v2_classifier.py \
  tests/unit/test_branch_b_v2_benchmark.py -q
```

Result:

```text
22 passed in 0.77s
```

Broader remote unit test:

```bash
cd /data/project/CNV/PGT-A/refactor_validation_20260419
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_v2_classifier.py \
  tests/unit/test_branch_b_v2_benchmark.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py \
  tests/unit/test_cnv_report.py -q
```

Result:

```text
45 passed in 1.26s
```

Remote Snakemake dry-runs succeeded for:

```text
branch_b_v2_benchmark branch_s_review cnv_report
```

Configs:

- `config_predict_y_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`
- `config_predict_h_20260608_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`
- `config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`

Forced materialization completed for all three configs:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake \
  -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile \
  --configfile <config> \
  --directory /data/project/CNV/PGT-A/refactor_validation_20260419 \
  --cores 8 \
  --forcerun cnv_branch_b_v2_classifier cnv_branch_b_v2_benchmark \
  branch_b_v2_benchmark
```

## 6. Current Evidence

Benchmark summary:

| cohort | samples | candidates | truth events | truth preserved | FN | hard-suppressed truth | status |
|---|---:|---:|---:|---:|---:|---:|---|
| Y1-Y8 | 8 | 97 | 10 | 10 | 0 | 0 | ready |
| H1-H16 | 16 | 105 | 10 | 10 | 0 | 0 | ready |
| 2026-06-15 | 5 | 165 | 0 | 0 | 0 | 0 | skipped_no_truth |

H6 chr21 remains retained:

```text
top_candidate_id=H6.A0003
top_a_abs_zscore=7.113507302991461
top_v2_disposition=background_unknown_review
top_v2_b_signal_context_label=B_SIGNAL_SUPPORTED_A_DIRECTION
top_v2_gc_rc_context_label=GC_RC_ATTENUATED
```

Materialized V2 disposition counts:

- Y1-Y8: `background_unknown_review=84`, `sca_branch_s_review=13`.
- H1-H16: `background_unknown_review=63`, `sca_branch_s_review=42`.
- 2026-06-15: `background_unknown_review=151`, `sca_branch_s_review=14`.

0615 per-sample burden:

| sample | background_unknown_review | sca_branch_s_review |
|---|---:|---:|
| JZ26125843-56-56 | 24 | 2 |
| JZ26125844-59-59 | 34 | 1 |
| JZ26125845-60-60 | 42 | 5 |
| JZ26125846-61-61 | 30 | 3 |
| JZ26125847-62-62 | 21 | 3 |

## 7. Current Conclusion

Branch B V2 is now more interpretable and truth-safe:

- locked Y/H truth is preserved;
- H6 chr21 is preserved;
- legacy/current-code Branch B decision fields are still ignored;
- GC/RC attenuation/discordance cannot hard-suppress a Branch A positive;
- 0615 remains context-only without TP/FN/FP claims.

This loop does not solve FP/review burden and does not promote Branch B V2,
Branch S, the gap2m overlay, or the shadow reference to final release.

## 8. Suggested Next Step

Next Branch B V2 loop:

1. Use `v2_disposition`, `v2_length_tier`, `v2_clean_support_label`,
   `v2_gc_rc_context_label`, and `v2_b_signal_context_label` to design
   truth-safe burden stratification.
2. Keep Y1-Y8/H1-H6 truth 10/10, FN=0, hard-suppressed truth=0.
3. Keep H6 chr21 preserved.
4. Do not convert background context, B-signal discordance, or length tier into
   hard filters without locked-truth ablation.
5. Continue Branch S separately toward `review_reportable_with_limitations`.

## 9. Do Not Misread

- Do not treat this loop as FP reduction.
- Do not treat `background_unknown_review` as benign.
- Do not treat `B_SIGNAL_DISCORDANT_WITH_A_DIRECTION` as an opposite CNV call.
- Do not use 0615 for TP/FN/FP performance.
- Do not use legacy/current-code Branch B fields as V2 evidence.
- Do not promote Branch B V2, Branch S, the gap2m overlay, or the shadow ref.
