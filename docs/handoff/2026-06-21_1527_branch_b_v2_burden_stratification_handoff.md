# Handoff: Branch B V2 Burden Stratification

Date: 2026-06-21 15:27 Asia/Shanghai

Status: `active_current_handoff`

## 1. Goal

Continue after the Branch B V2 truth-safe filter contract. The immediate goal
was to add burden stratification on top of the already materialized Branch A
`merge_gap_bp=2_000_000` overlay without changing Branch A, without using
legacy/current-code Branch B decisions, and without hard-suppressing locked
truth candidates.

Branch A remains WisecondorX/CBS-derived candidate discovery. Branch B V2
refines Branch A candidates only and does not create B-only report events.

## 2. Restored Project State

Context was restored from:

- `docs/CURRENT_CONTEXT_INDEX.md`
- `docs/handoff/2026-06-21_1454_branch_b_v2_truth_safe_filter_handoff.md`
- `AGENTS.md`
- `skills/conversation_handoff/SKILL.md`
- `skills/pgta_reference_modeling_analysis/SKILL.md`
- `CURRENT_STATE.md`
- `PLANS.md`

Active reference:

```text
h_r0_shadow_ref_20260619
```

Active Branch A overlay for Branch B V2 benchmarking:

```text
merge_gap_bp=2_000_000
strong_z=10.0
postprocess.output_dir=wisecondorx/cnv/postprocess_gap2m
```

Current working branch for this loop:

```text
codex/branch-b-v2-burden-stratification
```

## 3. Completed In This Loop

Branch B V2 classifier now emits burden-stratification fields:

- `v2_burden_reduction_version`
- `v2_burden_reduction_tier`
- `v2_burden_reduction_action`
- `v2_burden_reduction_reason`
- `v2_burden_evidence_tags`

Branch B V2 benchmark now records:

- top-candidate burden tier/action/reason/tag fields in truth metrics;
- sample-level burden counts for report/review/background/technical/Branch S;
- summary-level stratification counts for filter action, disposition, length
  tier, clean support, B-side signal context, GC/RC context, burden tier, and
  burden action.

New report:

- `docs/reports/branch_b_v2_burden_stratification_2026-06-21.md`

Core files updated:

- `CURRENT_STATE.md`
- `PLANS.md`
- `docs/CURRENT_CONTEXT_INDEX.md`

## 4. CNVseq / CNVpro Use

The following borrowed concepts are explicitly tagged and limited to evidence
context:

- `[CNVpro-inspired]` length tiers for report/review routing.
- `[CNVpro-confirmed]` acrocentric qter context for chr13/14/15/21/22 review.
- `[CNVseq-asset]` mask/mappability/PAR/sex-homology annotation context.
- `[CNVpro-like]` GC/RC context.
- `[Not used]` CNVcalling.R/cghFLasso as primary caller replacements.

None of these tags replaces WisecondorX/CBS, changes Branch A discovery, or
acts as a hard artifact filter.

## 5. Files Modified

- `pgta/predict/branch_b/v2_classifier.py`
- `pgta/predict/branch_b/v2_benchmark.py`
- `tests/unit/test_branch_b_v2_classifier.py`
- `tests/unit/test_branch_b_v2_benchmark.py`
- `tests/unit/test_current_context_index.py`
- `docs/reports/branch_b_v2_burden_stratification_2026-06-21.md`
- `CURRENT_STATE.md`
- `PLANS.md`
- `docs/CURRENT_CONTEXT_INDEX.md`
- `docs/handoff/2026-06-21_1527_branch_b_v2_burden_stratification_handoff.md`

## 6. Verification So Far

Remote TDD red check before implementation:

```text
5 failed, 25 passed
```

Remote unit validation after implementation:

```bash
cd /data/project/CNV/PGT-A/refactor_validation_20260419
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_v2_classifier.py \
  tests/unit/test_branch_b_v2_benchmark.py -q
```

Result:

```text
30 passed in 0.88s
```

Snakemake dry-runs succeeded for
`branch_b_v2_benchmark branch_s_review cnv_report` under all three active
gap2m configs:

- `config_predict_y_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`
- `config_predict_h_20260608_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`
- `config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`

Forced materialization completed for:

```text
cnv_branch_b_v2_classifier
cnv_branch_b_v2_benchmark
branch_b_v2_benchmark
```

Materialized remote benchmark result:

| cohort | candidates | truth events | truth preserved | FN | hard-suppressed truth | filter-suppressed candidates | background review | Branch S review |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Y1-Y8 | 97 | 10 | 10 | 0 | 0 | 0 | 84 | 13 |
| H1-H16 | 105 | 10 | 10 | 0 | 0 | 0 | 63 | 42 |
| 2026-06-15 | 165 | 0 | 0 | 0 | 0 | 0 | 151 | 14 |

H6 chr21 remains preserved:

```text
sample=H6
top_candidate_id=H6.A0003
top_a_abs_zscore=7.113507
top_v2_disposition=background_unknown_review
top_v2_burden_reduction_tier=background_unknown_review
top_v2_burden_reduction_action=stratify_background_unknown_review
```

## 7. Current Interpretation

This loop improves auditability and future report/review routing. It does not
yet prove total FP/review burden reduction, because current materialized
outputs remain only `background_unknown_review` and `branch_s_review`.

Branch B V2 remains development-only. Branch S remains
`review_reportable_with_limitations`, not final SCA.

## 8. Next Recommended Step

Proceed to burden-display and report-contract integration. Do not add hard
filters unless a locked-truth ablation shows Y1-Y8/H1-H6 truth preservation,
explicit H6 chr21 preservation, and no regression in weak-positive /
mosaic-sensitive events.

H7-H16 and 2026-06-15 should remain context/burden cohorts unless locked truth
labels are added.
