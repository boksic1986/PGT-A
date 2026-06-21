# Handoff: Branch B V2 Truth-Safe Filter Contract

Date: 2026-06-21 14:54 Asia/Shanghai

Status: `active_current_handoff`

## 1. Goal

Continue after merging `codex/branch-a-burden-merge-gap` into `main` and
synchronizing the remote mirror. The immediate goal is to start the next Branch
B V2 loop by making filter actions explicit, truth-safe, and separate from
legacy/current-code Branch B thresholds.

Branch A remains the WisecondorX/CBS-derived candidate discovery layer. Branch
B V2 refines Branch A candidates only and does not create B-only report events.

## 2. Restored Project State

Context was restored from:

- `docs/CURRENT_CONTEXT_INDEX.md`
- `docs/handoff/2026-06-21_1415_branch_b_v2_method_reset_handoff.md`
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

## 3. Completed Before This Handoff

The previous branch was closed out:

- feature branch: `codex/branch-a-burden-merge-gap`;
- merged into `main` by fast-forward;
- pushed `main` to `origin`;
- synchronized the remote mirror at
  `/data/project/CNV/PGT-A/refactor_validation_20260419`;
- verified mirror file hashes for `docs/CURRENT_CONTEXT_INDEX.md`,
  `tests/unit/test_current_context_index.py`,
  `pgta/predict/branch_b/v2_classifier.py`, and `Snakefile`;
- remote post-sync unit tests passed:
  `51 passed in 1.80s`.

Current working branch for this loop:

```text
codex/branch-b-v2-truth-safe-filter
```

## 4. Completed In This Loop

Branch B V2 classifier now emits explicit filter-action fields:

- `v2_filter_version`;
- `v2_filter_action`;
- `v2_filter_reason`;
- `v2_filter_scope`;
- `v2_filter_hard_suppression_allowed`.

Branch B V2 benchmark now records:

- `top_v2_filter_action` in truth-overlap metrics;
- `v2_filter_suppressed_count` per sample;
- `v2_filter_suppressed_candidate_count` in summary JSON.

Current filter contract:

- `suppress_workflow_contract_risk` is the only hard-suppression action.
- The current hard-suppression reason is
  `same_chrom_ref_leakage_contract_risk`.
- GC/RC attenuation, B-side signal discordance, unknown background, length
  tier, low clean support, and high-risk-region burden cannot hard-suppress a
  Branch A candidate.
- Low clean support/high-risk burden can only downgrade to
  `technical_risk_review`.
- Sex-chromosome candidates route to Branch S review and do not finalize SCA.

New report:

- `docs/reports/branch_b_v2_truth_safe_filter_2026-06-21.md`

Core files updated:

- `CURRENT_STATE.md`
- `PLANS.md`
- `docs/CURRENT_CONTEXT_INDEX.md`

## 5. Files Modified

- `pgta/predict/branch_b/v2_classifier.py`
- `pgta/predict/branch_b/v2_benchmark.py`
- `tests/unit/test_branch_b_v2_classifier.py`
- `tests/unit/test_branch_b_v2_benchmark.py`
- `tests/unit/test_current_context_index.py`
- `docs/reports/branch_b_v2_truth_safe_filter_2026-06-21.md`
- `CURRENT_STATE.md`
- `PLANS.md`
- `docs/CURRENT_CONTEXT_INDEX.md`
- `docs/handoff/2026-06-21_1454_branch_b_v2_truth_safe_filter_handoff.md`

## 6. Verification So Far

Remote TDD red check before implementation:

```text
9 failed, 17 passed
```

Remote unit validation after implementation, before this handoff was written:

```bash
cd /data/project/CNV/PGT-A/refactor_validation_20260419
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_v2_classifier.py \
  tests/unit/test_branch_b_v2_benchmark.py -q
```

Result:

```text
26 passed in 0.90s
```

Broader remote unit tests after synchronizing this handoff and context index:

```bash
cd /data/project/CNV/PGT-A/refactor_validation_20260419
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_v2_classifier.py \
  tests/unit/test_branch_b_v2_benchmark.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py \
  tests/unit/test_cnv_report.py \
  tests/unit/test_branch_a_candidates.py \
  tests/unit/test_current_context_index.py -q
```

Result:

```text
55 passed in 1.37s
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

| cohort | candidates | truth events | truth preserved | FN | hard-suppressed truth | filter-suppressed candidates | status |
|---|---:|---:|---:|---:|---:|---:|---|
| Y1-Y8 | 97 | 10 | 10 | 0 | 0 | 0 | ready |
| H1-H16 | 105 | 10 | 10 | 0 | 0 | 0 | ready |
| 2026-06-15 | 165 | 0 | 0 | 0 | 0 | 0 | skipped_no_truth |

Remote result-file check confirmed:

- classifier TSV contains `v2_filter_action`;
- truth metrics TSV contains `top_v2_filter_action`;
- sample summary TSV contains `v2_filter_suppressed_count`.

## 7. Current Conclusion

The filter-action contract is implemented and materialized on the active gap2m
benchmark path. It remains a truth-safe contract validation, not an FP-reduction
or final-report promotion claim.

Current interpretation:

- Branch B V2 is still `none_shadow_only`.
- Branch B V2 has not yet proven FP/review burden reduction.
- Branch S remains `review_reportable_with_limitations`, not final SCA.
- The gap2m overlay and `h_r0_shadow_ref_20260619` remain shadow/development
  assets, not production promotion.

## 8. Suggested Next Step

Commit this loop, merge/sync `main`, then continue the next Branch B V2 burden
reduction step using the now-materialized filter-action fields.

Promotion criteria for the next loop remain limited:

- Y1-Y8 truth preserved 10/10;
- H1-H6 truth preserved 10/10;
- H6 chr21 preserved;
- hard-suppressed truth remains 0;
- 0615 remains burden/context only.

## 9. Environment Constraints

All true validation must run on:

```text
ssh fengxian
```

Remote mirror:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419
```

Locked executables:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/python
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake
```

Do not use local test results as final validation evidence.

## 10. Core File Sync

- `CURRENT_STATE.md`: updated
- `PLANS.md`: updated
- `docs/CURRENT_CONTEXT_INDEX.md`: updated
- `REPO_MAP.md`: not updated; no repo entrypoint or directory structure change
- `AGENTS.md`: not updated; no hard-rule change
