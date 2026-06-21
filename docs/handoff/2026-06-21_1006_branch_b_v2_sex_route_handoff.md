# Handoff: Branch B V2 Sex-Route Refinement

Date: 2026-06-21 10:06 Asia/Shanghai

Status: `active_current_handoff`

## 1. Goal

Continue the current PGT-A / CNV-seq R&D mainline by refining Branch B V2
evidence/disposition after the gap2m benchmark. This loop separates
sex-chromosome candidates from autosomal Branch B positive-support review and
routes them to Branch S review.

This loop does not promote Branch B V2, Branch S/SCA, the gap2m overlay, the
shadow reference, or any P6/report package to final release.

## 2. Context Restored

Read before execution:

- `C:\Users\11217\.codex\attachments\b29e8565-0516-429a-91ff-991e2ad43c59\goal-objective.md`
- `docs/CURRENT_CONTEXT_INDEX.md`
- `docs/handoff/2026-06-21_0937_branch_b_v2_gap2m_benchmark_handoff.md`
- `AGENTS.md`
- `skills/conversation_handoff/SKILL.md`
- `skills/pgta_reference_modeling_analysis/SKILL.md`
- `CURRENT_STATE.md`
- `PLANS.md`
- `docs/reports/branch_b_v2_gap2m_benchmark_2026-06-21.md`
- `pgta/predict/branch_b/v2_classifier.py`
- `tests/unit/test_branch_b_v2_classifier.py`

Active reference:

```text
h_r0_shadow_ref_20260619
```

Active Branch A overlay:

```text
merge_gap_bp=2_000_000
strong_z=10.0
branch_a.output_dir=wisecondorx/cnv/a_branch_gap2m
postprocess.output_dir=wisecondorx/cnv/postprocess_gap2m
```

Context conflict resolved:

- `docs/CURRENT_CONTEXT_INDEX.md` previously still pointed
  `active_handoff` at the older `2026-06-21_0318` handoff.
- This loop updates it to this handoff so future context restoration resumes
  from the sex-route refinement rather than the older Branch A-only state.

## 3. Completed Changes

Code:

- Added `is_sex_chromosome_candidate()` in
  `pgta/predict/branch_b/v2_classifier.py`.
- Updated `classify_candidate_row()` so `chrX`/`chrY` candidates route to:

```text
v2_candidate_class=V2_SEX_CHROMOSOME_REVIEW
v2_classifier_action=V2_ROUTE_BRANCH_S_REVIEW
v2_classifier_reason=sex_chromosome_branch_s_review:<original evidence tier>
```

Preserved behavior:

- original evidence tier;
- evidence gate;
- review priority;
- `v2_final_report_impact=none_shadow_only`;
- WisecondorX predict;
- Branch A candidate generation;
- mosaic logic;
- sex calling;
- final SCA/report status.

Tests:

- Added `test_sex_chromosome_candidate_routes_to_branch_s_review_not_autosomal_positive_support`.

Docs:

- Added `docs/reports/branch_b_v2_sex_route_refinement_2026-06-21.md`.
- Updated `docs/reports/branch_b_v2_gap2m_benchmark_2026-06-21.md`.
- Updated `CURRENT_STATE.md`.
- Updated `PLANS.md`.
- Updated `docs/CURRENT_CONTEXT_INDEX.md`.

## 4. Remote Materialization

Remote mirror:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419
```

Forced materialization commands were run on `ssh fengxian` for all three gap2m
benchmark configs:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile <branch_b_v2_gap2m_config> --cores 24 --forcerun cnv_branch_b_v2_classifier cnv_branch_b_v2_benchmark branch_b_v2_benchmark
```

Background jobs:

| cohort | PID | log |
|---|---:|---|
| Y1-Y8 | 10950 | `logs/branch_b_v2_sex_route_y_20260621.log` |
| H1-H16 | 10951 | `logs/branch_b_v2_sex_route_h_20260621.log` |
| 2026-06-15 | 10952 | `logs/branch_b_v2_sex_route_0615_20260621.log` |

Completion evidence:

- Y1-Y8 log: `12 of 12 steps (100%) done`.
- H1-H16 log: `20 of 20 steps (100%) done`.
- 2026-06-15 log: `9 of 9 steps (100%) done`.

Resource check during completion review:

- load average about `1.99`;
- memory used about `4.6G`;
- no abnormal saturation observed.

## 5. Materialized Result Summary

Class-level burden after sex routing:

| cohort | candidates | V2 positive-support class | V2 sex-chromosome review | V2 no-call contract risk | truth preserved | FN | hard-suppressed truth |
|---|---:|---:|---:|---:|---:|---:|---:|
| Y1-Y8 | 97 | 78 | 13 | 6 | 10/10 | 0 | 0 |
| H1-H16 | 105 | 57 | 42 | 6 | 10/10 | 0 | 0 |
| H7-H16 subset | 57 | 24 | 31 | 2 | no truth | NA | NA |
| 2026-06-15 | 165 | 127 | 14 | 24 | no truth | NA | NA |

Truth notes:

- Y3 chrX loss routes to `V2_SEX_CHROMOSOME_REVIEW` and remains preserved.
- H5/H6 chrX loss routes to `V2_SEX_CHROMOSOME_REVIEW` and remains preserved.
- H6 chr21 gain remains autosomal `V2_POSITIVE_SUPPORT_REVIEW` with
  `top_a_abs_zscore=7.1135`.
- 2026-06-15 has no locked truth table and remains burden/context only.

Important interpretation:

- Existing benchmark summary field `v2_positive_support_candidate_count` is
  evidence-tier based and may include sex-chromosome positive-support evidence.
- For autosomal Branch B burden, read classifier class counts from
  `v2_classifier/*.candidate_classification.tsv` or
  `docs/reports/branch_b_v2_sex_route_refinement_2026-06-21.md`.

## 6. Verification

Remote unit tests:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 /biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest tests/unit/test_branch_b_v2_classifier.py tests/unit/test_branch_b_v2_benchmark.py tests/unit/test_branch_ab_phase12_workflow_contract.py -q'
```

Result:

```text
28 passed in 0.95s
```

Remote single-test verification:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 /biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest tests/unit/test_branch_b_v2_classifier.py::test_v2_summary_reports_review_burden_without_report_promotion -q'
```

Result:

```text
1 passed in 0.59s
```

Remote materialized result audit:

- Parsed `v2_benchmark/summary.json`.
- Parsed all `v2_classifier/*.candidate_classification.tsv`.
- Parsed Y/H `v2_benchmark/truth_metrics.tsv`.

Result:

- truth preservation remains 10/10 for Y1-Y8 and H1-H6;
- FN=0;
- hard-suppressed truth=0;
- sex-chromosome candidates route to Branch S review.

## 7. Current Conclusion

This loop completes a safe Branch B V2 refinement:

- autosomal Branch B positive-support burden is now separated from
  sex-chromosome review burden;
- locked Y/H truth preservation remains intact;
- final-report impact remains `none_shadow_only`;
- Branch S receives the routed sex-chromosome review signal, but is still not
  final SCA.

This does not solve all FP/review burden. H7-H16 still has 24 autosomal
positive-support review rows, and 2026-06-15 still has 127 autosomal
positive-support review rows. Those remain the next Branch B V2 refinement
target.

## 8. Suggested Next Step

Proceed to the next Branch B V2 burden gate:

1. Keep the sex-route contract fixed.
2. Analyze remaining autosomal `V2_POSITIVE_SUPPORT_REVIEW` rows in H7-H16 and
   2026-06-15.
3. Separate broad whole-arm/whole-chromosome signals from smaller noisy events.
4. Add only review-label evidence that preserves Y/H truth and H6 chr21; do not
   introduce hard artifact decisions without locked ablation.

## 9. Key Files

Local code/docs:

- `pgta/predict/branch_b/v2_classifier.py`
- `tests/unit/test_branch_b_v2_classifier.py`
- `docs/reports/branch_b_v2_sex_route_refinement_2026-06-21.md`
- `docs/reports/branch_b_v2_gap2m_benchmark_2026-06-21.md`
- `CURRENT_STATE.md`
- `PLANS.md`
- `docs/CURRENT_CONTEXT_INDEX.md`

Remote result roots:

- `/data/project/CNV/PGT-A/refactor_validation_20260419/results_y_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m`
- `/data/project/CNV/PGT-A/refactor_validation_20260419/results_h_20260608_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m`
- `/data/project/CNV/PGT-A/refactor_validation_20260419/results_20260615_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m`

## 10. Core File Sync

- `CURRENT_STATE.md`: updated with materialized sex-route results.
- `PLANS.md`: updated so the next Branch B target is remaining autosomal review
  burden.
- `docs/CURRENT_CONTEXT_INDEX.md`: updated to point to this handoff and this
  report.
- `AGENTS.md`: not updated; hard constraints did not change.
- `REPO_MAP.md`: not updated; stable repository structure did not change.

## 11. Do Not Misread

- Do not treat this as Branch B V2 final promotion.
- Do not treat Branch S as final SCA.
- Do not treat sex-route review as a PASS/FAIL SCA call.
- Do not treat 2026-06-15 as TP/FN/FP evidence.
- Do not use legacy/current-code Branch B kept/artifact/final-disposition fields
  as V2 decision evidence.
- Do not promote gap2m to default until downstream A/B/S/report gates are
  benchmarked and accepted.
