# 2026-06-19 12:25 Phase 4 Unknown Background Handoff

## 1. Goal

Fix and document a Branch B V2 shadow-classifier contract gap where Phase 1 fallback ledgers with `matched_negative_source=UNKNOWN_BACKGROUND` were not treated the same as Phase 3 matched-negative rows with `matched_negative_background_status=UNKNOWN_BACKGROUND`.

This run did not change WisecondorX predict, Branch A candidate generation, legacy Branch B final decisions, Branch S, or final report schema.

## 2. Context Source

- Latest active handoff before this run:
  `docs/handoff/2026-06-19_1201_phase1b_asset_decision_matrix_handoff.md`
- Repository constraints read:
  `AGENTS.md`
- Skills read:
  `skills/conversation_handoff/SKILL.md`
  `skills/pgta_reference_modeling_analysis/SKILL.md`
- Active plan document:
  `docs/reports/branch_ab_v2_rnd_constraints_2026-06-18.md`
- Main implementation files:
  `pgta/predict/branch_b/v2_classifier.py`
  `tests/unit/test_branch_b_v2_classifier.py`

## 3. Completed

- Added a regression test:
  `test_phase1_unknown_matched_negative_source_forces_review_not_artifact`.
- Updated V2 classifier logic so unknown matched-negative background is detected from either:
  - `matched_negative_background_status=UNKNOWN_BACKGROUND`;
  - `matched_negative_source=UNKNOWN_BACKGROUND`.
- Refreshed V2 classifier outputs on `fengxian` for:
  - Y1-Y8 under `results_build_ref_v2_mask_only`;
  - H1-H16 under `results_h_20260608_mask_only`;
  - 2026-06-15 five-sample run under `results_20260615_mask_only`.
- Added report:
  `docs/reports/branch_ab_v2_phase4_unknown_background_fix_2026-06-19.md`
- Updated:
  `docs/reports/branch_ab_v2_rnd_constraints_2026-06-18.md`

## 4. Current Conclusion

Current Phase 4 is intentionally conservative:

- Missing matched-negative evidence remains `UNKNOWN_BACKGROUND`.
- V2 classifier outputs `REVIEW_REQUIRED` / `SHADOW_REVIEW_ONLY` for all current materialized Y, H, and 2026-06-15 candidates.
- No V2 output is promoted to final report.
- Current Phase 4 cannot yet be used as FP-reduction evidence because no valid N0 matched-negative background is present.

## 5. Remote Materialization Summary

Remote summary directory:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419/reports/branch_ab_v2_v2_classifier_eval_20260619
```

| dataset | files | rows | truth-overlap rows | V2 result |
|---|---:|---:|---:|---|
| Y1-Y8 locked positives | 8 | 143 | 13 | all `REVIEW_REQUIRED` |
| H1-H16 positives+review | 16 | 307 | 14 | all `REVIEW_REQUIRED` |
| 2026-06-15 shadow-report | 5 | 170 | 0 | all `REVIEW_REQUIRED` |

## 6. Validation

All validation commands ran on `ssh fengxian`.

### RED test

Executable:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/python
```

Command:

```bash
cd /data/project/CNV/PGT-A/refactor_validation_20260419 &&
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_v2_classifier.py -q
```

Result before code change:

```text
FAILED test_phase1_unknown_matched_negative_source_forces_review_not_artifact
assert 'LIKELY_ARTIFACT_SHADOW' == 'REVIEW_REQUIRED'
1 failed, 4 passed in 0.53s
```

### Targeted GREEN test

Executable:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/python
```

Command:

```bash
cd /data/project/CNV/PGT-A/refactor_validation_20260419 &&
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_v2_classifier.py -q
```

Result:

```text
5 passed in 0.53s
```

### Branch A/B/S regression tests

Executable:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/python
```

Command:

```bash
cd /data/project/CNV/PGT-A/refactor_validation_20260419 &&
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_a_candidates.py \
  tests/unit/test_branch_b_evidence_ledger.py \
  tests/unit/test_branch_b_matched_negative.py \
  tests/unit/test_branch_b_v2_classifier.py \
  tests/unit/test_branch_s_shadow.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py -q
```

Result:

```text
23 passed in 0.74s
```

### Snakemake materialization

Executable:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake
```

Command pattern:

```bash
cd /data/project/CNV/PGT-A/refactor_validation_20260419 &&
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake \
  -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile \
  --configfile <predict_config.yaml> \
  --cores 4 \
  --rerun-triggers mtime \
  --forcerun cnv_branch_b_v2_classifier \
  cnv
```

Configs:

```text
config_predict_h_20260608_mask_only.yaml
config_predict_20260615_mask_only.yaml
config_predict_build_ref_v2_mask_only.yaml
```

Results:

- H1-H16: completed 19/19 Snakemake steps.
- 2026-06-15: completed 8/8 Snakemake steps.
- Y1-Y8: completed 11/11 Snakemake steps.

## 7. Known Risks

- Current V2 output is all-review because matched-negative background is unknown.
- H7-H16 remain N1/N2 candidates, not N0 clean negatives.
- Any future N0 or reference expansion must rerun WisecondorX predict, Branch A candidates, Branch B, evaluation, benchmark, and report before conclusions are updated.
- This run did not re-run full final reports because the change is shadow-only and final report impact is explicitly `none_shadow_only`.

## 8. Recommended Next Step

Build a valid negative-bank experiment before treating Phase 4 as an FP-reduction classifier:

1. define locked N0 negatives separately from H7-H16 review candidates;
2. rerun matched-negative percentile in shadow mode;
3. report per-sample review burden and known-positive overlap;
4. only then decide whether any V2 classifier label can reduce FP without increasing FN.
