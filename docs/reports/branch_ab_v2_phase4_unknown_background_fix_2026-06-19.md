# Branch B V2 Phase 4 Unknown Background Fix

Date: 2026-06-19

## Scope

This report documents a Branch B V2 shadow-classifier contract fix.

It does not change:

- WisecondorX predict;
- Branch A candidate generation;
- legacy Branch B final report decisions;
- Branch S;
- `cnv_report` schema or report promotion.

## Root Cause

`cnv_branch_b_v2_classifier` can consume two input shapes:

| input source | matched-negative field shape |
|---|---|
| Phase 3 matched-negative output | `matched_negative_background_status=UNKNOWN_BACKGROUND` |
| Phase 1 evidence ledger fallback | `matched_negative_source=UNKNOWN_BACKGROUND`, no `matched_negative_background_status` column |

The previous classifier only checked `matched_negative_background_status`.
Therefore H1-H16 and the 2026-06-15 shadow-report set, which currently use the Phase 1 fallback ledger, could classify legacy artifacts as `LIKELY_ARTIFACT_SHADOW` even though their matched-negative background was explicitly unknown.

That violated the project constraint:

```text
Missing evidence must stay UNKNOWN / NO_CALL, not clean or artifact support.
```

## Code Change

Updated:

```text
pgta/predict/branch_b/v2_classifier.py
tests/unit/test_branch_b_v2_classifier.py
```

New classifier behavior:

```text
matched_negative_background_status == UNKNOWN_BACKGROUND
OR matched_negative_source == UNKNOWN_BACKGROUND
  -> v2_candidate_class = REVIEW_REQUIRED
  -> v2_classifier_action = SHADOW_REVIEW_ONLY
  -> v2_classifier_reason = unknown_matched_negative_background
```

This remains shadow-only:

```text
v2_final_report_impact = none_shadow_only
```

## TDD Evidence

Remote executable:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/python
```

RED command:

```bash
cd /data/project/CNV/PGT-A/refactor_validation_20260419 &&
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_v2_classifier.py -q
```

RED result before code change:

```text
FAILED test_phase1_unknown_matched_negative_source_forces_review_not_artifact
assert 'LIKELY_ARTIFACT_SHADOW' == 'REVIEW_REQUIRED'
1 failed, 4 passed in 0.53s
```

GREEN result after code change:

```text
5 passed in 0.53s
```

Regression command:

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

Regression result:

```text
23 passed in 0.74s
```

## Remote Materialization

Remote executable:

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

Configs refreshed:

```text
config_predict_h_20260608_mask_only.yaml
config_predict_20260615_mask_only.yaml
config_predict_build_ref_v2_mask_only.yaml
```

Dry-run note:

- Without `--rerun-triggers mtime`, Snakemake provenance would rerun many upstream jobs.
- With `--rerun-triggers mtime`, the H dry-run narrowed to 16 `cnv_branch_b_v2_classifier` jobs plus runtime tracking.

Execution result:

- H1-H16: 16 V2 classifier files refreshed.
- 2026-06-15: 5 V2 classifier files refreshed.
- Y1-Y8: 8 V2 classifier files refreshed.

## Post-Fix Summary

Remote summary directory:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419/reports/branch_ab_v2_v2_classifier_eval_20260619
```

| dataset | files | rows | truth-overlap rows | V2 class result |
|---|---:|---:|---:|---|
| Y1-Y8 locked positives | 8 | 143 | 13 | all `REVIEW_REQUIRED` |
| H1-H16 positives+review | 16 | 307 | 14 | all `REVIEW_REQUIRED` |
| 2026-06-15 shadow-report | 5 | 170 | 0 | all `REVIEW_REQUIRED` |

Class/reason counts after the fix:

| dataset | `v2_candidate_class` | `v2_classifier_reason` |
|---|---|---|
| Y1-Y8 | `REVIEW_REQUIRED=143` | `unknown_matched_negative_background=143` |
| H1-H16 | `REVIEW_REQUIRED=307` | `unknown_matched_negative_background=307` |
| 2026-06-15 | `REVIEW_REQUIRED=170` | `unknown_matched_negative_background=170` |

Matched-negative field status:

| dataset | status |
|---|---|
| Y1-Y8 | `matched_negative_source=UNKNOWN_BACKGROUND`, `matched_negative_background_status=UNKNOWN_BACKGROUND` |
| H1-H16 | `matched_negative_source=UNKNOWN_BACKGROUND`, `matched_negative_background_status` absent |
| 2026-06-15 | `matched_negative_source=UNKNOWN_BACKGROUND`, `matched_negative_background_status` absent |

## Interpretation

The current Phase 4 shadow classifier is deliberately conservative because no usable N0 matched-negative background exists in these materialized runs.

The fix prevents a false sense of precision:

- unknown background no longer supports `LIKELY_ARTIFACT_SHADOW`;
- truth-overlap candidates in H1-H6 remain review-only rather than artifact-shadow;
- 2026-06-15 candidates remain shadow review-only and are not validation proof;
- final current-workflow reports are unchanged.

## Next Step

To make Phase 4 informative rather than all-review, the next required evidence is a valid negative-bank / matched-negative input:

1. define N0 only from locked clean negatives;
2. keep H7-H16 as N1/N2 unless manually reviewed and ablated;
3. rerun Branch A, evidence ledger, matched-negative, V2 classifier, evaluation, benchmark, and report for any new reference or negative-bank variant;
4. only then evaluate whether Branch B V2 can reduce FP without increasing FN.
