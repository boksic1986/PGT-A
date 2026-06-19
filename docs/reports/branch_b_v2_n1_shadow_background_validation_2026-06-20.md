# Branch B V2 N1 Shadow Background Validation

Date: 2026-06-20

## Scope

This report documents the first Branch B V2 matched-negative background update after PR #19.

The change does not promote Branch B V2 into final report decisions.

Current final-report impact remains:

```text
v2_final_report_impact = none_shadow_only
```

## Design

N0 remains the only formal empirical-null label.

N1 can now be used as an explicitly configured shadow background:

```yaml
negative_bank:
  shadow_background_labels: ["N1"]
```

When only N1 background is available, matched-negative output is:

```text
matched_negative_background_status = SHADOW_BACKGROUND
matched_negative_action = SHADOW_CONTEXT_ONLY
```

This is not equivalent to:

```text
matched_negative_background_status = OK
matched_negative_action = BACKGROUND_SUPPORTED
```

The purpose is to provide context from H7-H16 R0/N1 samples without pretending they are locked clean negatives.

## Files Changed

- `pgta/predict/branch_b/matched_negative.py`
- `pgta/predict/branch_b/v2_classifier.py`
- `rules/predict_layout.smk`
- `rules/predict_workflow.smk`
- `config_predict_y_h_r0_shadow_20260619.yaml`
- `config_predict_h_20260608_h_r0_shadow_20260619.yaml`
- `config_predict_20260615_h_r0_shadow_20260619.yaml`
- `tests/unit/test_branch_b_matched_negative.py`
- `tests/unit/test_branch_b_v2_classifier.py`
- `tests/unit/test_branch_ab_phase12_workflow_contract.py`

## TDD Evidence

RED command:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 /biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest tests/unit/test_branch_b_matched_negative.py tests/unit/test_branch_b_v2_classifier.py -q'
```

Expected RED failures:

- `build_matched_negative_percentiles()` did not accept `shadow_background_labels`.
- V2 classifier treated `SHADOW_BACKGROUND` as no matched background.

GREEN command:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 /biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest tests/unit/test_branch_b_matched_negative.py tests/unit/test_branch_b_v2_classifier.py tests/unit/test_negative_bank.py tests/unit/test_branch_ab_phase12_workflow_contract.py -q'
```

GREEN result:

```text
33 passed in 0.69s
```

Additional edge case covered:

- if `matched_negative_background_status` exists but is blank / pandas NA, V2 falls back to `matched_negative_source=UNKNOWN_BACKGROUND`;
- if `matched_negative_background_status=SHADOW_BACKGROUND`, the old source fallback cannot overwrite it.

## Workflow Validation

Dry-run command:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && for cfg in config_predict_y_h_r0_shadow_20260619.yaml config_predict_h_20260608_h_r0_shadow_20260619.yaml config_predict_20260615_h_r0_shadow_20260619.yaml; do /biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile "$cfg" --cores 1 --rerun-triggers mtime --forcerun cnv_branch_b_matched_negative cnv_branch_b_v2_classifier -n cnv; done'
```

Dry-run result:

- Y H-R0 shadow: 8 matched-negative + 8 V2 classifier jobs.
- H 2026-06-08 H-R0 shadow: 16 matched-negative + 16 V2 classifier jobs.
- 2026-06-15 H-R0 shadow: 5 matched-negative + 5 V2 classifier jobs.
- No `map_reads` jobs.
- No `wisecondorx_convert` jobs.
- No `wisecondorx_predict` jobs.

Actual refresh command:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && for cfg in config_predict_y_h_r0_shadow_20260619.yaml config_predict_h_20260608_h_r0_shadow_20260619.yaml config_predict_20260615_h_r0_shadow_20260619.yaml; do /biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile "$cfg" --cores 16 --rerun-triggers mtime --forcerun cnv_branch_b_matched_negative cnv_branch_b_v2_classifier cnv; done'
```

Classifier-only refresh after status fallback fix:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && for cfg in config_predict_y_h_r0_shadow_20260619.yaml config_predict_h_20260608_h_r0_shadow_20260619.yaml config_predict_20260615_h_r0_shadow_20260619.yaml; do /biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile "$cfg" --cores 16 --rerun-triggers mtime --forcerun cnv_branch_b_v2_classifier cnv; done'
```

Latest classifier-refresh log audit:

| config | V2 classifier jobs | map reads | convert | predict |
|---|---:|---:|---:|---:|
| Y H-R0 shadow | 8 | 0 | 0 | 0 |
| H 2026-06-08 H-R0 shadow | 16 | 0 | 0 | 0 |
| 2026-06-15 H-R0 shadow | 5 | 0 | 0 | 0 |

## Shadow Background Results

| cohort | candidates | shadow background | unknown background |
|---|---:|---:|---:|
| Y1-Y8 | 131 | 38 | 93 |
| H1-H16 | 221 | 145 | 76 |
| 2026-06-15 five samples | 201 | 32 | 169 |

V2 evidence tiers:

| cohort | shadow compatible | shadow outlier positive | unknown positive | unknown direction conflict |
|---|---:|---:|---:|---:|
| Y1-Y8 | 25 | 13 | 84 | 9 |
| H1-H16 | 107 | 38 | 66 | 10 |
| 2026-06-15 five samples | 16 | 16 | 141 | 28 |

All rows still have:

```text
v2_final_report_impact = none_shadow_only
```

## 2026-06-15 Current Report

Source:

`/data/project/CNV/PGT-A/refactor_validation_20260419/results_20260615_h_r0_shadow_ref_20260619/wisecondorx/cnv/report/cnv_summary.tsv`

| sample | QC | sex | Branch A candidates | Branch B kept | reportable | top event |
|---|---|---|---:|---:|---:|---|
| JZ26125843-56-56 | PASS | XX | 28 | 0 | 0 | none |
| JZ26125844-59-59 | PASS | XX | 38 | 3 | 3 | chr4:52500001-186000000 gain |
| JZ26125845-60-60 | PASS | XY | 57 | 26 | 26 | chr14:60750001-106500000 gain |
| JZ26125846-61-61 | PASS | XY | 44 | 3 | 3 | chr12:37500001-124500000 gain |
| JZ26125847-62-62 | PASS | XY | 34 | 1 | 1 | chr4:121500001-186000000 loss |

The current report output is unchanged by V2 shadow background because V2 is not connected to final report promotion.

## Interpretation

This update solves the immediate V2 background problem without violating the N0 constraint:

- H9/H10/H11/H12/H15/H16 can contribute N1 shadow context.
- N1 shadow context can mark candidates as background-compatible or outlier-positive.
- N1 shadow context cannot become a hard artifact/pass decision.
- Missing background still remains unknown.

For the 2026-06-15 set, N1 shadow context is informative but not sufficient to issue final negative/positive calls:

- 56 still has no Branch B reportable event.
- 59, 61, and 62 remain low-count reportable/review events under the current report logic.
- 60 remains high-burden and needs separate review.

## Remaining Gate

Before Branch B V2 can affect final reports, the project still needs either:

1. locked held-out N0 negatives, or
2. a documented cross-fit design where each query sample is evaluated against a reference/background that did not include it.

N1 shadow background is a review aid, not a promotion gate.
