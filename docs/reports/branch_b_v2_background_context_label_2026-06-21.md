# Branch B V2 Background Context Review Label

Date: 2026-06-21

Status: `materialized_review_label_only`

## 1. Purpose

This loop makes the unresolved Branch B V2 background condition explicit in
the classifier output:

```text
UNKNOWN_BACKGROUND + NO_NULL_SUPPORT
```

The change is intentionally review-only. It does not promote Branch B V2, does
not suppress candidates, and does not change final-report impact.

## 2. Scope

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

The classifier now emits two additional columns:

```text
v2_background_context_label
v2_background_context_reason
```

The sample-level summary now emits:

```text
background_context_label_counts
```

## 3. Label Contract

Current materialized labels:

| label | interpretation |
|---|---|
| `UNKNOWN_BACKGROUND_NO_NULL_SUPPORT` | No matched-negative background is available and calibration-null support is absent. This is unresolved review context, not benign evidence. |
| `SHADOW_BACKGROUND_NO_NULL_SUPPORT` | Shadow background context exists, but calibration-null support is absent. This is context-only and not a release/promotion gate. |

Other labels are implemented for future inputs, including matched-negative
informative background, limited-null support, and calibration-null-only
contexts.

## 4. What Did Not Change

This loop does not alter:

- WisecondorX predict;
- Branch A candidate discovery;
- `v2_candidate_class`;
- `v2_classifier_action`;
- `v2_evidence_tier`;
- `v2_direction_support_label`;
- `v2_final_report_impact`;
- Branch S status;
- report-release status.

All outputs remain:

```text
v2_final_report_impact=none_shadow_only
```

## 5. Remote Verification

Unit test command:

```bash
cd /data/project/CNV/PGT-A/refactor_validation_20260419
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_v2_classifier.py \
  tests/unit/test_branch_b_v2_benchmark.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py -q
```

Result:

```text
35 passed in 0.91s
```

Snakemake dry-run command pattern:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake \
  -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile \
  --configfile <config> \
  --cores 1 -n branch_b_v2_benchmark
```

Configs:

- `config_predict_y_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`
- `config_predict_h_20260608_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`
- `config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`

Result after materialization:

```text
Nothing to be done (all requested files are present and up to date).
```

Forced materialization used:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake \
  -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile \
  --configfile <config> \
  --cores 8 \
  --forcerun cnv_branch_b_v2_classifier cnv_branch_b_v2_benchmark \
  branch_b_v2_benchmark
```

Rules rerun:

- `cnv_branch_b_v2_classifier`
- `cnv_branch_b_v2_benchmark`
- `branch_b_v2_benchmark`
- dependent `collect_runtime_tracking`

All three materializations completed.

## 6. Materialized Counts

| cohort | rows | UNKNOWN background no null | SHADOW background no null |
|---|---:|---:|---:|
| Y1-Y8 truth cohort | 97 | 86 | 11 |
| H1-H16 | 105 | 69 | 36 |
| 2026-06-15 context | 165 | 155 | 10 |

Class counts remain unchanged from the prior direction-support materialization:

| cohort | positive support review | sex chromosome review | no-call contract risk |
|---|---:|---:|---:|
| Y1-Y8 | 78 | 13 | 6 |
| H1-H16 | 57 | 42 | 6 |
| 2026-06-15 | 127 | 14 | 24 |

Direction label counts also remain unchanged:

| cohort | B direction supported | A-only weak B direction | B direction conflict |
|---|---:|---:|---:|
| Y1-Y8 | 66 | 20 | 11 |
| H1-H16 | 68 | 26 | 11 |
| 2026-06-15 | 97 | 40 | 28 |

## 7. Truth Preservation

| cohort | truth events | truth preserved | FN | hard-suppressed truth |
|---|---:|---:|---:|---:|
| Y1-Y8 | 10 | 10 | 0 | 0 |
| H1-H16 | 10 | 10 | 0 | 0 |
| 2026-06-15 | 0 | 0 | 0 | 0 |

H6 chr21 remains preserved:

```text
candidate_id=H6.A0003
v2_candidate_class=V2_POSITIVE_SUPPORT_REVIEW
v2_background_context_label=UNKNOWN_BACKGROUND_NO_NULL_SUPPORT
v2_direction_support_label=B_DIRECTION_SUPPORTED
v2_final_report_impact=none_shadow_only
a_abs_zscore=7.113507302991461
```

## 8. Interpretation

The new background context label is useful because it prevents ambiguous
`UNKNOWN_BACKGROUND` rows from looking like evaluated negatives. It states that
the classifier has no matched-negative empirical background and no calibration
null support.

This is not FP reduction yet. It is a safer evidence contract that makes the
remaining review burden measurable and reportable without pretending that
unknown background is benign.

## 9. Next Gate

The next Branch B V2 step can use this explicit background context to design
truth-safe burden stratification, but it must still obey:

- no hard suppression of truth-overlap candidates;
- no use of legacy/current-code Branch B decision fields;
- no use of 2026-06-15 as TP/FN/FP evidence;
- no promotion of Branch B V2 or Branch S from this label alone.
