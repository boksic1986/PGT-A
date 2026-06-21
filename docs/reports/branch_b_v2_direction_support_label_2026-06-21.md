---
status: development_evidence
decision_use: branch_b_v2_review_label_contract
date: 2026-06-21
active_reference_id: h_r0_shadow_ref_20260619
branch_a_overlay: merge_gap_bp_2_000_000
---

# Branch B V2 Direction-Support Review Label

This report records the first implementation of Branch B-side direction support
as a non-final review label.

It does not change WisecondorX predict, Branch A candidate discovery, Branch B
V2 candidate class, hard-suppression behavior, final report impact, or Branch S
status.

## Contract

The Branch B V2 classifier now emits two additional review fields:

```text
v2_direction_support_label
v2_direction_support_reason
```

Allowed labels in this implementation:

| label | meaning | final decision impact |
|---|---|---|
| `B_DIRECTION_SUPPORTED` | Branch B-side direction support is present by `same_direction_fraction >= 0.50` or same-direction raw/corrected amplitude with absolute value at least 2. | none |
| `A_ONLY_WEAK_B_DIRECTION` | Branch A positive support is present, but Branch B-side direction support is weak or absent. | none |
| `B_DIRECTION_CONFLICT` | raw/corrected amplitude has opposite direction with absolute value at least 2. | none by label itself |
| `NO_POSITIVE_SUPPORT` | no positive-support candidate context is present for direction comparison. | none |

The label is review evidence only:

- it must not hard-suppress a candidate;
- it must not mark a candidate benign or artifact;
- it must not demote `V2_POSITIVE_SUPPORT_REVIEW` out of positive-support review;
- it must keep `v2_final_report_impact=none_shadow_only`.

## Why It Is Label-Only

The previous autosomal burden audit found that weak Branch B-side direction
support is not safe as a hard filter. Several locked truth top candidates would
be labeled `A_ONLY_WEAK_B_DIRECTION` if this split were used as a decision rule:

- Y2 chr14 gain;
- Y4 chr13 gain;
- H2 chr6 gain;
- H3 chr13 gain;
- H4 chr15 gain.

Therefore the field is implemented only as review context.

## Remote Materialized Results

Remote executable:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake
```

Remote code tree:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419
```

Forced materialization reran:

```text
cnv_branch_b_v2_classifier
cnv_branch_b_v2_benchmark
branch_b_v2_benchmark
```

for:

- `config_predict_y_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`;
- `config_predict_h_20260608_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`;
- `config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`.

## Direction Label Counts

| cohort | rows | B direction supported | A-only weak B direction | B direction conflict |
|---|---:|---:|---:|---:|
| Y1-Y8 truth cohort | 97 | 66 | 20 | 11 |
| H1-H16 | 105 | 68 | 26 | 11 |
| 2026-06-15 context | 165 | 97 | 40 | 28 |

Autosomal `V2_POSITIVE_SUPPORT_REVIEW` rows:

| cohort | autosomal positive rows | B direction supported | A-only weak B direction |
|---|---:|---:|---:|
| Y1-Y8 truth cohort | 78 | 58 | 20 |
| H1-H16 | 57 | 35 | 22 |
| 2026-06-15 context | 127 | 89 | 38 |

H6 chr21 remains:

```text
candidate_id=H6.A0003
v2_candidate_class=V2_POSITIVE_SUPPORT_REVIEW
v2_direction_support_label=B_DIRECTION_SUPPORTED
v2_final_report_impact=none_shadow_only
top_a_abs_zscore=7.113507302991461
```

## Truth Preservation

| cohort | truth events | truth preserved | FN | hard-suppressed truth |
|---|---:|---:|---:|---:|
| Y1-Y8 | 10 | 10 | 0 | 0 |
| H1-H16 | 10 | 10 | 0 | 0 |
| 2026-06-15 | 0 | 0 | 0 | 0 |

The 2026-06-15 cohort has no locked truth and remains burden/context only.

## Verification

Remote pytest:

```text
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_v2_classifier.py \
  tests/unit/test_branch_b_v2_benchmark.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py -q
```

Result:

```text
32 passed in 0.71s
```

TDD red check before implementation:

```text
3 failed, 11 passed
```

The failures were the expected missing
`v2_direction_support_label` / `direction_support_label_counts` fields.

An additional conflict-priority red check failed as expected before the final
priority fix:

```text
1 failed
```

That test confirmed that `B_DIRECTION_CONFLICT` must take priority over
fraction-based support when raw/corrected amplitude is strongly opposite in
direction.

## Current Conclusion

This loop adds review-visible evidence but does not solve FP/review burden by
itself. The main unresolved condition remains:

```text
UNKNOWN_BACKGROUND + NO_NULL_SUPPORT
```

Next Branch B work should focus on reducing positive-support review burden with
truth-safe background/evidence contracts, not by turning direction support into
a hard filter.
