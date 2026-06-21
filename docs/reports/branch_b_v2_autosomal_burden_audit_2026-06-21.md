# Branch B V2 Autosomal Burden Audit

Date: 2026-06-21

Status: `materialized_result_audit_no_classifier_change`

## Scope

This report audits the remaining autosomal `V2_POSITIVE_SUPPORT_REVIEW` burden
after the Branch B V2 sex-route refinement.

Active contract:

```text
reference_id=h_r0_shadow_ref_20260619
branch_a.merge_gap_bp=2_000_000
branch_a.strong_z=10.0
postprocess.output_dir=wisecondorx/cnv/postprocess_gap2m
final_report_impact=none_shadow_only
```

This audit does not change classifier logic. It is meant to identify which
next refinements are safe and which tempting filters would risk truth loss or
truth demotion.

## Remote Inputs

Remote mirror:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419
```

Audited inputs:

```text
results_y_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m/v2_classifier/*.candidate_classification.tsv
results_h_20260608_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m/v2_classifier/*.candidate_classification.tsv
results_20260615_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m/v2_classifier/*.candidate_classification.tsv
results_y_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m/v2_benchmark/truth_metrics.tsv
results_h_20260608_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m/v2_benchmark/truth_metrics.tsv
```

Rows audited:

```text
event_arm_class=autosome
v2_candidate_class=V2_POSITIVE_SUPPORT_REVIEW
```

## Remaining Autosomal Burden

| cohort | autosomal positive rows | samples with rows | strong | sensitive | z 5-10 | z 10-20 | z 20-50 | z >=50 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Y1-Y8 truth cohort | 78 | 8 | 29 | 49 | 49 | 13 | 8 | 8 |
| H1-H6 truth subset | 33 | 5 | 11 | 22 | 22 | 1 | 5 | 5 |
| H7-H16 context subset | 24 | 4 | 2 | 22 | 22 | 2 | 0 | 0 |
| 2026-06-15 context | 127 | 5 | 72 | 55 | 55 | 43 | 22 | 7 |

Length / CNVpro-like tier:

| cohort | large_ge10mb | large_ge4mb | reportable_ge2mb | review_ge1mb | median length |
|---|---:|---:|---:|---:|---:|
| Y1-Y8 truth cohort | 67 | 3 | 4 | 4 | 39.0 Mb |
| H1-H6 truth subset | 29 | 2 | 0 | 2 | 53.25 Mb |
| H7-H16 context subset | 22 | 1 | 0 | 1 | 39.0 Mb |
| 2026-06-15 context | 103 | 9 | 7 | 8 | 40.5 Mb |

The remaining burden is mostly broad-segment signal, not tiny-bin-only noise.

## Background And Technical Status

For all audited autosomal positive rows in all cohorts:

```text
matched_negative_background_status=UNKNOWN_BACKGROUND
calibration_null_status=NO_NULL_SUPPORT
refmap_status=OK
sample_noise_status=OK
cnvpro_like_evidence_status=SHADOW_EVIDENCE_ONLY
v2_evidence_tier=UNKNOWN_BACKGROUND_POSITIVE_SUPPORT
```

This means the current V2 classifier is still operating mainly in an unknown
background mode for autosomes. It cannot use matched-negative percentile or
calibration null evidence to safely call benign/artifact.

## Direction-Support Audit

The audit tested a possible direction-support split using Branch B-side
evidence only:

```text
B_DIRECTION_SUPPORTED:
  same_direction_fraction >= 0.5
  OR corrected/raw amplitude supports the same event direction with abs >= 2

A_ONLY_NO_B_DIRECTION_SUPPORT:
  Branch A z/support label exists, but Branch B-side direction support is weak
```

Impact:

| cohort | autosomal positive rows | B-direction supported | A-only / weak B-direction |
|---|---:|---:|---:|
| Y1-Y8 truth cohort | 78 | 58 | 20 |
| H1-H16 all | 57 | 35 | 22 |
| H7-H16 context subset | 24 | 15 | 9 |
| 2026-06-15 context | 127 | 89 | 38 |

This split is useful as review evidence, but it is not safe as a hard filter or
as a universal positive-support downgrade.

## Truth Safety Check

Several locked autosomal truth top candidates would be labeled
`A_ONLY_NO_B_DIRECTION_SUPPORT` by the direction-support audit:

| sample | truth event | top candidate | A z | same_direction_fraction | corrected amplitude | audit label |
|---|---|---|---:|---:|---:|---|
| Y2 | chr14 gain | Y2.A0001 | 94.785 | 0.152 | -0.289 | A-only |
| Y4 | chr13 gain | Y4.A0001 | 123.904 | 0.328 | 1.234 | A-only |
| H2 | chr6 gain | H2.A0005 | 54.714 | 0.424 | 0.461 | A-only |
| H3 | chr13 gain | H3.A0002 | 27.977 | 0.304 | 1.182 | A-only |
| H4 | chr15 gain | H4.A0001 | 111.546 | 0.382 | 0.276 | A-only |

H6 chr21 gain remains direction-supported:

```text
sample=H6
truth=chr21 gain
top_candidate=H6.A0003
A z=7.1135
same_direction_fraction=0.917
corrected_amplitude=0.404
audit label=B_DIRECTION_SUPPORTED
```

Conclusion:

```text
Do not convert weak Branch B direction support into hard suppression.
Do not remove positive review status solely because B-direction support is weak.
```

The most likely interpretation is that Branch B-side corrected amplitude or
direction summary can be attenuated or inconsistent even for broad true events,
especially under unknown background and no calibration-null support. This is
exactly why Branch B V2 must remain review/evidence mode until a better
background model exists.

## Context Cohort Patterns

H7-H16 context autosomal burden:

- 24 autosomal positive rows.
- Concentrated in four samples:
  - H13: 10
  - H14: 10
  - H8: 3
  - H7: 1
- Mostly sensitive rather than strong:
  - strong=2
  - sensitive=22
- Mostly broad:
  - large_ge10mb=22
- z-score range:
  - min=5.195
  - median=6.141
  - max=15.879

2026-06-15 context autosomal burden:

- 127 autosomal positive rows across all five samples.
- Sample burden:
  - JZ26125845-60-60: 34
  - JZ26125844-59-59: 27
  - JZ26125846-61-61: 27
  - JZ26125843-56-56: 20
  - JZ26125847-62-62: 19
- Many rows are strong:
  - strong=72
  - sensitive=55
- Mostly broad:
  - large_ge10mb=103
- z-score range:
  - min=5.003
  - median=11.345
  - max=122.763

Because 2026-06-15 has no locked truth, these rows cannot be converted into
TP/FN/FP metrics and cannot prove positive or negative sample status.

## Current Interpretation

The remaining Branch B V2 autosomal burden is not explained by the existing
region-risk or sample-noise fields:

- `refmap_status=OK`
- `sample_noise_status=OK`
- `hard_region_fraction` median is modest in all audited cohorts
- most rows are broad segments

The dominant blocker is background:

```text
UNKNOWN_BACKGROUND + NO_NULL_SUPPORT
```

Without a reliable matched-negative / calibration-null model, Branch B V2 can
separate review evidence classes, but it should not issue final benign/artifact
suppression for autosomal candidates.

## Safe Next Refinement

Recommended next implementation should be review-label only:

1. Keep `V2_POSITIVE_SUPPORT_REVIEW` for truth-safe preservation.
2. Add or report a non-final review sublabel for autosomal rows:
   - `B_DIRECTION_SUPPORTED`
   - `A_ONLY_WEAK_B_DIRECTION`
3. Do not hard-suppress or final-demote `A_ONLY_WEAK_B_DIRECTION`, because it
   includes locked truth top candidates.
4. Preserve H6 chr21 as a non-regression sentinel.
5. Continue using H7-H16 and 2026-06-15 only as burden/context unless locked
   truth labels are added.

## Do Not Promote

This audit does not support:

- final Branch B V2 promotion;
- a direction-support hard filter;
- matched-negative or calibration-null promotion;
- 2026-06-15 TP/FN/FP conclusions;
- final SCA or Branch S promotion;
- production promotion of `h_r0_shadow_ref_20260619`;
- default promotion of `merge_gap_bp=2_000_000`.
