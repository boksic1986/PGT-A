---
status: development_review_only
decision_use: 20260615_high_confidence_review_candidates
cohort: 2026-06-15
reference_id: h_r0_shadow_ref_20260619
generated_from: remote_read_only_current_results
---

# 0615 High-Confidence Report Candidate Review

Date: 2026-06-22

## Scope

This review uses the current remote materialized 2026-06-15 outputs only. It
does not change workflow code, thresholds, reference, Branch A, Branch B V2, or
Branch S.

Remote path:

`/data/project/CNV/PGT-A/refactor_validation_20260419/results_20260615_h_r0_shadow_ref_20260619`

Primary inputs:

- `wisecondorx/cnv/postprocess_gap2m/v2_benchmark/report_events.tsv`
- `wisecondorx/cnv/plots/{sample}.plot_bins.tsv`
- `wisecondorx/cnv/postprocess_gap2m/v2_classifier/*.candidate_classification.tsv`
- `wisecondorx/cnv/report/cnv_summary.tsv`
- `wisecondorx/cnv/postprocess_gap2m/branch_s/*.summary.json`

2026-06-15 has no locked truth labels. The rows below are therefore
high-confidence review/report candidates, not TP/FP/FN calls and not production
promotion evidence.

## Remote Read-Only Check

Executable:

`/biosoftware/miniconda/envs/snakemake_env/bin/python`

Result:

- `report_events.tsv`: 71 autosomal report rows.
- Plot bin TSV files: 5/5 samples present.
- Plot bin columns: `chrom`, `start`, `end`, `genome_pos`, `z`,
  `report_state`.
- Current plot contract: `z` is the calibrated-z plotting signal. No fallback
  signal is used for this review.

## Conservative Selection Criteria

A high-confidence autosomal candidate must satisfy all of the following:

- `v2_report_visibility = report_strong_event`.
- `v2_b_signal_context_label = B_SIGNAL_SUPPORTED_A_DIRECTION`.
- `v2_lowres_context_label = LOWRES_2MB_3MB_SAME_DIRECTION_SUPPORT`.
- `v2_ref_stability_context = REF_STABILITY_STABLE`.
- `v2_clean_support_label = CLEAN_SUPPORT_AVAILABLE`.
- `v2_gc_rc_context_label != GC_RC_ATTENUATED_SEVERE`.
- It must not overlap the 5/5 same-direction batch-shared region
  `chr4:67.50-101.25Mb gain`.
- It must not overlap a 4/5 long batch-shared region unless bin-level
  calibrated-z shows clear local same-direction support.
- For gain, event-bin median z and mean z must both be positive.
- For loss, event-bin median z and mean z must both be negative.
- Rows with median z close to zero and only sparse extreme bins are not
  considered high-confidence.

The mechanical field filter returns 15 rows. Five of those rows have only weak
bin-level trends with median/mean z close to zero. The conservative current
review table therefore keeps 10 rows.

## Current Batch-Shared Regions

These regions were re-derived from the current 0615 classifier outputs.

| region | state | samples | handling |
|---|---|---:|---|
| chr4:67.50-101.25Mb | gain | 5/5 | exclude from high-confidence |
| chr4:52.50-67.50Mb | gain | 4/5 | batch-shared review |
| chr4:101.25-121.50Mb | gain | 4/5 | batch-shared review |
| chr14:60.75-97.50Mb | gain | 4/5 | batch-shared review |

`chr4:67.50-101.25Mb gain` is the clearest batch-shared candidate and is not
allowed into the high-confidence table.

## Per-Sample Summary

`Branch S context` below is the Branch-S-routed candidate count from
`cnv_summary.tsv`, not a final SCA-positive count.

| sample | original report events | mechanical pass | conservative high-confidence | batch-shared review | Branch S context |
|---|---:|---:|---:|---:|---:|
| JZ26125843-56-56 | 7 | 1 | 0 | 0 | 2 |
| JZ26125844-59-59 | 11 | 0 | 0 | 0 | 1 |
| JZ26125845-60-60 | 23 | 10 | 10 | 4 | 5 |
| JZ26125846-61-61 | 18 | 3 | 0 | 2 | 3 |
| JZ26125847-62-62 | 12 | 1 | 0 | 2 | 3 |

## High-Confidence Autosomal Candidates

Current conservative high-confidence candidates are all in
`JZ26125845-60-60`.

| sample | region | state | A z | median z | mean z |
|---|---|---:|---:|---:|---:|
| JZ26125845-60-60 | chr1:30.75-118.50Mb | gain | 48.15 | 0.394 | 0.298 |
| JZ26125845-60-60 | chr2:6.00-87.00Mb | gain | 18.89 | 0.414 | 0.274 |
| JZ26125845-60-60 | chr7:143.25-159.00Mb | gain | 40.19 | 0.806 | 0.069 |
| JZ26125845-60-60 | chr8:3.00-42.75Mb | loss | -122.76 | -0.372 | -0.442 |
| JZ26125845-60-60 | chr9:5.25-39.00Mb | loss | -47.59 | -0.288 | -0.302 |
| JZ26125845-60-60 | chr10:18.75-36.00Mb | gain | 28.59 | 0.381 | 0.179 |
| JZ26125845-60-60 | chr12:1.50-33.75Mb | loss | -66.78 | -0.290 | -0.433 |
| JZ26125845-60-60 | chr17:1.50-21.00Mb | gain | 19.84 | 0.325 | 0.307 |
| JZ26125845-60-60 | chr18:23.25-72.00Mb | loss | -17.35 | -0.290 | -0.601 |
| JZ26125845-60-60 | chr21:15.00-42.00Mb | loss | -66.08 | -0.573 | -0.670 |

Interpretation:

- These rows have strong Branch A support, Branch B same-direction support,
  2Mb and 3Mb low-resolution same-direction support, stable reference context,
  available clean support, and clear local bin-level calibrated-z direction.
- They are review/report candidates only because the cohort has no locked truth.

## Mechanical-Pass Rows Not Kept As High-Confidence

These rows passed the categorical evidence filter but do not show strong enough
bin-level calibrated-z local trend under the conservative review standard.

| sample | region | state | A z | median z | mean z | reason |
|---|---|---:|---:|---:|---:|---|
| JZ26125843-56-56 | chr6:4.50-44.25Mb | gain | 12.33 | 0.004 | 0.016 | bin-level trend too weak |
| JZ26125846-61-61 | chr4:26.25-48.75Mb | gain | 11.49 | 0.015 | 0.231 | median z close to zero |
| JZ26125846-61-61 | chr8:104.25-123.75Mb | loss | -13.11 | -0.069 | -0.073 | bin-level trend too weak |
| JZ26125846-61-61 | chr9:5.25-39.00Mb | gain | 16.45 | 0.032 | 0.024 | bin-level trend too weak |
| JZ26125847-62-62 | chr4:121.50-186.00Mb | loss | -12.83 | -0.086 | -0.099 | bin-level trend too weak |

## Batch-Shared Review Rows

| sample | region | state | A z | batch region | event overlap | median z | mean z | handling |
|---|---|---:|---:|---|---:|---:|---:|---|
| JZ26125845-60-60 | chr14:60.75-106.50Mb | gain | 45.76 | chr14:60.75-97.50Mb gain 4/5 | 0.803 | 0.933 | 0.915 | batch review, not high-confidence because lowres is not 2Mb+3Mb same-direction |
| JZ26125845-60-60 | chr4:52.50-186.00Mb | gain | 23.84 | chr4:67.50-101.25Mb gain 5/5 | 0.253 | 0.266 | 0.069 | excluded from high-confidence |
| JZ26125846-61-61 | chr4:52.50-101.25Mb | gain | 9.35 | chr4:67.50-101.25Mb gain 5/5 | 0.692 | 0.013 | 1.266 | excluded from high-confidence |
| JZ26125847-62-62 | chr4:67.50-121.50Mb | gain | 14.44 | chr4:67.50-101.25Mb gain 5/5 | 0.625 | 0.027 | 0.117 | excluded from high-confidence |

The same events may also overlap adjacent 4/5 chr4 shared segments. The 5/5
chr4 overlap is the dominant exclusion reason.

## Main Exclusion Reasons

Across the 71 current report events, the main reasons for not calling a row
high-confidence are:

| reason | count |
|---|---:|
| bin median/mean z not same direction | 26 |
| lowres not 2Mb+3Mb same direction | 20 |
| sparse-bin-driven / weak bin trend | 20 |
| not `report_strong_event` | 19 |
| reference context not stable | 11 |
| severe GC/RC attenuation | 6 |
| overlaps 5/5 batch-shared region | 3 |

These reasons are review labels for the 0615 cohort. They are not new workflow
thresholds and must not be used as production filtering rules without locked
truth ablation.

## Current Recommendation

For the 2026-06-15 cohort under the current development contract:

- Prioritize `JZ26125845-60-60` for high-confidence autosomal review/report
  follow-up.
- Do not put autosomal events from 56, 59, 61, or 62 into a high-confidence
  report table yet.
- Keep chr4 shared-gain events as batch/context review, not high-confidence
  candidates.
- Keep Branch S/SCA as a separate development-only context section. It is not
  final SCA evidence.

This result should be used to guide review of the current report output. It is
not a new Branch B V2 rule set and does not promote the shadow reference,
Branch B V2, Branch S, or 0615 report to production-final.
