---
status: development_review_only
decision_use: cn_plot_and_branch_s_current_evidence
created_at: 2026-06-22
active_reference_id: h_r0_shadow_ref_20260619
branch_a_overlay: merge_gap_bp_2000000
---

# Copy Number Centering And Branch S Fix

## Scope

This loop fixes two report-review issues only:

1. CN plots were globally shifted below 2.0 because ratio-derived CN was not
   sample-centered.
2. Branch S could over-call XY chrX gain from sparse high-z bins, while weak
   XX X-loss truth-like events such as G3/G5 stayed internal-only.

No Branch A, autosomal Branch B filter, reference, mapping, or report-event
classification logic was changed.

## CN Plot Change

The CN plot still uses sample/reference ratio-derived CN, not calibrated-z CN.
The active formula is now:

```text
raw_log2R = log2((sample_cpm + 0.001) / (ref_cpm + 0.001))
center_shift = median(raw_log2R over non-gap autosomal bins)
centered_log2R = raw_log2R - center_shift
CN = 2 * 2^centered_log2R
```

`plot_bins_cn.tsv` now carries `raw_log2r`, `log2r`,
`copy_number_centering_log2_shift`, `copy_number`, and
`copy_number_source=normalized_signal_ref_median_log2r_autosome_centered`.

Remote materialized medians after the fix:

| cohort | samples | autosome median CN |
|---|---:|---|
| G1-G8 | 8 | all 2.000 |
| 2026-06-15 | 5 | all 2.000 |

This is visualization/support-summary centering only. It does not feed back into
WisecondorX, Branch A, Branch B V2 filtering, or TP/FP/FN evaluation.

## Branch S Change

Branch S segment-level corroboration no longer uses arithmetic mean. It now uses
median/robust-median support only, because sparse extreme chrX/chrY bins can
inflate the mean and create false SCA review calls.

Weak `X_LOSS` in `sex_call=XX` with Branch A support `10 <= abs(z) < 30` is now
report-review visible. This protects low/mosaic SCA-like X-loss examples while
keeping XY short X-loss or XY X-gain artifacts out of report review.

Current G1-G8 Branch S results after rematerialization:

| sample | sex | SCA state | tier | report layer | reason |
|---|---|---|---|---|---|
| G1 | XX | none_detected | SCA_NO_CALL | sca_no_call | insufficient_sca_evidence |
| G2 | XY | none_detected | SCA_NO_CALL | sca_filtered_or_sex_consistent_event | branch_a_only_uncorroborated_by_nonpar_median |
| G3 | XX | X_LOSS | SCA_REVIEW_WEAK | sca_report_review_event | weak report-visible Branch A support |
| G4 | XX | none_detected | SCA_NO_CALL | sca_no_call | insufficient_sca_evidence |
| G5 | XX | X_LOSS | SCA_REVIEW_WEAK | sca_report_review_event | weak report-visible Branch A support |
| G6 | XY | none_detected | SCA_NO_CALL | sca_filtered_or_sex_consistent_event | branch_a_only_uncorroborated_by_nonpar_median |
| G7 | XX | none_detected | SCA_NO_CALL | sca_no_call | insufficient_sca_evidence |
| G8 | XY | none_detected | SCA_NO_CALL | sca_filtered_or_sex_consistent_event | branch_a_only_uncorroborated_by_nonpar_median |

0615 remains no-truth context only. After this fix, its XY X-gain-like Branch A
signals are no-call/sex-consistent audit rather than strong SCA report review.

## Remote Validation

All commands ran on `fengxian` in:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419
```

Remote tests:

```text
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_plot.py \
  tests/unit/test_branch_s_shadow.py \
  tests/unit/test_cnv_report.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py -q
```

Result:

```text
51 passed in 1.50s
```

Remote dry-run:

- G1-G8, 0615, H1-H16, and Y1-Y8 active lowres/gap2m configs all parsed and
  planned only Branch S, plot, report, and runtime refresh jobs.
- No mapping or reference rebuild was requested.

Remote materialization:

- G1-G8: completed `21 of 21 steps (100%)`.
- 0615: completed `15 of 15 steps (100%)`.

Local synced plot paths:

```text
D:\Pipeline\PGT-A\reports\g1_g8_cnv_plots\
D:\Pipeline\PGT-A\reports\0615_cnv_plots\
```

## Interpretation

The corrected CN plot can now be used to compare report-event horizontal CN
trends with centered per-bin CN scatter. A report event whose z evidence is
strong but centered CN support is weak/mixed should be treated as a review item,
not automatically true or false.

Branch S remains development/review-only. The fix improves G1-G8 visibility and
reduces obvious XY false SCA burden, but it is not final SCA promotion.
