# Handoff: Plot Event Manifest And Low-Confidence Ablation

Date: 2026-06-22 21:45

## Current Context

Use this handoff after:

`docs/handoff/2026-06-22_1938_sex_specific_ref_cnv_plot_handoff.md`

Active reference remains:

`h_r0_shadow_ref_20260619`

This handoff covers only plot/support/report-layer visibility. It does not
change Branch A, Branch B V2 classifier/filtering, Branch S classification,
reference, mapping, sex calling, or TP/FN/FP definitions.

## Completed

- Added `plot_event_manifest.tsv` output to `cnv_branch_ab_plot`.
- Added manifest dependency to `cnv_report_summary`.
- z plot, CN plot, and `plot_event_support.tsv` now use a common manifest for
  event coordinates and visibility.
- Same event can be represented as one manifest/report row while plot lines are
  split around centromere/structure gaps.
- Added `plot_support_class` and visibility ablation:
  - `Z_SUPPORTED_CN_NOT_SUPPORTED` autosomal report events become
    `internal_review_event_candidate` and `review_plot_only`.
  - `CN_DIRECTION_WEAK_OR_MIXED` remains visible and is not demoted.
- Support TSV now carries `event_id`, `candidate_id`, `plot_visibility`,
  `plot_layer_class`, `plot_support_class`, and
  `merged_source_event_ids`.

## Files Modified

- `pgta/predict/branch_b/plot.py`
- `rules/predict_layout.smk`
- `rules/predict_workflow.smk`
- `tests/unit/test_branch_b_plot.py`
- `docs/reports/plot_event_manifest_low_confidence_ablation_2026-06-22.md`
- `docs/handoff/2026-06-22_2145_plot_event_manifest_low_confidence_handoff.md`
- `docs/CURRENT_CONTEXT_INDEX.md`
- `CURRENT_STATE.md`
- `PLANS.md`

## Remote Validation

Remote path:

`/data/project/CNV/PGT-A/refactor_validation_20260419`

Pytest:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_plot.py \
  tests/unit/test_cnv_report.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py \
  tests/unit/test_branch_s_shadow.py -q
```

Result:

`69 passed in 4.17s`

Snakemake dry-run:

Y/H/G/0615 active lowres/gap2m configs all parsed for:

`branch_s_review cnv_report`

Materialized outputs:

- Y1-Y8
- H1-H16
- G1-G8
- 0615 context cohort

## Acceptance

- Y truth: 10 events, FN=0, hard-suppressed truth=0.
- H truth: 10 events, FN=0, hard-suppressed truth=0.
- G truth: 10 events, FN=0, hard-suppressed truth=0.
- H6 chr21 remains protected because `CN_DIRECTION_WEAK_OR_MIXED` is not
  demoted.
- G2 truth remains preserved.
- 0615 remains context-only and is not used for TP/FN/FP.

## Current Result Counts

Low-confidence autosomal report events demoted to review-only plot visibility:

- Y: 1
- H: 1
- G: 0 report-layer events; 2 low-confidence rows were already internal review
- 0615: 3 context-only events

## Next Recommended Step

Review the final plot overlays and `plot_event_manifest.tsv` for the demoted
events first. If the visibility ablation looks reasonable, the next method
step is a formal Branch B V2 report-layer ablation table comparing:

- report events before plot visibility ablation,
- final-report-plot events after ablation,
- internal review candidates,
- truth overlap protection status.

Do not promote this visibility ablation to production filtering without that
table and locked truth review.

