# Plot Event Manifest Low-Confidence Ablation

Date: 2026-06-22

## Scope

This loop only changes plot/support/report-layer contracts. It does not modify
Branch A, Branch B V2 classifier/filtering, Branch S classification, reference
build, mapping, sex calling, or TP/FN/FP definitions.

## Implemented Contract

- Added per-sample `plot_event_manifest.tsv` as the single source of truth for
  z/CN plot overlays and `plot_event_support.tsv` event coordinates.
- Manifest rows include `event_id`, `candidate_id`, `sample_id`, `chrom`,
  `start`, `end`, `state`, `event_layer`, `plot_visibility`, and
  `plot_support_class`.
- Same-direction events split by centromere/structure gaps can remain one table
  event while z/CN plots draw separate line chunks and leave the gap blank.
- Added plot support classes:
  - `Z_AND_CN_SUPPORTED`
  - `Z_SUPPORTED_CN_WEAK`
  - `Z_SUPPORTED_CN_NOT_SUPPORTED`
  - `Z_AND_CN_NOT_SUPPORTED`
- Autosomal `report_event` rows with `Z_SUPPORTED_CN_NOT_SUPPORTED` are
  demoted only at plot/report visibility level to
  `internal_review_event_candidate` and `review_plot_only`.
- `CN_DIRECTION_WEAK_OR_MIXED` is not demoted, protecting weak/mosaic-sensitive
  events such as H6 chr21, Y6 chr7, and G2 chr8.
- Branch S rows remain separate `branch_s_review` rows and are not mixed into
  autosomal CNV report-table decisions.

## Remote Validation

All validation was run on `fengxian` under:

`/data/project/CNV/PGT-A/refactor_validation_20260419`

Commands:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_plot.py \
  tests/unit/test_cnv_report.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py \
  tests/unit/test_branch_s_shadow.py -q
```

Result:

`69 passed in 4.17s`

Dry-runs for Y, H, G, and 0615 active lowres/gap2m configs passed for:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake \
  -s Snakefile --configfile <config> --cores 1 -n \
  branch_s_review cnv_report
```

Materialized plot/report outputs for:

- Y1-Y8
- H1-H16
- G1-G8
- 2026-06-15 context cohort

## Materialized Result Summary

Manifest files:

- Y: 8/8
- H: 16/16
- G: 8/8
- 0615: 5/5

Low-confidence autosomal report events demoted to review-only plot visibility:

- Y: 1
- H: 1
- G: 0 report-layer events; 2 low-confidence rows were already internal review
- 0615: 3 context-only events

Truth gate:

- Y truth: 10 events, FN=0, hard-suppressed truth=0
- H truth: 10 events, FN=0, hard-suppressed truth=0
- G truth: 10 events, FN=0, hard-suppressed truth=0

0615 remains burden/context only. No TP/FN/FP is computed for 0615.

## Boundary

This ablation hides low-confidence events from final plot overlays only. It does
not delete evidence, hard-filter events, or change report-event classifier
outputs. Demoted events remain in `plot_event_manifest.tsv`,
`plot_event_support.tsv`, and review/audit outputs.

