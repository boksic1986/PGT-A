# Handoff: Sex-Aware CN/Z Plot And Branch S Overlay

Date: 2026-06-22 15:35

Status: active handoff

## What Changed

- Combined genome plots now consume gender TSV and Branch S summary/scores/evidence.
- Branch S summary JSON is supported directly by the plot command.
- Branch S `sca_report_review_event` intervals are drawn on chrX/chrY in both
  z and CN plots.
- CN scatter is sex-aware for chrX/chrY.
- chrY ratio-derived CN is marked not interpretable when the chromosome-level
  reference denominator is too low, preventing huge false dup points.

No Branch A, autosomal Branch B V2 filtering, Branch S classifier, reference,
mapping, or report-event classification was changed.

## Files Modified

- `pgta/predict/branch_b/plot.py`
- `rules/predict_workflow.smk`
- `tests/unit/test_branch_b_plot.py`
- `tests/unit/test_branch_ab_phase12_workflow_contract.py`
- `docs/reports/sex_aware_cn_z_plot_branch_s_overlay_2026-06-22.md`
- `docs/CURRENT_CONTEXT_INDEX.md`
- `CURRENT_STATE.md`
- `PLANS.md`

## Remote Validation

Remote path:

`/data/project/CNV/PGT-A/refactor_validation_20260419`

Pytest:

```bash
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_plot.py \
  tests/unit/test_branch_s_shadow.py \
  tests/unit/test_cnv_report.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py -q
```

Result:

`56 passed in 1.86s`

Dry-runs:

- `config_predict_g1_g8_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`
- `config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`

Targets:

`cnv_branch_s_shadow cnv_branch_ab_plot cnv_report`

Result:

- Both parsed successfully.
- No mapping or reference rebuild jobs were requested.

Materialized:

- G1-G8: `20 of 20 steps (100%) done`.
- 2026-06-15: `14 of 14 steps (100%) done`.

## Key Results

G1-G8:

- G3 and G5 X-loss Branch S events are visible on chrX in both z and CN plots.
- G3/G5 CN TSV marks chrX non-PAR bins as `event_layer=branch_s_review`.
- G7 XX chrY is `sex_aware_absent_expected` and neutral.
- G8 XY chrX uses expected CN = 1; it is no longer interpreted against a
  diploid chrX baseline.
- G8 chrY is `sex_chrom_ref_ratio_not_interpretable` for all chrY bins because
  chrY `chrom_ref_cpm_median=0`; huge chrY dup points are removed.

2026-06-15:

- 5/5 samples have refreshed z SVG, CN SVG, and CN TSV outputs.
- XY sample chrY bins are marked `sex_chrom_ref_ratio_not_interpretable`.
- 0615 remains no-truth burden/context only.

## Local Synced Outputs

- `D:\Pipeline\PGT-A\reports\g1_g8_cnv_plots\`
- `D:\Pipeline\PGT-A\reports\0615_cnv_plots\`

## Current Boundaries

- z plot y-axis is still bin-level `calibrated_z`.
- Branch A `a_zscore` and Branch S `state_score` are event-level support labels,
  not bin-level z values.
- Sex-aware CN is visualization-only and must not become a hard filter without
  a separate truth-safe benchmark.
- Branch S remains review/development-only, not final SCA promotion.

## Next Step

Review G1-G8 plots first because they have locked truth. Use:

- z plot for calibrated-z visual context;
- CN plot for sex-aware copy-number visual context;
- `plot_bins_cn.tsv` for bin-level CN interpretation status;
- Branch S summary/scores/evidence for SCA review context.

Do not use 0615 to derive thresholds. It has no locked truth and remains
burden/context only.
