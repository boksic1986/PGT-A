# Sex-Aware CN/Z Plot And Branch S Overlay

Date: 2026-06-22

Status: development visualization contract updated; not a calling or filtering
change.

## Scope

This update fixes the plot-layer contract for combined genome CNV figures. It
does not change Branch A, autosomal Branch B V2 filtering, Branch S
classification, reference build, mapping, or report-event classification.

## Root Causes Fixed

1. Branch S summaries are JSON files, while the plot command previously treated
   all optional inputs as TSV. As a result, G3/G5 `sca_report_review_event`
   X-loss calls existed in Branch S outputs but were not drawn on the combined
   genome plots.
2. Sex chromosome CN scatter used a fixed diploid interpretation. XY chrX bins
   could therefore be colored as deletion even when their CN was near the
   expected haploid baseline.
3. chrY ratio-derived CN is not reliable when the reference median is zero or
   near zero. G8 chrY produced huge CN values because a few non-zero reference
   bins were divided by a chromosome whose median reference denominator was
   zero.

## Implemented Contract

- `pgta/predict/branch_b/plot.py` now accepts and uses:
  - `--gender-tsv`
  - `--branch-s-summary`
  - `--branch-s-scores`
  - `--branch-s-evidence`
- JSON Branch S summary files are read as one-row tables.
- z plots still use bin-level `calibrated_z` for points.
- Branch A/Branch S event-level scores are shown only as event-level overlay
  labels/tooltips; they are not substituted for bin-level z values.
- Branch S `sca_report_review_event` intervals are overlaid on chrX/chrY in the
  combined z and CN plots.
- CN scatter is sex-aware:
  - autosomes: expected CN = 2;
  - XX chrX: expected CN = 2;
  - XY chrX: expected CN = 1;
  - XY chrY: expected CN = 1 only when chrY reference context is interpretable;
  - XX chrY: expected absent, colored neutral.
- chrY CN ratio is marked `sex_chrom_ref_ratio_not_interpretable` when the
  chromosome-level `chrom_ref_cpm_median` is below the interpretation threshold.
  Those bins are neutral in the scatter and do not produce huge dup points.

## Output Schema Additions

`*.plot_bins_cn.tsv` now includes:

- `sex_call`
- `expected_copy_number`
- `copy_number_delta`
- `cn_scatter_state_sex_aware`
- `sex_chrom_region_class`
- `copy_number_interpretation_status`
- `event_layer`
- `chrom_ref_cpm_median`

## Remote Validation

Remote path:

`/data/project/CNV/PGT-A/refactor_validation_20260419`

Pytest command:

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

Dry-run commands:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile config_predict_g1_g8_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml \
  --cores 1 -n --forcerun cnv_branch_s_shadow cnv_branch_ab_plot cnv_report
```

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml \
  --cores 1 -n --forcerun cnv_branch_s_shadow cnv_branch_ab_plot cnv_report
```

Both dry-runs parsed successfully and planned only Branch S, plot, report, and
runtime refresh jobs. They did not request mapping or reference rebuild jobs.

Materialized:

- G1-G8: `20 of 20 steps (100%) done`.
- 2026-06-15: `14 of 14 steps (100%) done`.

## Key Acceptance Checks

G1-G8:

- G3 and G5 now have Branch S X-loss overlays in both z and CN plots.
- G3/G5 `event_layer=branch_s_review` covers the chrX non-PAR bins in
  `plot_bins_cn.tsv`.
- G7 XX chrY is marked `sex_aware_absent_expected` and stays neutral.
- G8 XY chrX is evaluated against expected CN = 1, not expected CN = 2.
- G8 chrY is marked `sex_chrom_ref_ratio_not_interpretable` for all chrY bins
  because chrY `chrom_ref_cpm_median=0`; it no longer renders huge dup points.

2026-06-15:

- 5/5 samples have refreshed `.final_cnv.svg`, `.final_cnv_cn.svg`, and
  `.plot_bins_cn.tsv`.
- XY sample chrY bins are marked
  `sex_chrom_ref_ratio_not_interpretable` rather than plotted as huge CN dup.
- 0615 remains no-truth burden/context only; no TP/FN/FP conclusion is made.

## Local Synced Outputs

- `D:\Pipeline\PGT-A\reports\g1_g8_cnv_plots\`
- `D:\Pipeline\PGT-A\reports\0615_cnv_plots\`

## Interpretation Notes

- A strong Branch A event can still have a near-flat bin-level calibrated-z
  line because the z plot is bin-level `calibrated_z`, while Branch A support is
  event-level WisecondorX/CBS evidence.
- Sex-aware CN is visualization-only. It must not be used as a new Branch B or
  Branch S hard filter without a separate truth-safe ablation.
- chrY copy number from the current mixed/low-reference denominator remains
  context-limited. Branch S should use its own evidence table and score summary
  for SCA review; chrY ratio-CN scatter is not a standalone caller.
