# Branch A Ref-Z Plot And Sex-Chrom CN Trend

Date: 2026-06-22

Status: development visualization contract updated; not a calling, filtering,
reference, or Branch S promotion change.

## Scope

This update fixes the z-axis signal used by the combined genome CNV plot. The
previous plot used downstream Branch B residual `calibrated_z` as the bin-level
z signal. That was useful for residual context, but it did not show the Branch A
read-depth/reference deviation that drives WisecondorX/CBS evidence.

This loop changes only plot/report-support visualization:

- Branch A candidates, Branch B V2 report visibility, Branch S classification,
  reference, mapping, and report-event classification are unchanged.
- `branch_a_ref_z` is visualization-only and must not become a filter without a
  separate truth-safe benchmark.
- Branch S CN trend lines are visualization-only and do not promote Branch S to
  final SCA calling.

## Implemented Contract

### z plot

Per-bin z now uses reference stability:

```text
sample_signal = normalized_signal
ref_signal = ref_bin_stability.ref_median
ref_scale = max(1.4826 * ref_mad, min_ref_scale)
branch_a_ref_z = (sample_signal - ref_signal) / ref_scale
```

`min_ref_scale` is the 10th percentile of autosomal non-gap
`1.4826 * ref_mad`. If a unit-test fixture has no autosomal bins, the code can
fall back to all non-gap bins, but real full-genome runs use the autosomal
scale floor.

Hard/structure gap bins, invalid reference bins, and invalid scale bins are not
drawn. `plot_bins.tsv` records them with `z_source` instead of silently
filling values.

Event trend lines no longer use median z. They use direction-aware quantiles:

- gain/dup: `Q75(branch_a_ref_z)`;
- loss/del: `Q25(branch_a_ref_z)`.

This prevents a merged interval with a tolerated 2Mb normal gap from being
visually compressed toward zero. The gap bins remain normal points; only the
event-level horizontal line uses the direction quantile.

`plot_bins.tsv` now includes:

- `z`
- `branch_a_ref_z`
- `residual_calibrated_z`
- `z_source`
- `ref_z_scale`
- `ref_z_scale_source`

### CN plot

The existing sex-aware CN scatter is retained. Branch S
`sca_report_review_event` intervals now also draw a sex-aware CN trend line when
the sex-chromosome bins are interpretable:

- expected CN = 2: draw only when median CN deviates by at least `0.10`;
- expected CN = 1: draw only when median CN deviates by at least `0.05`;
- expected CN = 0 or chrY reference ratio not interpretable: no CN trend line.

This keeps G3/G5 X-loss visible on chrX while avoiding giant chrY trend lines
when the chrY reference denominator is zero or near zero.

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

`58 passed in 2.06s`

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

Result:

- both parsed successfully;
- no mapping or reference rebuild jobs were requested;
- plot rule now always receives `CNV_B_REF_STABILITY_BINS`.

Materialization:

- G1-G8: `20 of 20 steps (100%) done`.
- 2026-06-15: `14 of 14 steps (100%) done`.

## Materialized Output Checks

G1-G8 remote plot directory:

`results_g1_g8_h_r0_shadow_ref_20260621/wisecondorx/cnv/plots/`

2026-06-15 remote plot directory:

`results_20260615_h_r0_shadow_ref_20260619/wisecondorx/cnv/plots/`

Counts:

| cohort | final_cnv.svg | final_cnv_cn.svg | plot_bins.tsv | plot_bins_cn.tsv |
|---|---:|---:|---:|---:|
| G1-G8 | 8 | 8 | 8 | 8 |
| 2026-06-15 | 5 | 5 | 5 | 5 |

Example bin-source check:

- G1 plot bins: `3480` valid `branch_a_ref_median_mad_z` bins and `664`
  masked/structure-gap unavailable bins.
- 0615 sample 56 plot bins: same valid/unavailable bin counts under the same
  genome bin contract.

G1-G8 locked-truth ref-z direction check:

| sample | event | valid bins | median ref-z | direction quantile | same-direction bins |
|---|---|---:|---:|---:|---:|
| G1 | chr18 gain | 85 | 4.568 | Q75 5.130 | 85 |
| G2 | chr8 gain | 3 | 13.747 | Q75 82.394 | 3 |
| G2 | chr11 loss | 6 | -10.496 | Q25 -13.002 | 5 |
| G3 | chr16 gain | 94 | 3.435 | Q75 4.872 | 94 |
| G3 | chrX loss | 175 | -1.109 | Q25 -1.311 | 170 |
| G4 | chr5 loss | 38 | -15.317 | Q25 -17.182 | 36 |
| G5 | chrX loss | 175 | -0.849 | Q25 -1.021 | 169 |
| G6 | chr9 gain | 134 | 4.454 | Q75 5.327 | 134 |
| G7 | chr22 gain | 41 | 3.462 | Q75 4.129 | 41 |
| G8 | chr16 gain | 94 | 1.759 | Q75 3.560 | 79 |

Interpretation:

- Autosomal and large whole-chromosome truth events now show visible
  Branch A-oriented ref-z deviation.
- G3/G5 chrX loss remains lower amplitude, consistent with mosaic/SCA context;
  its visualization support is combined with Branch S CN trend and Branch S
  summaries rather than treated as an autosomal strong z event.

Branch S / sex-chrom checks:

- G3 and G5 CN SVGs contain Branch S CN trend lines.
- G3 and G5 z SVGs contain Branch S review regions.
- G8 chrY remains `sex_chrom_ref_ratio_not_interpretable`; no giant chrY CN
  value or Branch S CN trend is drawn from an uninterpretable denominator.

## Local Synced Outputs

- `D:\Pipeline\PGT-A\reports\g1_g8_cnv_plots\`
- `D:\Pipeline\PGT-A\reports\0615_cnv_plots\`

Local counts after sync:

- G1-G8: `16` SVG files and `24` TSV files.
- 2026-06-15: `10` SVG files and `15` TSV files.

## Boundaries

- `residual_calibrated_z` is preserved in `plot_bins.tsv` for audit context,
  but it is no longer the z plot y-axis.
- `branch_a_ref_z` and sex-aware CN trend lines are visualization-only.
- No report-event counts, Branch B V2 filtering decisions, Branch S classes, or
  TP/FN/FP metrics are changed by this loop.
- 2026-06-15 remains no-truth burden/context only.
