# Copy Number Ratio CNV Plot Contract

Date: 2026-06-22

Status: development visualization contract, not production promotion

## Scope

This report supersedes the earlier `calibrated_z_mosaic30_cn_proxy` copy-number
plot contract. The change is limited to visualization and event-support
summaries. It does not change Branch A, Branch B V2 filtering, Branch S,
reference build, mapping, WisecondorX predict, report-event classification, or
truth metrics.

## Problem Corrected

The previous CN scatter used:

`CN_proxy = clip(2 + calibrated_z * 0.05, 0, 4)`

That was not a read-depth ratio or log2R-derived copy number. It could create
contradictions such as a deletion report interval containing blue duplication
points, because `calibrated_z` and event direction do not represent the same
quantity.

The active rule is now:

- `calibrated_z` is only for the z plot.
- CN scatter is derived from sample/reference normalized depth ratio.
- If ratio-derived CN cannot be computed, the CN plot must not fall back to z.

## CN Formula

The formula follows the CNVpro-style interpretation:

`log2R = log2(sample_depth / ref_median)`

`CopyNumber = 2 * 2^log2R`

In this workflow the available bin-level values are on the normalized-signal
scale, where `normalized_signal = log1p(CPM)`. Therefore the implemented plot
contract uses:

`sample_cpm = expm1(normalized_signal)`

`ref_cpm = expm1(ref_bin_stability.ref_median)`

`log2R = log2((sample_cpm + 0.001) / (ref_cpm + 0.001))`

`CN = 2 * 2^log2R`

`plot_bins_cn.tsv` does not clip CN values. SVG display clips only for drawing
and marks out-of-range points with `data-copy-number-out-of-range="true"`.

## Outputs

For each sample:

- `{sample}.final_cnv.svg`: calibrated-z plot, unchanged.
- `{sample}.plot_bins.tsv`: calibrated-z plot bins, unchanged.
- `{sample}.final_cnv_cn.svg`: ratio-derived copy-number plot.
- `{sample}.plot_bins_cn.tsv`: ratio-derived CN bin table.
- `{sample}.plot_event_support.tsv`: per-report-event support summary.

`plot_bins_cn.tsv` includes:

- `chrom`
- `start`
- `end`
- `genome_pos`
- `z`
- `log2r`
- `copy_number`
- `cn_scatter_state`
- `report_state`
- `event_report_state`
- `is_structure_gap_blank`
- `copy_number_source`

`plot_event_support.tsv` includes:

- `valid_bin_count`
- `cn_support_bin_count`
- `cn_same_direction_fraction`
- `median_bin_cn`
- `mean_bin_cn`
- `median_log2r`
- `median_calibrated_z`
- `z_support_bin_count`
- `centromere_gap_bin_count`
- `cn_direction_consistency_status`

CN support is based on ratio-derived CN/log2R, not calibrated z. z support is
kept as a separate audit field.

## Visualization Contract

- CN scatter points are colored by bin-level CN only:
  - CN `< 1.7`: deletion red
  - CN `1.7..2.3`: neutral grey
  - CN `> 2.3`: duplication blue
- Final report event intervals remain visible as separate region shading and
  horizontal event-level CN trend lines.
- Structure/centromere blank bins are not plotted as scatter and are excluded
  from support counts.
- The CN plot remains a visualization and support-audit artifact, not a
  standalone CN caller.

## Remote Validation

Remote path:

`/data/project/CNV/PGT-A/refactor_validation_20260419`

Remote pytest:

```bash
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_plot.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py \
  tests/unit/test_cnv_report.py -q
```

Result:

`38 passed in 1.19s`

Dry-runs:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile config_predict_g1_g8_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml \
  --cores 1 -n --forcerun cnv_branch_ab_plot cnv_report
```

Planned only 8 `cnv_branch_ab_plot` jobs plus report/runtime refresh. No
mapping or reference rebuild jobs were requested.

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml \
  --cores 1 -n --forcerun cnv_branch_ab_plot cnv_report
```

Planned only 5 `cnv_branch_ab_plot` jobs plus report/runtime refresh. No
mapping or reference rebuild jobs were requested.

Materialization:

- G1-G8 completed `12 of 12 steps (100%) done`.
- 2026-06-15 completed `9 of 9 steps (100%) done`.

## G1-G8 Truth-Backed Check

G1-G8 remains the preferred visual check because it has locked truth.

Current truth status:

- locked truth: 10
- preserved: 10
- FN: 0
- hard-suppressed truth: 0

Important interpretation:

- G2 chr8 truth remains `internal_review_event`, not filtered. It is not
  highlighted in the autosomal final report CN main plot because the current
  CN main plot only highlights final autosomal report events.
- G2 chr11 truth remains `report_strong_event` and has CN direction support in
  the new support table.

Selected G1-G8 report-event CN support examples:

| sample | event | visibility | CN support |
|---|---|---|---|
| G1 | chr18 gain | report_strong_event | CN_DIRECTION_SUPPORTED, fraction 1.000 |
| G2 | chr11 loss | report_strong_event | CN_DIRECTION_SUPPORTED, fraction 1.000 |
| G2 | chr4 gain | report_strong_event | CN_DIRECTION_WEAK_OR_MIXED, fraction 0.128 |
| G2 | chr17 gain | report_strong_event | CN_DIRECTION_WEAK_OR_MIXED, fraction 0.031 |
| G4 | chr5 loss | report_strong_event | CN_DIRECTION_SUPPORTED, fraction 1.000 |
| G6 | chr9 gain | report_strong_event | CN_DIRECTION_SUPPORTED, fraction 1.000 |
| G8 | chr16 gain | report_strong_event | CN_DIRECTION_SUPPORTED, fraction 0.528 |

This demonstrates why the new support table is needed: some report events have
strong Branch A evidence but weak ratio-derived CN support and should be
reviewed rather than visually treated as clearly supported.

## 2026-06-15 Context

0615 still has no locked truth labels. The regenerated plots and support tables
are review/context artifacts only; they cannot be used to compute TP/FP/FN or
derive production thresholds.

Local synced outputs:

- `D:\Pipeline\PGT-A\reports\g1_g8_cnv_plots\`
- `D:\Pipeline\PGT-A\reports\0615_cnv_plots\`

## Next Step

Review G1-G8 first. Use the z plot, ratio-derived CN plot, and
`plot_event_support.tsv` together. For any future Branch B V2 demotion/filter
proposal, ablate against Y1-Y8, H1-H6, and G1-G8 with FN=0, H6 chr21 visible,
and G2 truth not filtered.
