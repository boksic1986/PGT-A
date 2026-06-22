# Plot Event Support Ref-Z Consistency

Date: 2026-06-22

Status: `development_only_visualization_support_ledger`

Active reference: `h_r0_shadow_ref_20260619`

Active handoff: `docs/handoff/2026-06-22_1715_plot_event_support_ref_z_handoff.md`

## Scope

This report records a visualization/support-ledger correction for CNV plots and
`plot_event_support.tsv`.

This loop does not modify:

- Branch A calling or merge-gap behavior
- Branch B V2 filtering or report-layer classification
- Branch S classifier
- reference build
- mapping
- report-event counts
- TP/FN/FP metrics

## Problem

`plot_event_support.tsv` had a misleading field:

- `median_calibrated_z` originally meant Branch B residual calibrated-z.
- The plot z signal later changed to Branch A reference-normalized z.
- The field name therefore became stale and easy to over-interpret.

Whole-event median z also diluted signals when a `merge_gap_bp=2Mb` event
contained normal bins. This is expected: the merged event may include normal
gap bins, so whole-segment median can move toward zero even when part of the
event has a clear same-direction Branch A signal.

## Active Support Contract

`plot_bins.tsv`:

- `z` and `display_ref_z` are the plot y-axis values.
- `branch_a_ref_z` is retained as raw Branch A reference-normalized z.
- `residual_calibrated_z` is retained for audit.

Autosomes:

- `display_ref_z = branch_a_ref_z - autosomal_neutral_ref_z_median`
- the neutral median is calculated from autosomal non-event, non-gap, finite
  bins.

Sex chromosomes:

- `display_ref_z = branch_a_ref_z`
- autosomal centering is not applied to chrX/chrY.

Reason: applying the autosomal neutral background to sex chromosomes can flip
Branch S X-loss support direction when the autosomal background is shifted.

`plot_event_support.tsv`:

- `median_calibrated_z` is retained as deprecated compatibility/audit output.
- `median_residual_calibrated_z` is the explicit residual calibrated-z audit
  field.
- The main support field is `same_direction_median_display_ref_z`.
- Same-direction bin counts and fractions are explicit.
- Branch S rows are included with `event_layer=branch_s_review`.

## Added Fields

- `median_raw_branch_a_ref_z`
- `median_display_ref_z`
- `same_direction_median_display_ref_z`
- `same_direction_ref_z_bin_count`
- `same_direction_ref_z_fraction`
- `opposite_direction_ref_z_bin_count`
- `near_zero_ref_z_bin_count`
- `structure_gap_bin_count`
- `centromere_gap_bin_count`
- `median_residual_calibrated_z`
- `median_calibrated_z_status`
- Branch S context fields:
  - `event_layer`
  - `branch_s_state`
  - `sex_call`
  - `sex_chrom_region_class`
  - `expected_copy_number`
  - `median_bin_cn`
  - `same_direction_cn_bin_count`
  - `cn_same_direction_fraction`
  - `support_interpretation_status`

## Remote Validation

Remote path:
`/data/project/CNV/PGT-A/refactor_validation_20260419`

Python executable:
`/biosoftware/miniconda/envs/snakemake_env/bin/python`

Snakemake executable:
`/biosoftware/miniconda/envs/snakemake_env/bin/snakemake`

### Unit Tests

Command:

```bash
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_plot.py \
  tests/unit/test_branch_s_shadow.py \
  tests/unit/test_cnv_report.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py -q
```

Result:

```text
60 passed in 2.57s
```

The new regression test failed before implementation, confirming that autosomal
centering could flip Branch S X-loss support direction.

### Dry-Runs

G1-G8:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake \
  -s Snakefile \
  --configfile config_predict_g1_g8_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml \
  --cores 1 -n --forcerun cnv_branch_s_shadow cnv_branch_ab_plot cnv_report
```

Result: parsed successfully; planned only Branch S, plot, report, and runtime
tracking jobs. No mapping/reference rebuild was requested.

2026-06-15:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake \
  -s Snakefile \
  --configfile config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml \
  --cores 1 -n --forcerun cnv_branch_s_shadow cnv_branch_ab_plot cnv_report
```

Result: parsed successfully; planned only Branch S, plot, report, and runtime
tracking jobs. No mapping/reference rebuild was requested.

### Materialization

G1-G8:

```text
20 of 20 steps (100%) done
log: .snakemake/log/2026-06-22T170626.626839.snakemake.log
```

2026-06-15:

```text
14 of 14 steps (100%) done
log: .snakemake/log/2026-06-22T170651.620498.snakemake.log
```

## Output Acceptance

G1-G8:

- 8 z SVG
- 8 CN SVG
- 8 z bin TSV
- 8 CN bin TSV
- 8 event support TSV

2026-06-15:

- 5 z SVG
- 5 CN SVG
- 5 z bin TSV
- 5 CN bin TSV
- 5 event support TSV

G3 X-loss Branch S support:

```text
same_direction_ref_z_bin_count=170
same_direction_ref_z_fraction=0.977011
support_interpretation_status=Z_DIRECTION_SUPPORTED
```

G5 X-loss Branch S support:

```text
same_direction_ref_z_bin_count=169
same_direction_ref_z_fraction=0.971264
support_interpretation_status=Z_DIRECTION_SUPPORTED
```

## Local Synced Outputs

- `D:\Pipeline\PGT-A\reports\g1_g8_cnv_plots\`
- `D:\Pipeline\PGT-A\reports\0615_cnv_plots\`

These local report files remain untracked and are not part of the commit.

## Interpretation Boundary

Use `same_direction_median_display_ref_z` and
`same_direction_ref_z_fraction` for visual support review.

Do not use deprecated `median_calibrated_z` as the main event support field.

This loop does not prove new performance, reduce FP burden, or promote Branch S
to final SCA calling. It only makes the support table consistent with the z/CN
plot contract.
