# Sex-Specific Ref CNV Plot Update

Date: 2026-06-22

Status: development visualization/support update only.

## Scope

This update fixes CNV plot display and plot support ledgers. It does not change
Branch A, Branch B V2 filtering, Branch S classification, reference build,
report event classification, or TP/FN/FP metrics.

## Problem

The previous plot iteration still had several misleading display behaviors:

- XY chrY could be displayed as a fixed `CN=1` / `z=0` guide. That was a guide,
  not a measured value.
- Raw mixed-reference chrY ratios could produce huge false CN values because
  the chrY denominator can be near zero.
- Extreme z values were clipped into a horizontal ceiling in SVG, which looked
  like real high-level plateaus.
- Some review-needed events, for example H4 chr15 internal-review gain, were
  not overlaid because the plot used final report events only.
- Branch S chrX events needed to align across z plot, CN plot, and
  `plot_event_support.tsv`.

## Implementation

### Reference stability groups

`pgta/predict/branch_b/ref_stability.py` now supports optional reference sample
sex labels and emits three groups when available:

- `mixed`
- `XX`
- `XY`

Autosomes use the `mixed` group. chrX/chrY prefer the matching sex-specific
group for plot z and CN.

The active lowres/gap2m configs now point
`core.wisecondorx.cnv.lowres_evidence.reference_sample_sex_map_config` to
`config_reference_h_r0_shadow_20260619.yaml`.

### Plot z

`branch_a_ref_z` remains the raw ref-normalized audit value in TSV.

SVG display no longer clips outliers to a fixed ceiling. If the display z value
is outside the axis range, the bin is hidden from SVG and marked:

```text
z_plot_status=out_of_range_hidden
```

This prevents artificial horizontal rows at the plot top.

### Sex chromosomes

XY chrY is no longer drawn as a fixed guide. If the XY sex-specific chrY
denominator is not interpretable, ordinary chrY CNV scatter is hidden and TSV
records:

```text
copy_number_interpretation_status=sex_chrom_cn_unavailable
copy_number_source=sex_chrom_ref_ratio_not_interpretable
```

XX chrY is retained in TSV as expected absence context and not colored as
dup/del.

### Review overlays

The review plot overlay now includes:

- final autosomal report events;
- internal-review autosomal events;
- Branch S review events.

This does not promote internal-review or Branch S events into the final report
table. The event layer is recorded in SVG titles and `plot_event_support.tsv`.

### Layout

z and CN SVGs are both 2560 px wide. The genome axis uses light chromosome
background blocks with small chromosome gaps and no vertical chromosome
separator lines. CNVpro-inspired display concepts are used only as review
visualization ideas: sex-specific reference, PAR/non-PAR context,
haploid/diploid centers, and segment-style horizontal lines. CNVcalling.R and
cghFLasso are not introduced.

## Remote Validation

Executed on:

```text
ssh fengxian
/data/project/CNV/PGT-A/refactor_validation_20260419
```

Unit tests:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_plot.py \
  tests/unit/test_ref_stability.py \
  tests/unit/test_branch_s_shadow.py \
  tests/unit/test_cnv_report.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py -q
```

Result:

```text
69 passed in 3.07s
```

Dry-runs:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile <active lowres/gap2m config> \
  --cores 1 -n branch_s_review cnv_report
```

Passed for:

- `config_predict_y_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`
- `config_predict_h_20260608_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`
- `config_predict_g1_g8_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`
- `config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`

Materialization was run in background Snakemake jobs. No `monitor/runtime.db`
was present, so no historical runtime estimate was available.

Logs:

- `logs/plot_sexaware_20260622/y.log`: `36 of 36 steps (100%) done`
- `logs/plot_sexaware_20260622/h.log`: `68 of 68 steps (100%) done`
- `logs/plot_sexaware_20260622/g.log`: `36 of 36 steps (100%) done`
- `logs/plot_sexaware_20260622/d0615.log`: `24 of 24 steps (100%) done`

## Materialized Acceptance

Output counts:

| cohort | z svg | CN svg | support TSV | report rows |
|---|---:|---:|---:|---:|
| Y1-Y8 | 8 | 8 | 8 | 8 |
| H1-H16 | 16 | 16 | 16 | 16 |
| G1-G8 | 8 | 8 | 8 | 8 |
| 0615 | 5 | 5 | 5 | 5 |

Truth preservation:

| cohort | truth events | preserved | FN |
|---|---:|---:|---:|
| Y1-Y8 | 10 | 10 | 0 |
| H1-H6 | 10 | 10 | 0 |
| G1-G8 | 10 | 10 | 0 |
| 0615 | 0 | 0 | 0 |

Key checks:

- H6 chr21 remains visible.
- G2 truth remains preserved.
- H4 chr15 internal-review gain is visible in z/CN review plots.
- H5, G3, and G5 chrX X-loss Branch S rows are present in
  `plot_event_support.tsv` and visible in z/CN overlays.
- G2/G8 XY chrY no longer has fixed CN=1 or fixed z=0.
- G7/H5/H4 XX chrY is expected absent and is not colored as dup/del.
- `ref_stability/summary.json` contains `ref_group_counts` for `mixed`, `XX`,
  and `XY` groups, each with 4144 bins in the active runs.

Local synced outputs:

- `D:\Pipeline\PGT-A\reports\truth_y1_y8_cnv_plots\`
- `D:\Pipeline\PGT-A\reports\truth_h1_h6_cnv_plots\`
- `D:\Pipeline\PGT-A\reports\truth_g1_g8_cnv_plots\`
- `D:\Pipeline\PGT-A\reports\0615_cnv_plots\`

## Remaining Risks

- Some chrY bins remain unavailable for ordinary CN scatter because the current
  sex-specific chrY denominator is still not interpretable enough. This is now
  exposed instead of hidden behind a fake guide.
- Final SCA/chrY quantification still needs a separate Branch S validation
  gate with locked SCA truth, sex-specific reference evaluation, and explicit
  PAR/non-PAR handling.
- The plot is a review visualization. It must not be used as an independent
  caller or filter without a locked-truth ablation.

## Next Gate

Review the synced Y/H/G truth plots first. If report burden reduction resumes,
any candidate-level rule must be ablated against Y1-Y8, H1-H6, and G1-G8 with
FN=0, H6 chr21 visible, and G2 truth not filtered. For SCA, design a separate
truth gate rather than promoting the current visualization output.
