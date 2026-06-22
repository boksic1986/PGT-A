# Handoff: Sex-Specific Ref CNV Plot And Review Overlay

Date: 2026-06-22 19:38

## Context

Active repository: `D:\Pipeline\PGT-A`

Remote mirror:
`/data/project/CNV/PGT-A/refactor_validation_20260419`

Active branch:
`codex/cnv-plot-wisecondor-style`

This handoff supersedes
`docs/handoff/2026-06-22_1815_sex_aware_chry_robust_plot_handoff.md` for plot
interpretation. The previous fixed XY chrY `CN=1` / `z=0` guide is no longer
valid.

## Files Modified

- `pgta/predict/branch_b/plot.py`
- `pgta/predict/branch_b/ref_stability.py`
- `rules/predict_layout.smk`
- `rules/predict_workflow.smk`
- `tests/unit/test_branch_b_plot.py`
- `tests/unit/test_ref_stability.py`
- `config_predict_y_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`
- `config_predict_h_20260608_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`
- `config_predict_g1_g8_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`
- `config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`
- `docs/CURRENT_CONTEXT_INDEX.md`
- `CURRENT_STATE.md`
- `PLANS.md`
- `docs/reports/sex_specific_ref_cnv_plot_2026-06-22.md`
- `docs/handoff/2026-06-22_1938_sex_specific_ref_cnv_plot_handoff.md`

## What Changed

### Sex-specific ref stability

`ref_stability.py` can now emit `mixed`, `XX`, and `XY` reference groups when
reference sample sex labels are available. The active lowres/gap2m configs
load sample sexes from:

```text
config_reference_h_r0_shadow_20260619.yaml
```

Autosomes use `mixed`; chrX/chrY prefer the matching `XX` or `XY` group for
plot z and CN.

### chrY display

XY chrY is no longer drawn as a fixed guide. If the sex-specific denominator is
not interpretable, ordinary chrY CNV scatter is not drawn and TSV records
`sex_chrom_cn_unavailable`.

XX chrY is retained as expected absence context and is not colored as dup/del.

### z display

Raw `branch_a_ref_z` remains in TSV. Out-of-range z points are hidden in SVG
instead of clipped into a horizontal top/bottom line. Hidden points are marked
with:

```text
z_plot_status=out_of_range_hidden
```

### Review overlays

Plots now overlay:

- report events;
- internal review events;
- Branch S review events.

This is display/support only. It does not promote any event in `cnv_report`.

## Remote Commands And Results

Unit test command:

```bash
cd /data/project/CNV/PGT-A/refactor_validation_20260419
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

Dry-run command pattern:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile <active lowres/gap2m config> \
  --cores 1 -n branch_s_review cnv_report
```

Passed for Y, H, G, and 0615 active lowres/gap2m configs. No mapping or
reference rebuild jobs were requested.

Materialization command pattern:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile <active lowres/gap2m config> \
  --cores 8 branch_s_review cnv_report
```

No `monitor/runtime.db` existed, so there was no historical runtime estimate.

Materialization logs:

- `logs/plot_sexaware_20260622/y.log`: `36 of 36 steps (100%) done`
- `logs/plot_sexaware_20260622/h.log`: `68 of 68 steps (100%) done`
- `logs/plot_sexaware_20260622/g.log`: `36 of 36 steps (100%) done`
- `logs/plot_sexaware_20260622/d0615.log`: `24 of 24 steps (100%) done`

## Materialized Results

Output counts:

- Y1-Y8: 8 z SVG, 8 CN SVG, 8 support TSV.
- H1-H16: 16 z SVG, 16 CN SVG, 16 support TSV.
- G1-G8: 8 z SVG, 8 CN SVG, 8 support TSV.
- 0615: 5 z SVG, 5 CN SVG, 5 support TSV.

Truth preservation:

- Y1-Y8: 10/10, FN=0.
- H1-H6: 10/10, FN=0.
- G1-G8: 10/10, FN=0.
- 0615: context only, no truth.

Specific checks:

- H4 chr15 internal-review gain is visible in z/CN review plots.
- H5, G3, and G5 chrX X-loss Branch S rows are present in
  `plot_event_support.tsv`.
- G2/G8 XY chrY is not fixed to CN=1 and is not fixed to z=0. Current chrY bins
  are unavailable for ordinary CN scatter when the denominator is not
  interpretable.
- G7/H5/H4 XX chrY is expected absent and not interpreted as dup/del.

Local synced output directories:

- `D:\Pipeline\PGT-A\reports\truth_y1_y8_cnv_plots\`
- `D:\Pipeline\PGT-A\reports\truth_h1_h6_cnv_plots\`
- `D:\Pipeline\PGT-A\reports\truth_g1_g8_cnv_plots\`
- `D:\Pipeline\PGT-A\reports\0615_cnv_plots\`

## Boundaries

- This loop is visualization/support-ledger only.
- Branch A candidate generation is unchanged.
- Branch B V2 filtering/report-layer classification is unchanged.
- Branch S classification is unchanged and not promoted to final SCA.
- Reference build and mapping are unchanged.
- TP/FN/FP metrics are unchanged.

## Next Recommended Step

Review the synced Y/H/G truth plots before using the plot output to design new
rules. If report burden reduction resumes, perform candidate-level ablation on
Y1-Y8, H1-H6, and G1-G8, preserving FN=0, H6 chr21, and G2 truth. If SCA needs
finalization, create a separate Branch S truth gate with sex-specific
reference/PAR/non-PAR validation rather than promoting the current review
visualization.
