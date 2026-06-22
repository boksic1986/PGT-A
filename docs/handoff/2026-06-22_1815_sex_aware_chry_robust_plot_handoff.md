# 2026-06-22 18:15 Sex-Aware chrY Robust Plot Handoff

Status: current handoff

Decision use: active context entrypoint

Reference: `h_r0_shadow_ref_20260619`

Report: `docs/reports/sex_aware_chry_robust_plot_truth_sync_2026-06-22.md`

## What Changed

This handoff records a visualization/support-ledger correction only.

Updated files:

- `pgta/predict/branch_b/plot.py`
- `tests/unit/test_branch_b_plot.py`
- `docs/reports/sex_aware_chry_robust_plot_truth_sync_2026-06-22.md`
- `docs/CURRENT_CONTEXT_INDEX.md`
- `CURRENT_STATE.md`
- `PLANS.md`

No Branch A, Branch B V2 filtering, Branch S classifier, reference, mapping, or
report-event classification logic was changed.

## Active Plot Contract

`plot_bins.tsv`:

- `branch_a_ref_z` keeps raw ref-normalized z for audit.
- `display_ref_z` remains the review z signal.
- `plot_z` is the SVG display signal clipped to `[-8, 8]`.
- `z_plot_clipped=true` marks display-only clipping.
- XX chrY rows are retained but hidden from ordinary scatter:
  `chrY_display_mode=xx_absent_expected_hidden`.
- XY chrY rows use a Y-presence guide rather than raw chrY ref-z:
  `chrY_display_mode=xy_y_presence_guide_not_ref_z`.

`plot_bins_cn.tsv`:

- XX chrY rows are retained but copy number is hidden:
  `chrY_display_mode=xx_absent_expected_hidden`.
- XY chrY rows with uninterpretable ref denominator are shown as a neutral
  presence guide at expected CN=1:
  `chrY_display_mode=xy_y_presence_guide_not_ratio_cn`.
- BAM sex evidence columns are carried into plot TSVs:
  `bam_y_relative_depth`, `bam_y_to_x_ratio`, `bam_inferred_sex`.

SVG:

- CN plots no longer use grey structure-gap blocks. Structural/centromere gaps
  are represented by absent scatter points.
- z legend no longer includes `event ref-z trend`.
- chrY guide points have `chry-presence-guide` class and are neutral/context,
  not dup/del calls.

## Remote Commands Executed

Remote path:
`/data/project/CNV/PGT-A/refactor_validation_20260419`

### Tests

```bash
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_plot.py \
  tests/unit/test_cnv_report.py \
  tests/unit/test_branch_s_shadow.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py -q
```

Result:

```text
63 passed in 2.63s
```

### Dry-Runs

Ran `cnv_branch_s_shadow cnv_branch_ab_plot cnv_report -n` on:

- `config_predict_y_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`
- `config_predict_h_20260608_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`
- `config_predict_g1_g8_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`
- `config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`

Result: success; no mapping/reference rebuild requested.

### Materialization

Materialized plot/report outputs:

- Y1-Y8: `logs/manual/sex_aware_chry_plot_Y_20260622_175637.log`
- H1-H16: `logs/manual/sex_aware_chry_plot_H_20260622_175650.log`
- G1-G8: `logs/manual/sex_aware_chry_plot_G_20260622_175712.log`
- 0615 context: `logs/manual/sex_aware_chry_plot_0615_20260622_180102.log`

Results:

- Y1-Y8 and G1-G8 completed `20 of 20 steps (100%) done`.
- 0615 context completed `14 of 14 steps (100%) done`.
- H1-H16 completed; H1-H6 used for truth-gate review.

## Acceptance Evidence

- Y1-Y8 truth preserved: `10/10`.
- H1-H6 truth preserved: `10/10`.
- G1-G8 truth preserved: `10/10`.
- H6 chr21 gain remains `report_weak_event` and visible.
- G3/G5 X-loss Branch S support rows remain visible and
  `Z_DIRECTION_SUPPORTED`.
- XX samples no longer display ordinary chrY z/CN scatter.
- XY samples display neutral Y-presence guide instead of giant chrY dup points.

Local synchronized outputs:

- `D:\Pipeline\PGT-A\reports\truth_y1_y8_cnv_plots\`
- `D:\Pipeline\PGT-A\reports\truth_h1_h6_cnv_plots\`
- `D:\Pipeline\PGT-A\reports\truth_g1_g8_cnv_plots\`
- `D:\Pipeline\PGT-A\reports\0615_cnv_plots\` (context only)

## Boundaries

- This is not a new Branch B filter.
- This is not a Branch S final promotion.
- chrY guide is sex-presence visualization, not chrY CNV calling.
- `plot_z` clipping is display-only and not used by caller/filter/evaluation.
- 0615 remains burden/context only because it has no truth.

## Next Recommended Step

Use the synchronized truth plot folders to visually review Y1-Y8, H1-H6, and
G1-G8. If report burden reduction resumes, do it separately from this plot
contract and keep the locked truth gate at FN=0 with H6 chr21 and G2 preserved.
