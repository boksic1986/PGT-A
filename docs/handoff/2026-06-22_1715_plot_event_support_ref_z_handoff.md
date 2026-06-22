# 2026-06-22 17:15 Plot Event Support Ref-Z Handoff

Status: current handoff

Decision use: active context entrypoint

Reference: `h_r0_shadow_ref_20260619`

Report: `docs/reports/plot_event_support_ref_z_consistency_2026-06-22.md`

## What Changed

This handoff records a visualization/support-ledger correction only.

Updated files:

- `pgta/predict/branch_b/plot.py`
- `tests/unit/test_branch_b_plot.py`
- `docs/CURRENT_CONTEXT_INDEX.md`
- `CURRENT_STATE.md`
- `PLANS.md`
- `docs/reports/plot_event_support_ref_z_consistency_2026-06-22.md`

No Branch A, Branch B V2 filtering, Branch S classifier, reference, mapping, or
report-event classification logic was changed.

## Current Plot/Support Contract

`plot_bins.tsv`:

- z plot points use `display_ref_z`.
- autosomes use autosomal neutral-background-centered `display_ref_z`.
- chrX/chrY use raw `branch_a_ref_z` as `display_ref_z`; autosomal centering is
  not applied to sex chromosomes.

Reason: autosomal centering can flip Branch S X-loss support direction when the
autosomal background is shifted.

`plot_event_support.tsv`:

- includes autosomal report rows and Branch S rows.
- Branch S rows use `event_layer=branch_s_review`.
- `median_calibrated_z` is retained only as a deprecated compatibility/audit
  column.
- `median_residual_calibrated_z` is the explicit residual calibrated-z audit
  column.
- Main z support is `same_direction_median_display_ref_z`, plus:
  - `same_direction_ref_z_bin_count`
  - `same_direction_ref_z_fraction`
  - `median_raw_branch_a_ref_z`
  - `median_display_ref_z`

Whole-event median remains background context. It is not the main support
metric because tolerated 2Mb merge gaps can include normal bins and dilute the
whole-event median toward zero.

## Remote Commands Executed

Remote path:
`/data/project/CNV/PGT-A/refactor_validation_20260419`

### Tests

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

### Dry-Runs

G1-G8:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake \
  -s Snakefile \
  --configfile config_predict_g1_g8_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml \
  --cores 1 -n --forcerun cnv_branch_s_shadow cnv_branch_ab_plot cnv_report
```

Result: success; no mapping/reference rebuild requested.

2026-06-15:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake \
  -s Snakefile \
  --configfile config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml \
  --cores 1 -n --forcerun cnv_branch_s_shadow cnv_branch_ab_plot cnv_report
```

Result: success; no mapping/reference rebuild requested.

### Materialization

G1-G8:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake \
  -s Snakefile \
  --configfile config_predict_g1_g8_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml \
  --cores 8 --forcerun cnv_branch_s_shadow cnv_branch_ab_plot cnv_report
```

Result:

```text
20 of 20 steps (100%) done
log: .snakemake/log/2026-06-22T170626.626839.snakemake.log
```

2026-06-15:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake \
  -s Snakefile \
  --configfile config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml \
  --cores 8 --forcerun cnv_branch_s_shadow cnv_branch_ab_plot cnv_report
```

Result:

```text
14 of 14 steps (100%) done
log: .snakemake/log/2026-06-22T170651.620498.snakemake.log
```

## Acceptance Evidence

Output counts:

- G1-G8: 8/8 z SVG, 8/8 CN SVG, 8/8 z TSV, 8/8 CN TSV, 8/8 support TSV.
- 0615: 5/5 z SVG, 5/5 CN SVG, 5/5 z TSV, 5/5 CN TSV, 5/5 support TSV.

Branch S support rows:

- G3 chrX X-loss:
  - `same_direction_ref_z_bin_count=170`
  - `same_direction_ref_z_fraction=0.977011`
  - `support_interpretation_status=Z_DIRECTION_SUPPORTED`
- G5 chrX X-loss:
  - `same_direction_ref_z_bin_count=169`
  - `same_direction_ref_z_fraction=0.971264`
  - `support_interpretation_status=Z_DIRECTION_SUPPORTED`

Local synced outputs:

- `D:\Pipeline\PGT-A\reports\g1_g8_cnv_plots\`
- `D:\Pipeline\PGT-A\reports\0615_cnv_plots\`

These local output directories are not tracked and should not be committed.

## Boundaries

- This is not a new Branch B filter.
- This is not a Branch S final promotion.
- 0615 remains burden/context only because it has no truth.
- Plot/support fields are for review and report explanation, not caller logic.

## Next Recommended Step

Review G1-G8 first because it has locked truth. Use:

- z plot
- CN plot
- `plot_event_support.tsv`
- Branch S summary/scores/evidence

For event support, prefer:

- `same_direction_median_display_ref_z`
- `same_direction_ref_z_fraction`

Do not use deprecated `median_calibrated_z` as the primary interpretation
field.
