# Handoff: Branch A Ref-Z Plot And Sex-Chrom CN Trend

Date: 2026-06-22 16:18

Status: active handoff

## 1. Goal

Fix the combined CNV plot so the z-axis displays Branch A-oriented
read-depth/reference deviation (`branch_a_ref_z`) instead of downstream Branch B
residual `calibrated_z`. Also add sex-aware Branch S CN horizontal trend lines
for interpretable chrX/chrY review events.

## 2. Confirmed Project State

- Current reference remains `h_r0_shadow_ref_20260619`, a fixed shadow baseline
  for R&D only.
- Branch A, Branch B V2 filtering, Branch S classifier, reference, mapping, and
  report-event classification were not changed.
- 0615 still has no locked truth and remains burden/context only.
- G1-G8 remains the first truth-backed plot review cohort for this loop.

## 3. Completed Changes

Modified:

- `pgta/predict/branch_b/plot.py`
- `rules/predict_workflow.smk`
- `tests/unit/test_branch_b_plot.py`
- `tests/unit/test_branch_ab_phase12_workflow_contract.py`

Documentation:

- `docs/reports/branch_a_ref_z_plot_sex_chrom_cn_trend_2026-06-22.md`
- `docs/CURRENT_CONTEXT_INDEX.md`
- `CURRENT_STATE.md`
- `PLANS.md`
- this handoff

Key implementation:

- z plot points now use `branch_a_ref_z`.
- `plot_bins.tsv` keeps `residual_calibrated_z` as audit context.
- report-event z trend lines use Q75 for gain/dup and Q25 for loss/del.
- plot rule always passes `CNV_B_REF_STABILITY_BINS`.
- Branch S CN trend lines are drawn only when sex-aware CN bins are
  interpretable and deviate beyond the configured visual thresholds.

## 4. Current Conclusion

The plot-layer contract is now aligned with Branch A visual review:

- autosomal/whole-chromosome G truth events show visible same-direction ref-z
  deviation;
- G3/G5 X-loss are visible through Branch S z regions and CN trend lines;
- G8 chrY stays uninterpretable instead of producing a giant false dup CN
  trend.

This is visualization-only. It does not change calling or filtering.

## 5. Validation

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

`58 passed in 2.06s`

Dry-runs:

- G1-G8: `cnv_branch_s_shadow cnv_branch_ab_plot cnv_report -n` parsed
  successfully and did not request mapping/reference rebuild.
- 2026-06-15: `cnv_branch_s_shadow cnv_branch_ab_plot cnv_report -n` parsed
  successfully and did not request mapping/reference rebuild.

Materialized:

- G1-G8: `20 of 20 steps (100%) done`.
- 2026-06-15: `14 of 14 steps (100%) done`.

Output checks:

- G1-G8: 8 z SVG, 8 CN SVG, 8 z TSV, 8 CN TSV.
- 2026-06-15: 5 z SVG, 5 CN SVG, 5 z TSV, 5 CN TSV.
- G3/G5 CN SVGs contain Branch S CN trend lines.
- G8 chrY `plot_bins_cn.tsv` marks all chrY bins
  `sex_chrom_ref_ratio_not_interpretable`.

G truth ref-z direction checks are documented in:

`docs/reports/branch_a_ref_z_plot_sex_chrom_cn_trend_2026-06-22.md`

## 6. Local Outputs

Synced local plot directories:

- `D:\Pipeline\PGT-A\reports\g1_g8_cnv_plots\`
- `D:\Pipeline\PGT-A\reports\0615_cnv_plots\`

These local report files are not meant to be committed.

## 7. Remaining Risks

- `branch_a_ref_z` depends on `ref_bin_stability.tsv`; if a config lacks this
  file, plotting should fail rather than fallback.
- The G3/G5 chrX ref-z signal is modest; SCA interpretation still requires
  Branch S summaries and broader SCA validation before production promotion.
- Branch S remains review/development-only.
- 0615 must not be used to derive thresholds because it has no locked truth.

## 8. Recommended Next Step

Review G1-G8 plots first, because they have locked truth. Use:

- z SVG for Branch A ref-z visual signal;
- CN SVG for sex-aware CN context;
- `plot_bins.tsv` for `branch_a_ref_z` and residual `calibrated_z`;
- `plot_bins_cn.tsv` for sex-aware CN interpretation status;
- Branch S summary/scores/evidence for SCA review context.

After manual G review, inspect 0615 only as burden/context and do not derive
new thresholds from it.

## 9. Command Log

Remote pytest:

```bash
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_plot.py \
  tests/unit/test_branch_s_shadow.py \
  tests/unit/test_cnv_report.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py -q
```

Result: `58 passed in 2.06s`.

Remote dry-run and materialization used:

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

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile config_predict_g1_g8_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml \
  --cores 8 --forcerun cnv_branch_s_shadow cnv_branch_ab_plot cnv_report
```

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml \
  --cores 8 --forcerun cnv_branch_s_shadow cnv_branch_ab_plot cnv_report
```

Materialization results:

- G1-G8: `20 of 20 steps (100%) done`.
- 2026-06-15: `14 of 14 steps (100%) done`.

## 10. Core File Sync

- `docs/CURRENT_CONTEXT_INDEX.md`: updated to make this handoff/report the
  active context entrypoint.
- `CURRENT_STATE.md`: updated with the new ref-z plot contract and validation.
- `PLANS.md`: updated with the next review gate.
- `AGENTS.md`: not updated; no repository-level hard rule changed.
- `REPO_MAP.md`: not updated; no stable repository structure changed.
