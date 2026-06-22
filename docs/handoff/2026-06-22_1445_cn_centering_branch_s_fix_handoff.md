# Handoff: CN Centering And Branch S Fix

Date: 2026-06-22 14:45

Status: active handoff

## What Changed

- CN plot ratio-derived copy number is now autosome median centered.
- `plot_bins_cn.tsv` includes `raw_log2r`, centered `log2r`,
  `copy_number_centering_log2_shift`, and centered `copy_number`.
- Branch S segment corroboration no longer uses mean-based support.
- Weak XX X-loss with Branch A support `10 <= abs(z) < 30` is report-review
  visible, fixing G3/G5 SCA visibility.

No Branch A, autosomal Branch B V2 filtering, reference, mapping, or report
event classification was changed.

## Files Modified

- `pgta/predict/branch_b/plot.py`
- `pgta/predict/branch_s.py`
- `tests/unit/test_branch_b_plot.py`
- `tests/unit/test_branch_s_shadow.py`
- `docs/reports/copy_number_centering_and_sca_fix_2026-06-22.md`
- `docs/CURRENT_CONTEXT_INDEX.md`
- `CURRENT_STATE.md`
- `PLANS.md`

The branch also contains earlier CN plot workflow/report contract changes in:

- `rules/predict_layout.smk`
- `rules/predict_workflow.smk`
- `rules/target_assembly.smk`
- `tests/unit/test_branch_ab_phase12_workflow_contract.py`
- `tests/unit/test_current_context_index.py`

## Remote Validation

Remote path:

`/data/project/CNV/PGT-A/refactor_validation_20260419`

Tests:

```bash
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_plot.py \
  tests/unit/test_branch_s_shadow.py \
  tests/unit/test_cnv_report.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py -q
```

Result:

`51 passed in 1.50s`

Dry-runs:

- `config_predict_g1_g8_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`
- `config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`
- `config_predict_h_20260608_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`
- `config_predict_y_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`

Targets:

`cnv_branch_s_shadow cnv_branch_ab_plot cnv_report_summary branch_s_review cnv_report`

Result:

- All parsed successfully.
- No mapping/reference rebuild jobs requested.

Materialized:

- G1-G8: `21 of 21 steps (100%)`.
- 0615: `15 of 15 steps (100%)`.

## Key Results

CN plot autosome medians after centering:

- G1-G8: all samples autosome median CN = `2.000`.
- 0615: all samples autosome median CN = `2.000`.

G1-G8 Branch S after fix:

- G3: `X_LOSS`, `SCA_REVIEW_WEAK`, `sca_report_review_event`.
- G5: `X_LOSS`, `SCA_REVIEW_WEAK`, `sca_report_review_event`.
- G2/G6/G8 XY chrX Branch A signals are no longer strong SCA calls; they are
  `SCA_NO_CALL` with `sca_filtered_or_sex_consistent_event`.

0615 Branch S after fix:

- JZ26125846-61-61 and JZ26125847-62-62 XY X-gain-like signals are no longer
  strong SCA calls.
- 0615 remains no-truth burden/context only.

## Local Synced Outputs

- `D:\Pipeline\PGT-A\reports\g1_g8_cnv_plots\`
- `D:\Pipeline\PGT-A\reports\0615_cnv_plots\`

## Next Step

Review G1-G8 plots first because they have locked truth. Use z plot, centered CN
plot, and `plot_event_support.tsv` together. Do not use centered CN alone as a
filtering rule until it is ablated against Y1-Y8, H1-H6, and G1-G8 with FN=0.
