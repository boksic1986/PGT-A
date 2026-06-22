# Handoff: Ratio-Derived Copy Number CNV Plot

Date: 2026-06-22 14:25

Status: active handoff

## Context

This handoff supersedes:

- `docs/handoff/2026-06-22_1335_copy_number_cnv_plot_cn_threshold_scatter_handoff.md`

The previous `calibrated_z_mosaic30_cn_proxy` CN plot is no longer the active
copy-number interpretation. It was a visualization proxy and could create
direction contradictions. Do not use it for interpretation, filtering, or
candidate support.

## Active Scope

This loop changed only visualization and support summaries:

- z plot remains `calibrated_z` only.
- CN plot now uses sample/reference ratio-derived CN.
- Branch A unchanged.
- Branch B V2 filtering unchanged.
- Branch S unchanged.
- Reference unchanged.
- Mapping unchanged.
- Report-event classification unchanged.

## Implementation Summary

Changed files:

- `pgta/predict/branch_b/plot.py`
- `rules/predict_layout.smk`
- `rules/predict_workflow.smk`
- `rules/target_assembly.smk`
- `tests/unit/test_branch_b_plot.py`
- `tests/unit/test_branch_ab_phase12_workflow_contract.py`
- `docs/CURRENT_CONTEXT_INDEX.md`
- `CURRENT_STATE.md`
- `PLANS.md`
- `docs/reports/copy_number_ratio_cnv_plot_2026-06-22.md`

New plot outputs:

- `{sample}.final_cnv_cn.svg`
- `{sample}.plot_bins_cn.tsv`
- `{sample}.plot_event_support.tsv`

CN formula:

`sample_cpm = expm1(normalized_signal)`

`ref_cpm = expm1(ref_bin_stability.ref_median)`

`log2R = log2((sample_cpm + 0.001) / (ref_cpm + 0.001))`

`CN = 2 * 2^log2R`

No fallback to `calibrated_z` is allowed for CN.

## Validation

All validation below was executed on remote `fengxian` in:

`/data/project/CNV/PGT-A/refactor_validation_20260419`

Remote pytest command:

```bash
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_plot.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py \
  tests/unit/test_cnv_report.py -q
```

Result:

`38 passed in 1.19s`

G1-G8 dry-run:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile config_predict_g1_g8_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml \
  --cores 1 -n --forcerun cnv_branch_ab_plot cnv_report
```

Result:

- Planned 8 `cnv_branch_ab_plot` jobs plus report/runtime refresh.
- No mapping or reference rebuild jobs were requested.

0615 dry-run:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml \
  --cores 1 -n --forcerun cnv_branch_ab_plot cnv_report
```

Result:

- Planned 5 `cnv_branch_ab_plot` jobs plus report/runtime refresh.
- No mapping or reference rebuild jobs were requested.

G1-G8 materialization:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile config_predict_g1_g8_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml \
  --cores 4 --forcerun cnv_branch_ab_plot cnv_report
```

Result:

- `12 of 12 steps (100%) done`
- log `.snakemake/log/2026-06-22T140807.441193.snakemake.log`

0615 materialization:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml \
  --cores 4 --forcerun cnv_branch_ab_plot cnv_report
```

Result:

- `9 of 9 steps (100%) done`
- log `.snakemake/log/2026-06-22T141139.761259.snakemake.log`

## Result Paths

Remote G1-G8:

`results_g1_g8_h_r0_shadow_ref_20260621/wisecondorx/cnv/plots/`

Remote 0615:

`results_20260615_h_r0_shadow_ref_20260619/wisecondorx/cnv/plots/`

Local synced G1-G8:

`D:\Pipeline\PGT-A\reports\g1_g8_cnv_plots\`

Local synced 0615:

`D:\Pipeline\PGT-A\reports\0615_cnv_plots\`

## Interpretation Notes

G1-G8 is the first truth-backed visual check:

- locked truth: 10
- preserved: 10
- FN: 0
- hard-suppressed truth: 0

G2 chr8 truth remains `internal_review_event`, not filtered. It is not
highlighted in the autosomal final report CN main plot because the main plot
only highlights final autosomal report events. Use the ledger/truth metrics for
internal-review truth events.

Some report events have strong Branch A evidence but weak ratio-derived CN
support. Example rows from G1-G8:

- G2 chr4 gain: CN support fraction 0.128.
- G2 chr17 gain: CN support fraction 0.031.
- G2 chr21 gain: CN support fraction 0.118.

These are now visible as `CN_DIRECTION_WEAK_OR_MIXED` in
`plot_event_support.tsv` and should be reviewed before report release.

Extreme CN points can still occur from real sample/reference depth-ratio
outliers, often in high-ref-MAD or difficult regions. They are no longer z
proxy artifacts. Use `plot_event_support.tsv`, ref stability, and z plots
together; do not treat a single out-of-range CN point as an independent call.

## Next Recommended Step

Review G1-G8 plots first, then 0615:

1. Start with G2 because it has both truth and high burden.
2. Compare z plot, CN plot, report events, and `plot_event_support.tsv`.
3. For future Branch B V2 filtering, only propose candidate-level rules after
   ablation against Y1-Y8, H1-H6, and G1-G8 with FN=0.
