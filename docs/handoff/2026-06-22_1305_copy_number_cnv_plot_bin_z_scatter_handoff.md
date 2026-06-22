# Handoff: Copy Number CNV Plot Bin-Z Scatter Fix

Date: 2026-06-22 13:05

## 1. Goal

Fix the remaining CN plot scatter issue reported during 0615 sample review:
the CN scatter points inside final report events were not reflecting each 1Mb
bin's own `calibrated_z` pattern. This made some events look like a uniform CN
block even when the underlying z bins were heterogeneous.

This is visualization-only. It does not change Branch A, Branch B V2, Branch S,
report-event classification, filtering, reference build, mapping, or thresholds.

## 2. Confirmed Project State

- Active reference remains `h_r0_shadow_ref_20260619`.
- Active 0615 config remains
  `config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`.
- Previous CN plot already had:
  - white plot background;
  - no alternating chromosome backgrounds;
  - centromere-only grey blank fallback;
  - `cn-bin-scatter` points with radius `2.00`;
  - chromosome separators and 50Mb ticks.
- The remaining defect was scatter value calculation inside report events:
  most event bins could be uniform event CN instead of a per-bin z-scaled proxy.

## 3. Completed Changes

Modified:

- `pgta/predict/branch_b/plot.py`
- `tests/unit/test_branch_b_plot.py`
- `tests/unit/test_current_context_index.py`
- `docs/CURRENT_CONTEXT_INDEX.md`
- `CURRENT_STATE.md`
- `PLANS.md`
- `docs/reports/report_main_convergence_cnv_plot_2026-06-22.md`

Behavior after the fix:

- `plot_bins_cn.tsv` includes a `z` column.
- final report-event bins use each bin's own `calibrated_z` to compute visual
  CN proxy:

```text
CN_proxy = 2 + calibrated_z * abs(event_cn - 2) / max(median(abs(event_z)), 0.25)
```

- event bins using this rule are labeled:
  `copy_number_source=event_calibrated_z_scaled_to_cn_proxy`.
- neutral bins outside final report events remain `CN=2`.
- centromere blank bins remain `copy_number_source=structure_gap_blank`.
- event-level horizontal CN trend lines still represent the merged event CN and
  remain separate from the per-bin scatter.
- the CN proxy is for visual review only and must not be used for filtering,
  SCA calls, or performance metrics.

## 4. Remote Validation

Remote path:

`/data/project/CNV/PGT-A/refactor_validation_20260419`

Executables:

- python: `/biosoftware/miniconda/envs/snakemake_env/bin/python`
- snakemake: `/biosoftware/miniconda/envs/snakemake_env/bin/snakemake`

Unit tests:

```text
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_plot.py \
  tests/unit/test_cnv_report.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py \
  tests/unit/test_current_context_index.py -q
```

Result:

```text
39 passed in 1.16s
```

0615 dry-run:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml \
  --cores 1 -n --forcerun cnv_branch_ab_plot cnv_report
```

Result:

```text
Only 5 cnv_branch_ab_plot jobs plus report/runtime refresh were planned.
No mapping or reference build jobs were requested.
```

0615 materialization:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml \
  --cores 4 --forcerun cnv_branch_ab_plot cnv_report
```

Result:

```text
9 of 9 steps (100%) done
Complete log: .snakemake/log/2026-06-22T125623.897147.snakemake.log
```

## 5. Materialized Output Check

0615 CN SVG checks:

| sample | scatter count | scatter radius | centromere blanks | chrom separators | chrom background |
|---|---:|---:|---:|---:|---|
| JZ26125843-56-56 | 4024 | 2.00 | 120 | 24 | false |
| JZ26125844-59-59 | 4024 | 2.00 | 120 | 24 | false |
| JZ26125845-60-60 | 4024 | 2.00 | 120 | 24 | false |
| JZ26125846-61-61 | 4024 | 2.00 | 120 | 24 | false |
| JZ26125847-62-62 | 4024 | 2.00 | 120 | 24 | false |

0615 CN TSV checks:

| sample | has z | event bins | z-scaled event bins | event blank bins | event CN unique values |
|---|---|---:|---:|---:|---:|
| JZ26125843-56-56 | true | 503 | 501 | 2 | 501 |
| JZ26125844-59-59 | true | 923 | 922 | 1 | 922 |
| JZ26125845-60-60 | true | 1274 | 1273 | 1 | 1273 |
| JZ26125846-61-61 | true | 1193 | 1191 | 2 | 1191 |
| JZ26125847-62-62 | true | 1016 | 1014 | 2 | 1014 |

Local synced outputs:

- `D:\Pipeline\PGT-A\reports\0615_cnv_plots\*.final_cnv_cn.svg`
- `D:\Pipeline\PGT-A\reports\0615_cnv_plots\*.plot_bins_cn.tsv`

## 6. Current Conclusion

The CN plot now draws event-region scatter from each bin's own
`calibrated_z`. The event-level CN horizontal line is still preserved as a
region-level summary, but the scatter layer now shows within-event bin
heterogeneity. This directly addresses the review issue that 1Mb bin-level z
variation was not visible in the CN plot.

## 7. Suggested Next Step

Resume manual review from sample 56 with the updated z/CN plot pair. Treat
extreme CN proxy points as high-z visual cues, not literal integer copy-number
calls. Do not use the CN proxy as a filter without a separate locked-truth
ablation.

## 8. Core File Sync

- `CURRENT_STATE.md`: updated.
- `PLANS.md`: updated.
- `docs/CURRENT_CONTEXT_INDEX.md`: updated to this handoff.
- `REPO_MAP.md`: not updated; no stable entrypoint or directory structure
  changed.
- `AGENTS.md`: not updated; no hard repository constraint changed.
