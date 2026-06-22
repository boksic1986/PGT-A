# Handoff: Copy Number CNV Plot Background Fix

Date: 2026-06-22 11:30

## 1. Goal

Fix the copy-number CNV plot readability regression where the whole plot panel
was rendered with a dark background. The requested behavior is a light
WisecondorX-like plot with grey background only at structural gap /
centromere-telomere bins.

This is visualization-only. It does not change Branch A, Branch B V2, Branch S,
report-event classification, filtering, mapping, reference build, or any
threshold.

## 2. Confirmed Project State

- Active reference remains `h_r0_shadow_ref_20260619`.
- Active 0615 config remains
  `config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`.
- Existing calibrated-z plot remains unchanged.
- Copy-number plot values remain an event-anchored CN proxy for visualization
  only and must not be used for filtering or performance metrics.

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

- CN SVG keeps width `2560`.
- Overall SVG and plotting panel are white/light.
- Chromosome backgrounds match the calibrated-z plot style:
  `#f8fafc` / `#eef2f7`.
- Only structural gap / centromere-telomere bins are drawn as grey blanks
  with `#cbd5e1`.
- CN trend lines are horizontal event-level segments colored by event
  direction:
  - `dup`: `#1d4ed8`
  - `del`: `#ef4444`
- Old all-red CN trend line color `#dc2626` is no longer used for CN trend
  lines.
- Legend remains `dup` / `del` only.

## 4. Remote Validation

Remote path:

`/data/project/CNV/PGT-A/refactor_validation_20260419`

Executables:

- python: `/biosoftware/miniconda/envs/snakemake_env/bin/python`
- snakemake: `/biosoftware/miniconda/envs/snakemake_env/bin/snakemake`

TDD red run after changing tests and before implementation:

```text
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_plot.py -q
```

Result:

```text
2 failed, 2 passed
```

Final unit tests:

```text
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_plot.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py \
  tests/unit/test_cnv_report.py \
  tests/unit/test_current_context_index.py -q
```

Result:

```text
37 passed in 1.01s
```

0615 dry-run:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml \
  --cores 1 -n --forcerun cnv_branch_ab_plot cnv_report
```

Result:

```text
Only 5 cnv_branch_ab_plot jobs, cnv_report_summary, cnv_report,
collect_runtime_tracking, and all were planned. No mapping or reference build
jobs were requested.
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
Complete log: .snakemake/log/2026-06-22T112510.010306.snakemake.log
```

## 5. Materialized Output Check

0615 outputs:

- CN SVG count: `5`
- CN bin table count: `5`

CN SVG checks:

| sample | width 2560 | light backgrounds | grey gap blanks | trend chunks | dup trend | del trend | old red CN trend | polyline |
|---|---|---|---|---:|---:|---:|---:|---|
| JZ26125843-56-56 | true | true | true | 15 | 11 | 4 | 0 | false |
| JZ26125844-59-59 | true | true | true | 26 | 16 | 10 | 0 | false |
| JZ26125845-60-60 | true | true | true | 47 | 28 | 19 | 0 | false |
| JZ26125846-61-61 | true | true | true | 42 | 20 | 22 | 0 | false |
| JZ26125847-62-62 | true | true | true | 34 | 23 | 11 | 0 | false |

Local synced outputs:

- `D:\Pipeline\PGT-A\reports\0615_cnv_plots\*.final_cnv_cn.svg`
- `D:\Pipeline\PGT-A\reports\0615_cnv_plots\*.plot_bins_cn.tsv`

## 6. Current Conclusion

The CN plot background regression is fixed. The plot now keeps the readable
light WisecondorX-like background and uses grey only for structural blank zones.
Trend lines now carry the same dup/del semantic colors requested for the CN
plot.

## 7. Suggested Next Step

Resume manual review from sample 56 using the updated z/CN plot pair. Use the
CN plot for readability only; do not use CN proxy values as an independent
filter or caller.

## 8. Core File Sync

- `CURRENT_STATE.md`: updated to mark the light-background CN plot as current.
- `PLANS.md`: updated to remove the old dark-background current-plan wording.
- `docs/CURRENT_CONTEXT_INDEX.md`: updated to this handoff.
- `REPO_MAP.md`: not updated; no stable entrypoint or directory structure
  changed.
- `AGENTS.md`: not updated; no hard repository constraint changed.
