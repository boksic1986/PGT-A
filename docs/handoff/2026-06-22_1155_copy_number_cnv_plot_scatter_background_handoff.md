# Handoff: Copy Number CNV Plot Scatter And White Background Fix

Date: 2026-06-22 11:55

## 1. Goal

Fix the CN plot readability issue where chromosome-wide alternating background
fills made the plot look like cytoband shading. Keep the CN plot background
white, keep grey only for explicit blank regions, and make every rendered
copy-number bin visible as a scatter point.

This is visualization-only. It does not change Branch A, Branch B V2, Branch S,
report-event classification, filtering, mapping, reference build, or any
threshold.

## 2. Confirmed Project State

- Active reference remains `h_r0_shadow_ref_20260619`.
- Active 0615 config remains
  `config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`.
- Existing calibrated-z plot remains unchanged.
- Copy-number plot values remain event-anchored CN proxy values for visual
  review only.

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
- CN plot panel is white.
- CN SVG no longer emits `chrom-background` rectangles.
- Alternating chromosome fills `#f8fafc` / `#eef2f7` are absent.
- Grey fill `#cbd5e1` is reserved for blank regions.
- `structural_gap_mask()` prioritizes centromere-only columns when present:
  `is_near_centromere`, `near_centromere`, `is_centromere`,
  `is_centromere_bin`, centromere overlap fraction columns, or
  `nearest_centromere_distance_bp <= 5Mb`.
- Current 0615 calibrated-bin inputs do not yet contain centromere-only fields,
  so materialized 0615 CN plots use the existing
  `is_gap_centromere_telomere` / `gap_centromere_telomere_overlap_fraction`
  fallback for grey blank regions.
- Every non-blank bin is rendered as a `cn-bin-scatter` point.
- Scatter colors:
  - neutral: `#64748b`
  - dup: `#1d4ed8`
  - del: `#ef4444`
- CN trend lines remain horizontal event-level segments colored by event
  direction:
  - dup: `#1d4ed8`
  - del: `#ef4444`
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

Expected failures:

- old implementation still emitted `chrom-background`;
- old implementation had no `cn-bin-scatter` class on CN scatter points.

Final unit tests:

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
Complete log: .snakemake/log/2026-06-22T115301.719875.snakemake.log
```

## 5. Materialized Output Check

0615 outputs:

- CN SVG count: `5`
- CN bin table count: `5`

CN SVG checks:

| sample | width 2560 | no chrom background | no alternating fill | scatter count | neutral/dup/del scatter | grey blanks | old red CN trend | polyline |
|---|---|---|---|---:|---|---:|---:|---|
| JZ26125843-56-56 | true | true | true | 3600 | true | 544 | 0 | false |
| JZ26125844-59-59 | true | true | true | 3600 | true | 544 | 0 | false |
| JZ26125845-60-60 | true | true | true | 3600 | true | 544 | 0 | false |
| JZ26125846-61-61 | true | true | true | 3600 | true | 544 | 0 | false |
| JZ26125847-62-62 | true | true | true | 3600 | true | 544 | 0 | false |

Local synced outputs:

- `D:\Pipeline\PGT-A\reports\0615_cnv_plots\*.final_cnv_cn.svg`
- `D:\Pipeline\PGT-A\reports\0615_cnv_plots\*.plot_bins_cn.tsv`

## 6. Current Conclusion

The CN plot now uses a white plotting background and explicit per-bin CN scatter
points. Alternating chromosome background shading has been removed. The current
0615 plots still show grey blank regions based on the existing structure
gap/centromere/telomere fallback because centromere-only columns are not yet
present in the calibrated-bin plot input.

## 7. Suggested Next Step

Resume manual review from sample 56 with the updated z/CN plot pair. If a later
task requires strictly centromere-only grey blanks in materialized plots, first
propagate centromere-only annotation columns into `CNV_B_CALIBRATED_BINS`.

## 8. Core File Sync

- `CURRENT_STATE.md`: updated.
- `PLANS.md`: updated.
- `docs/CURRENT_CONTEXT_INDEX.md`: updated to this handoff.
- `REPO_MAP.md`: not updated; no stable entrypoint or directory structure
  changed.
- `AGENTS.md`: not updated; no hard repository constraint changed.
