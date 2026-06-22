# Handoff: Copy Number CNV Plot Centromere Scatter Fix

Date: 2026-06-22 12:34

## 1. Goal

Fix the remaining CN plot readability issues reported after the scatter/white
background pass:

- background still looked non-white because too many combined
  gap/centromere/telomere bins were shaded;
- per-bin CN scatter points were present in SVG but too small to read clearly;
- chromosome boundaries needed stronger separation.

This is visualization-only. It does not change Branch A, Branch B V2, Branch S,
report-event classification, filtering, reference build, mapping, or thresholds.

## 2. Confirmed Project State

- Active reference remains `h_r0_shadow_ref_20260619`.
- Active 0615 config remains
  `config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`.
- Current 0615 calibrated-bin inputs do not include centromere-only columns such
  as `is_near_centromere` or `nearest_centromere_distance_bp`.
- The previous materialized CN plots had `3600` scatter points and `544`
  structure blanks per sample, because the combined gap/telomere/centromere
  fallback blanked too many bins.

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

- CN plot panel remains white.
- CN SVG still does not emit `chrom-background`.
- Alternating chromosome fills `#f8fafc` / `#eef2f7` remain absent.
- CN grey blanks now use centromere-only fields when present.
- If centromere-only fields are absent, CN grey blanks use hg19 centromere
  coordinates as visualization fallback.
- Combined `is_gap_centromere_telomere` /
  `gap_centromere_telomere_overlap_fraction` fields are no longer used for CN
  plot grey background shading.
- CN scatter points use `class="cn-bin-scatter"`, radius `2.00`, and a white
  stroke for visibility.
- CN chromosome layout uses a larger inter-chromosome gap and explicit
  `chrom-separator` lines.
- CN trend lines remain horizontal event-level segments colored by event
  direction: `dup=#1d4ed8`, `del=#ef4444`.

## 4. Remote Validation

Remote path:

`/data/project/CNV/PGT-A/refactor_validation_20260419`

Executables:

- python: `/biosoftware/miniconda/envs/snakemake_env/bin/python`
- snakemake: `/biosoftware/miniconda/envs/snakemake_env/bin/snakemake`

TDD red run:

```text
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_plot.py -q
```

Result:

```text
2 failed, 3 passed
```

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
38 passed in 1.05s
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
Complete log: .snakemake/log/2026-06-22T123237.156721.snakemake.log
```

## 5. Materialized Output Check

0615 outputs:

- CN SVG count: `5`
- CN bin table count: `5`

CN SVG checks:

| sample | scatter count | scatter radius | centromere blanks | chrom separators | chrom background | alternating fill |
|---|---:|---:|---:|---:|---|---|
| JZ26125843-56-56 | 4024 | 2.00 | 120 | 24 | false | false |
| JZ26125844-59-59 | 4024 | 2.00 | 120 | 24 | false | false |
| JZ26125845-60-60 | 4024 | 2.00 | 120 | 24 | false | false |
| JZ26125846-61-61 | 4024 | 2.00 | 120 | 24 | false | false |
| JZ26125847-62-62 | 4024 | 2.00 | 120 | 24 | false | false |

Local synced outputs:

- `D:\Pipeline\PGT-A\reports\0615_cnv_plots\*.final_cnv_cn.svg`
- `D:\Pipeline\PGT-A\reports\0615_cnv_plots\*.plot_bins_cn.tsv`

## 6. Current Conclusion

The CN plot now has a white background except centromere blank regions. The
previous broad structure fallback was too visually aggressive and has been
removed from the CN background layer. Per-bin CN scatter is now materially more
visible and chromosome boundaries are more explicit.

## 7. Suggested Next Step

Resume manual review from sample 56 with the updated z/CN plot pair. If later
workflow outputs need centromere-only columns instead of the hg19 coordinate
fallback, propagate `is_near_centromere` / `nearest_centromere_distance_bp` into
`CNV_B_CALIBRATED_BINS` as a separate annotation-contract task.

## 8. Core File Sync

- `CURRENT_STATE.md`: updated.
- `PLANS.md`: updated.
- `docs/CURRENT_CONTEXT_INDEX.md`: updated to this handoff.
- `REPO_MAP.md`: not updated; no stable entrypoint or directory structure
  changed.
- `AGENTS.md`: not updated; no hard repository constraint changed.
