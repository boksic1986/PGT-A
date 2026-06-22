# Handoff: Copy Number CNV Plot V2

Date: 2026-06-22 11:00

## 1. Goal

Revise the copy-number CNV plot so manual review can inspect both event-level
CN and bin-level variation. This is visualization-only and does not change
Branch A, Branch B V2, Branch S, report-event classification, filtering,
reference build, or mapping.

## 2. Confirmed Project State

- Active reference remains `h_r0_shadow_ref_20260619`.
- Active 0615 config remains
  `config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`.
- Existing calibrated-z plot remains unchanged.
- The calibrated bin table has no true per-bin copy-number field; therefore
  bin-level CN scatter is an event-anchored visualization proxy, not an
  independent CN caller.

## 3. Completed Changes

Modified:

- `pgta/predict/branch_b/plot.py`
- `tests/unit/test_branch_b_plot.py`
- `docs/reports/report_main_convergence_cnv_plot_2026-06-22.md`
- `CURRENT_STATE.md`
- `PLANS.md`
- `docs/CURRENT_CONTEXT_INDEX.md`

CN V2 plot behavior:

- CN SVG width is now `2560`.
- CN plot uses a dark plotting panel and wider chromosome gaps.
- each chromosome has 50Mb tick labels.
- structural gap / centromere-telomere bins are blanked when
  `is_gap_centromere_telomere` is true or
  `gap_centromere_telomere_overlap_fraction >= 0.5`.
- red report CN trend lines are split across non-gap contiguous chunks.
- final `dup` event scatter points use dark blue.
- final `del` event scatter points use red.
- legend only shows `dup` and `del`.
- no genome-wide smooth CN line or polyline is drawn.
- `plot_bins_cn.tsv` records `copy_number_source`, including
  `event_scaled_calibrated_z_proxy`,
  `event_cn_uniform_median_z_uninformative`,
  `neutral_diploid_baseline`, and `structure_gap_blank`.

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

Result before implementation:

```text
2 failed, 2 passed
```

Final unit tests:

```text
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_plot.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py \
  tests/unit/test_cnv_report.py -q
```

Result:

```text
34 passed in 0.89s
```

Dry-run:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml \
  --cores 1 -n --forcerun cnv_branch_ab_plot cnv_report
```

Result:

```text
Only 5 cnv_branch_ab_plot jobs plus report/runtime refresh were planned.
No mapping or reference rebuild jobs were requested.
```

Materialization:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml \
  --cores 4 --forcerun cnv_branch_ab_plot cnv_report
```

Result:

```text
9 of 9 steps (100%) done
Complete log: .snakemake/log/2026-06-22T105631.231123.snakemake.log
```

## 5. Materialized Output Check

0615 outputs:

- CN SVG count: `5`
- CN bin table count: `5`
- report summary links CN SVGs for 5/5 samples

CN SVG checks:

| sample | width 2560 | trend chunks | gap blanks | 50Mb ticks | dup/del-only legend | polyline |
|---|---|---:|---:|---|---|---|
| JZ26125843-56-56 | true | 15 | 544 | true | true | false |
| JZ26125844-59-59 | true | 26 | 544 | true | true | false |
| JZ26125845-60-60 | true | 47 | 544 | true | true | false |
| JZ26125846-61-61 | true | 42 | 544 | true | true | false |
| JZ26125847-62-62 | true | 34 | 544 | true | true | false |

CN bin table checks:

| sample | rows | structural-gap blank bins | event-scaled proxy bins |
|---|---:|---:|---:|
| JZ26125843-56-56 | 4144 | 544 | 0 |
| JZ26125844-59-59 | 4144 | 544 | 0 |
| JZ26125845-60-60 | 4144 | 544 | 1003 |
| JZ26125846-61-61 | 4144 | 544 | 0 |
| JZ26125847-62-62 | 4144 | 544 | 0 |

Local synced outputs:

- `D:\Pipeline\PGT-A\reports\0615_cnv_plots\*.final_cnv_cn.svg`
- `D:\Pipeline\PGT-A\reports\0615_cnv_plots\*.plot_bins_cn.tsv`

## 6. Current Conclusion

CN V2 is ready for manual 0615 review. It improves visual readability by
showing event-level CN trends, bin-level proxy scatter, structural blank zones,
and wider chromosome spacing.

The CN proxy is not an independent bin-level CN caller. It must not be used for
Branch B V2 filtering, SCA decisions, report-event status changes, or
TP/FP/FN metrics.

## 7. Suggested Next Step

Review sample 56 with the updated z/CN plot pair. Use the red event-level CN
trend against the per-bin CN proxy scatter to identify whether report events
have coherent bin-level support or sparse-bin behavior.

## 8. Key Files And Paths

- plot code: `pgta/predict/branch_b/plot.py`
- plot tests: `tests/unit/test_branch_b_plot.py`
- report: `docs/reports/report_main_convergence_cnv_plot_2026-06-22.md`
- active context: `docs/CURRENT_CONTEXT_INDEX.md`
- local 0615 plots: `D:\Pipeline\PGT-A\reports\0615_cnv_plots`

## 9. Core File Sync

- `CURRENT_STATE.md`: updated with CN V2 contract and validation.
- `PLANS.md`: updated with CN V2 manual review next gate.
- `docs/CURRENT_CONTEXT_INDEX.md`: updated to this handoff.
- `REPO_MAP.md`: not updated; no stable entrypoint or directory structure
  changed.
- `AGENTS.md`: not updated; no hard repository constraint changed.
