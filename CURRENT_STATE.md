# CURRENT_STATE.md

## 2026-06-22 Sex-Aware CN/Z Plot And Branch S Overlay

Current handoff:
`docs/handoff/2026-06-22_1535_sex_aware_cn_z_plot_branch_s_overlay_handoff.md`.

Current report:
`docs/reports/sex_aware_cn_z_plot_branch_s_overlay_2026-06-22.md`.

This loop supersedes the plot-layer limitations in
`docs/handoff/2026-06-22_1445_cn_centering_branch_s_fix_handoff.md`. It updates
only visualization inputs and sex-aware CN scatter interpretation. It does not
modify Branch A, autosomal Branch B V2 filtering, Branch S classifier,
reference, mapping, or report-event classification.

Active plot contract:

- z plot remains bin-level `calibrated_z`; Branch A `a_zscore` and Branch S
  `state_score` are event-level support labels only.
- Branch S summary JSON, score TSV, and evidence TSV are now consumed by
  `cnv_branch_ab_plot`.
- Branch S `sca_report_review_event` intervals are drawn on chrX/chrY in the
  combined z and CN plots.
- CN scatter is sex-aware:
  - autosomes expected CN = 2;
  - XX chrX expected CN = 2;
  - XY chrX expected CN = 1;
  - XY chrY expected CN = 1 only if chrY reference context is interpretable;
  - XX chrY is absent/neutral.
- `plot_bins_cn.tsv` now includes `sex_call`, `expected_copy_number`,
  `copy_number_delta`, `cn_scatter_state_sex_aware`,
  `copy_number_interpretation_status`, `sex_chrom_region_class`,
  `event_layer`, and `chrom_ref_cpm_median`.
- chrY bins with chromosome-level low/zero reference denominator are marked
  `sex_chrom_ref_ratio_not_interpretable` and do not produce huge dup points.

Remote validation:

- Remote pytest:
  `tests/unit/test_branch_b_plot.py`,
  `tests/unit/test_branch_s_shadow.py`,
  `tests/unit/test_cnv_report.py`, and
  `tests/unit/test_branch_ab_phase12_workflow_contract.py`: `56 passed in
  1.86s`.
- G1-G8 dry-run for `cnv_branch_s_shadow cnv_branch_ab_plot cnv_report`
  parsed successfully and did not request mapping/reference rebuild.
- 0615 dry-run for `cnv_branch_s_shadow cnv_branch_ab_plot cnv_report` parsed
  successfully and did not request mapping/reference rebuild.
- G1-G8 materialization completed: `20 of 20 steps (100%) done`.
- 0615 materialization completed: `14 of 14 steps (100%) done`.

Materialized acceptance:

- G3/G5 X-loss Branch S events are visible on chrX in both z and CN plots.
- G3/G5 CN TSV marks chrX non-PAR bins as `event_layer=branch_s_review`.
- G7 XX chrY is `sex_aware_absent_expected` and neutral.
- G8 XY chrX is interpreted against expected CN = 1, not diploid CN = 2.
- G8 chrY is `sex_chrom_ref_ratio_not_interpretable` for all chrY bins because
  chrY `chrom_ref_cpm_median=0`.
- 0615 5/5 samples have refreshed z SVG, CN SVG, and CN TSV outputs. 0615
  remains no-truth burden/context only.

Local synced outputs:

- `D:\Pipeline\PGT-A\reports\g1_g8_cnv_plots\`
- `D:\Pipeline\PGT-A\reports\0615_cnv_plots\`

## 2026-06-22 CN Centering And Branch S Fix

Current handoff:
`docs/handoff/2026-06-22_1445_cn_centering_branch_s_fix_handoff.md`.

Current report:
`docs/reports/copy_number_centering_and_sca_fix_2026-06-22.md`.

This loop supersedes the uncentered ratio-CN plot interpretation from
`docs/handoff/2026-06-22_1425_copy_number_ratio_cnv_plot_handoff.md`. It updates
only report visualization/support summaries and Branch S report-layer evidence.
It does not modify Branch A, autosomal Branch B V2 filtering, reference,
mapping, or report-event classification.

Active CN estimate:

- z plot remains unchanged and uses `calibrated_z` only.
- CN plot no longer estimates CN from z.
- `plot_bins_cn.tsv` derives bin CN from normalized depth ratio and then
  recenters each sample by the non-gap autosomal median:
  `sample_cpm = expm1(normalized_signal)`,
  `ref_cpm = expm1(ref_bin_stability.ref_median)`,
  `raw_log2R = log2((sample_cpm + 0.001) / (ref_cpm + 0.001))`,
  `centered_log2R = raw_log2R - median(raw_log2R over non-gap autosomal bins)`,
  `CN = 2 * 2^centered_log2R`.
- `copy_number_source=normalized_signal_ref_median_log2r_autosome_centered`
  marks valid centered ratio-derived bins.
- `copy_number_source=ref_median_unavailable` marks bins where the reference
  denominator is unavailable or zero.
- `copy_number_source=structure_gap_blank` marks centromere/structure blanks.
- CN TSV values are not clipped. SVG display clips out-of-range points only for
  drawing and marks them with `data-copy-number-out-of-range="true"`.
- Scatter color is based on bin CN only: `<1.7` deletion red, `1.7..2.3`
  neutral grey, `>2.3` duplication blue.
- Final report intervals remain visible through separate region shading and
  horizontal event-level CN trend lines.

New support output:

- `{sample}.plot_event_support.tsv` is produced beside each plot.
- It records `valid_bin_count`, `cn_support_bin_count`,
  `cn_same_direction_fraction`, `median_bin_cn`, `mean_bin_cn`, `median_log2r`,
  `median_calibrated_z`, `z_support_bin_count`, `centromere_gap_bin_count`, and
  `cn_direction_consistency_status`.
- CN support is based on ratio-derived CN/log2R, not calibrated z. z support is
  kept as a separate audit field.

Remote validation:

- Remote pytest:
  `tests/unit/test_branch_b_plot.py`,
  `tests/unit/test_branch_s_shadow.py`,
  `tests/unit/test_branch_ab_phase12_workflow_contract.py`, and
  `tests/unit/test_cnv_report.py`: `51 passed in 1.50s`.
- G1-G8, 0615, H1-H16, and Y1-Y8 dry-runs planned only Branch S, plot, report,
  and runtime refresh jobs; no mapping or reference rebuild jobs were requested.
- G1-G8 materialization completed: `21 of 21 steps (100%) done`.
- 0615 materialization completed: `15 of 15 steps (100%) done`.

Materialized acceptance:

- G1-G8: 8/8 z plots, 8/8 CN plots, 8/8 CN bin TSVs, 8/8 event-support TSVs.
- G1-G8 autosomal median CN after centering: all samples `2.000`.
- G1-G8 locked truth remains 10/10 preserved, FN=0, hard-suppressed truth=0.
- G2 chr8 truth remains `internal_review_event`, not filtered; it is not
  highlighted in the autosomal final report CN main plot because the main plot
  only highlights final autosomal report events.
- 2026-06-15: 5/5 z plots, 5/5 CN plots, 5/5 CN bin TSVs, 5/5 event-support
  TSVs. 0615 remains no-truth burden/context only.
- 2026-06-15 autosomal median CN after centering: all samples `2.000`.

Branch S acceptance:

- segment-level support no longer uses mean; median/robust-median support is
  required to avoid sparse-bin XY X-gain artifacts.
- G3 and G5 are now visible as `X_LOSS`, `SCA_REVIEW_WEAK`,
  `sca_report_review_event`.
- G2/G6/G8 XY X-gain-like Branch A signals are `SCA_NO_CALL` with
  `sca_filtered_or_sex_consistent_event`.
- 0615 XY X-gain-like Branch A signals are also no-call/sex-consistent audit;
  0615 remains no-truth context only.

Interpretation notes:

- The previous contradiction where del regions could show blue points was caused
  by z-derived CN proxy. The active CN scatter now uses ratio-derived bin CN,
  so discordance between event direction and scatter is exposed as evidence
  rather than hidden.
- Extreme CN points can still occur when a bin has a very high sample/reference
  depth ratio, often in high-ref-MAD or structurally difficult regions. These
  are shown as out-of-range display points and should be reviewed with
  `plot_event_support.tsv`, not treated as independent true CN calls.

Local synced paths:

- `D:\Pipeline\PGT-A\reports\g1_g8_cnv_plots\`
- `D:\Pipeline\PGT-A\reports\0615_cnv_plots\`

## 2026-06-22 Copy Number CNV Plot Bin-CN Threshold Scatter Fix

Current handoff:
`docs/handoff/2026-06-22_1335_copy_number_cnv_plot_cn_threshold_scatter_handoff.md`.

Current report:
`docs/reports/report_main_convergence_cnv_plot_2026-06-22.md`.

This loop updates only copy-number plot visualization. It does not modify
Branch A, Branch B V2, Branch S, report-event classification, filtering,
reference, mapping, or any threshold.

CN scatter correction:

- z plot remains unchanged: `{sample}.final_cnv.svg` and
  `{sample}.plot_bins.tsv`.
- CN plot remains `{sample}.final_cnv_cn.svg`, with
  `{sample}.plot_bins_cn.tsv`.
- CN plot background remains white; no chromosome alternating background is
  drawn.
- grey blanks remain restricted to centromere visualization fallback regions.
- chromosome separators and 50Mb ticks remain present.
- every non-centromere 1Mb bin is drawn as `class="cn-bin-scatter"` with radius
  `1.35`.
- every non-centromere bin uses its own `calibrated_z` to compute a bounded
  visual CN proxy:
  `CN_proxy = clip(2 + calibrated_z * 0.05, 0, 4)`.
- `plot_bins_cn.tsv` includes `z`, `copy_number`, `report_state`,
  `event_report_state`, and `copy_number_source`.
- `copy_number_source=calibrated_z_mosaic30_cn_proxy` marks bins using the
  current per-bin proxy.
- scatter color is based only on the bin CN proxy, not on the event interval:
  `CN < 1.7` is `del`, `1.7 <= CN <= 2.3` is `neutral`, and `CN > 2.3` is
  `dup`.
- `event_report_state` records overlap with final report events for audit, but
  it does not color the scatter.
- the event-level horizontal CN trend line is still the region-level CN
  estimate and remains separate from the per-bin scatter.
- the CN proxy is for visual review only; it is not a bin-level CN caller and
  must not be used for filtering or performance metrics.
- the previously observed extreme CN values came from the superseded
  event-anchored scaling formula. Example: sample 56 had a bin with
  `calibrated_z=37.920977`, which was amplified to `CN=14.923469`; this was a
  visualization proxy artifact, not a real copy-number estimate.

Remote validation:

- Remote pytest:
  `tests/unit/test_branch_b_plot.py`,
  `tests/unit/test_branch_ab_phase12_workflow_contract.py`,
  `tests/unit/test_cnv_report.py`, and
  `tests/unit/test_current_context_index.py`: `40 passed in 1.19s`.
- 0615 dry-run planned only 5 `cnv_branch_ab_plot` jobs plus report/runtime
  refresh; no mapping or reference rebuild jobs were requested.
- 0615 materialization completed: `9 of 9 steps (100%) done`,
  log `.snakemake/log/2026-06-22T132746.835251.snakemake.log`.

Materialized 0615 acceptance:

- 5/5 `.final_cnv_cn.svg` files exist.
- 5/5 `.plot_bins_cn.tsv` files exist.
- 5/5 CN SVGs contain 4024 `cn-bin-scatter` points, 120 centromere blanks, and
  24 chromosome separators; `chrom-background` is absent and scatter radius is
  `1.35`.
- 5/5 CN TSVs contain `z`, `copy_number`, `report_state`,
  `event_report_state`, and `copy_number_source`.
- 5/5 CN TSVs have 4024 bins with
  `copy_number_source=calibrated_z_mosaic30_cn_proxy` and 120 bins with
  `copy_number_source=structure_gap_blank`.
- Current 0615 CN proxy ranges:
  - `JZ26125843-56-56`: CN 1.894-4.000; neutral=3678, dup=346.
  - `JZ26125844-59-59`: CN 1.911-4.000; neutral=3679, dup=345.
  - `JZ26125845-60-60`: CN 0.000-2.188; neutral=3766, del=258.
  - `JZ26125846-61-61`: CN 1.823-4.000; neutral=3682, dup=342.
  - `JZ26125847-62-62`: CN 1.891-4.000; neutral=3679, dup=345.

Local synced path:

- `D:\Pipeline\PGT-A\reports\0615_cnv_plots\*.final_cnv_cn.svg`
- `D:\Pipeline\PGT-A\reports\0615_cnv_plots\*.plot_bins_cn.tsv`

## 2026-06-22 Copy Number CNV Plot Centromere Background Fix

Current handoff:
`docs/handoff/2026-06-22_1234_copy_number_cnv_plot_centromere_scatter_handoff.md`.

Current report:
`docs/reports/report_main_convergence_cnv_plot_2026-06-22.md`.

This loop updates only copy-number plot visualization. It does not modify
Branch A, Branch B V2, Branch S, report-event classification, filtering,
reference, mapping, or any threshold.

CN V2 visualization contract:

- z plot remains unchanged: `{sample}.final_cnv.svg` and
  `{sample}.plot_bins.tsv`.
- CN plot remains `{sample}.final_cnv_cn.svg`, with
  `{sample}.plot_bins_cn.tsv`.
- CN SVG width is now `2560` to improve chromosome readability.
- CN plot uses a white plotting panel like the calibrated-z plot, wider
  chromosome gaps, explicit chromosome separators, and 50Mb intra-chromosome
  ticks.
- CN plot no longer draws alternating chromosome background rectangles:
  `chrom-background`, `#f8fafc`, and `#eef2f7` are absent from the CN SVG.
- grey background is restricted to centromere blanks. If centromere-only columns
  are present, they have priority; current 0615 calibrated-bin inputs do not yet
  carry centromere-only columns, so the materialized plot uses hg19 centromere
  coordinates as visualization fallback. Combined
  `is_gap_centromere_telomere` / `gap_centromere_telomere_overlap_fraction`
  fields are no longer used for CN plot background shading.
- report CN trend remains a horizontal event-level line, is split across
  non-gap contiguous chunks, and is colored by report direction:
  `dup=#1d4ed8`, `del=#ef4444`.
- CN scatter points are shown per non-centromere bin with SVG class
  `cn-bin-scatter`, radius `2.00`, and a white stroke for visibility:
  - event outside: neutral baseline `CN=2`
  - final dup event: dark blue
  - final del event: red
- CN scatter inside final report events uses an event-anchored proxy:
  `2 + calibrated_z * (event_cn - 2) / median_event_z` when median z is
  informative and direction-consistent.
- if event median z is uninformative, event bins fall back to uniform event CN
  and `copy_number_source=event_cn_uniform_median_z_uninformative`.
- this proxy is for visual review only; it is not a bin-level CN caller and
  must not be used for filtering or performance metrics.

Remote validation:

- TDD red run: `test_branch_b_plot.py` failed on the old implementation because
  scatter points were too small and the combined gap/telomere fallback still
  shaded non-centromere bins.
- Remote pytest:
  `tests/unit/test_branch_b_plot.py`,
  `tests/unit/test_branch_ab_phase12_workflow_contract.py`,
  `tests/unit/test_cnv_report.py`, and
  `tests/unit/test_current_context_index.py`: `38 passed in 1.05s`.
- 0615 dry-run planned only 5 `cnv_branch_ab_plot` jobs plus report/runtime
  refresh; no mapping or reference rebuild jobs were requested.
- 0615 materialization completed: `9 of 9 steps (100%) done`,
  log `.snakemake/log/2026-06-22T123237.156721.snakemake.log`.

Materialized 0615 acceptance:

- 5/5 `.final_cnv_cn.svg` files exist.
- 5/5 `.plot_bins_cn.tsv` files exist.
- 5/5 CN SVGs are width `2560`, contain 50Mb ticks, no alternating
  chromosome backgrounds, centromere-only grey blank regions,
  `cn-bin-scatter` points, explicit chromosome separators, dup/del-only legend,
  no polyline, and state-colored CN trend chunks.
- 5/5 CN SVGs contain 4024 `cn-bin-scatter` points after excluding 120
  centromere blank bins from 4144 input bins.
- CN TSVs contain 4144 rows each and 120 centromere blank bins each.
- Report summary still links 5/5 samples to `.final_cnv_cn.svg`.

Local synced path:

- `D:\Pipeline\PGT-A\reports\0615_cnv_plots\*.final_cnv_cn.svg`
- `D:\Pipeline\PGT-A\reports\0615_cnv_plots\*.plot_bins_cn.tsv`

## 2026-06-22 Copy Number CNV Plot Supplement

Current handoff:
`docs/handoff/2026-06-22_1025_copy_number_cnv_plot_handoff.md`.

Current report:
`docs/reports/report_main_convergence_cnv_plot_2026-06-22.md`.

This loop adds copy-number plots as a visualization supplement only. It does
not modify Branch A, Branch B V2, Branch S, report-event classification,
filtering, reference, mapping, or any threshold.

Implemented contract:

- existing calibrated-z plot remains `{sample}.final_cnv.svg`
- new copy-number plot is `{sample}.final_cnv_cn.svg`
- new CN bin table is `{sample}.plot_bins_cn.tsv`
- CN plot uses event-level copy number:
  - neutral bins are `CN=2`
  - bins inside final autosomal report events use `copy_number_estimate`
  - if needed, event CN can fall back to `sex_adjusted_copy_number`, then
    `CN = 2 * (1 + a_ratio)`
  - CN is not derived from `calibrated_z` or `normalized_signal`
- CN plot legend is limited to `dup`, `del`, `neutral bin`,
  and `report CN trend`
- report summary includes `copy_number_plot_svg`

Remote validation:

- `tests/unit/test_branch_b_plot.py`,
  `tests/unit/test_branch_ab_phase12_workflow_contract.py`, and
  `tests/unit/test_cnv_report.py`: `33 passed in 1.07s`.
- 0615 dry-run planned only 5 `cnv_branch_ab_plot` jobs plus report/runtime
  refresh; no mapping or reference rebuild jobs were requested.
- 0615 materialization completed: `9 of 9 steps (100%) done`,
  log `.snakemake/log/2026-06-22T102013.293346.snakemake.log`.

Materialized 0615 acceptance:

- 5/5 `.final_cnv_cn.svg` files exist.
- 5/5 `.plot_bins_cn.tsv` files exist.
- 5/5 report rows link to `.final_cnv_cn.svg`.
- CN SVGs contain `Copy number`, red `report CN trend`, yellow `dup`, blue
  `del`, grey neutral bins, and no polyline.

Local synced path:

- `D:\Pipeline\PGT-A\reports\0615_cnv_plots\*.final_cnv_cn.svg`
- `D:\Pipeline\PGT-A\reports\0615_cnv_plots\*.plot_bins_cn.tsv`

## 2026-06-22 0615 High-Confidence Report Candidate Review

Current handoff:
`docs/handoff/2026-06-22_0907_0615_high_confidence_report_handoff.md`.

Current report:
`docs/reports/0615_high_confidence_report_candidates_2026-06-22.md`.

This loop only records a read-only review of current 2026-06-15 materialized
report outputs. It does not modify workflow code, thresholds, Branch A, Branch
B V2, Branch S, reference, or remote result files.

Remote input check:

- `report_events.tsv`: 71 autosomal report rows.
- Plot bin TSV files: 5/5 samples present.
- Current plot bin TSVs expose the report plotting signal as `z` and states as
  `dup`, `del`, and `neutral`.

Current conservative high-confidence autosomal review candidates:

- `JZ26125845-60-60`: 10 rows.
- `JZ26125843-56-56`: 0 rows.
- `JZ26125844-59-59`: 0 rows.
- `JZ26125846-61-61`: 0 rows.
- `JZ26125847-62-62`: 0 rows.

Current batch-shared review context:

- `chr4:67.50-101.25Mb gain` is present in 5/5 samples and is excluded from
  high-confidence interpretation.
- `chr4:52.50-67.50Mb gain`, `chr4:101.25-121.50Mb gain`, and
  `chr14:60.75-97.50Mb gain` are 4/5 shared review regions.

This result is development-only review context. 2026-06-15 has no locked truth,
so no TP/FP/FN conclusion is allowed and these observations must not be used to
derive production Branch B V2 filters.

## 2026-06-22 Report Main Convergence And CNV Plot

Current handoff:
`docs/handoff/2026-06-22_1025_copy_number_cnv_plot_handoff.md`.

Current report:
`docs/reports/report_main_convergence_cnv_plot_2026-06-22.md`.

This loop changed only report-layer behavior, Branch B V2 report visibility, and
CNV plot output. It did not modify Branch A, rebuild reference, rerun mapping,
or promote Branch B V2 / Branch S / the shadow reference to production-final.

Implemented contract:

- Branch B V2 report visibility labels are now:
  `report_strong_event`, `report_weak_event`, `internal_review_event`,
  `filtered_event`, and `branch_s_event`.
- The autosomal final report main table consumes V2 benchmark
  `report_events.tsv`, which contains `report_strong_event` and
  `report_weak_event`.
- Internal review rows remain separate from the main table.
- Filtered rows remain audit-only.
- Branch S rows remain in the sex-chromosome/SCA section.
- `cnv_report` now consumes V2 sample summary whenever V2 benchmark is
  available, so zero-autosomal-report-event samples remain visible in
  `cnv_summary.tsv`.
- WisecondorX-style CNV plots now use only `calibrated_z`; missing/non-finite
  bins are skipped and missing `calibrated_z` is an error.
- Plot TSV outputs are written as
  `wisecondorx/cnv/plots/{sample}.plot_bins.tsv` with states limited to
  `dup`, `del`, and `neutral`.
- Plot SVGs no longer draw genome-wide or chromosome-wide smooth polylines.
  They draw red horizontal `report-z-trend` lines only over final autosomal
  report event intervals. Duplication bins are yellow, deletion bins are blue,
  and neutral bins remain grey.
- Copy-number plot SVGs are generated as `{sample}.final_cnv_cn.svg` using
  event-level CN only. They are report visualization supplements and do not
  change any calling/filtering logic.

Remote validation:

- Relevant unit tests: `79 passed`.
- Four active lowres-enabled gap2m configs dry-ran successfully for
  `branch_b_v2_benchmark branch_s_review cnv_report`.
- Reports were rematerialized by forcing the real report rule:
  `--forcerun cnv_report_summary cnv_report`.
- The 0615 plot style refresh was rematerialized with:
  `--forcerun cnv_branch_ab_plot cnv_report`.

Materialized acceptance:

- Y1-Y8: truth 10/10, FN=0, truth filtered=0, report=40, strong=21,
  weak=19, plots 8/8, report rows 8/8.
- H1-H16: truth 10/10, FN=0, truth filtered=0, H6 chr21 remains
  `report_weak_event`, report=23, strong=6, weak=17, plots 16/16, report rows
  16/16.
- G1-G8: truth 10/10, FN=0, truth filtered=0, G2 truth is not filtered,
  report=26, strong=15, weak=11, plots 8/8, report rows 8/8.
- 2026-06-15: no locked truth, burden/context only, report=71, strong=52,
  weak=19, plots 5/5, report rows 5/5.

Current limitation:

- 2026-06-15 final report burden remains high. Further reduction must start
  from candidate-level evidence review and must not use sample-level event
  counts or 2026-06-15 burden counts to reverse-engineer filters.

## 2026-06-22 Lowres Branch B/S Integration

Current handoff:
`docs/handoff/2026-06-22_0930_lowres_branch_bs_integration_handoff.md`.

Current report:
`docs/reports/branch_b_s_lowres_integration_2026-06-22.md`.

This loop connects the completed 2Mb/3Mb shadow references to Branch B and
Branch S as auxiliary evidence only. It does not modify Branch A, does not
rebuild the active 1Mb reference, does not promote `merge_gap_bp=2_000_000` to
default, and does not promote Branch B V2 or Branch S to production-final.

Implemented contract changes:

- `lowres_evidence.reference_npz_glob` is now accepted and expanded for the
  Branch B ref-MAD input contract, avoiding repeated static NPZ lists in cohort
  configs.
- Branch S now receives optional 2Mb/3Mb lowres event inputs when lowres
  evidence is enabled.
- Branch S now separates whole-chromosome X/Y non-PAR median evidence from
  local segment-level X/Y non-PAR evidence:
  - global non-PAR median remains the whole-SCA trend evidence;
  - local segment non-PAR median/mean can preserve small sex-chromosome CNV
    review even when global non-PAR median is neutral;
  - PAR evidence remains secondary context only.
- Branch S lowres context fields are emitted with
  `sex_chrom_lowres_final_impact=context_only_not_filter`.

Remote validation:

- `tests/unit/test_branch_s_shadow.py`,
  `tests/unit/test_branch_ab_phase12_workflow_contract.py`,
  `tests/unit/test_lowres_evidence.py`,
  `tests/unit/test_ref_stability.py`, and
  `tests/unit/test_branch_b_v2_classifier.py`: `68 passed in 1.53s`.
- 2Mb/3Mb predict configs dry-run parsed for Y/H/G/2026-06-15 and planned
  existing-BAM WisecondorX convert/predict jobs, not BWA/mapping.

Lowres predict materialization completed in the remote mirror:

- PID: `4690`.
- log:
  `/data/project/CNV/PGT-A/refactor_validation_20260419/logs/driver/lowres_predict_2mb_3mb_20260622.log`.
- command:
  `/data/project/CNV/PGT-A/refactor_validation_20260419/logs/driver/lowres_predict_2mb_3mb_20260622.command.sh`.
- `monitor/runtime.db` is absent, so no historical duration estimate is
  available.
- Final output counts:
  - Y: 8/8 at 2Mb and 3Mb.
  - H: 16/16 at 2Mb and 3Mb.
  - G: 8/8 at 2Mb and 3Mb.
  - 2026-06-15: 5/5 at 2Mb and 3Mb.
- Lowres-enabled main-chain dry-runs and materialization passed for Y/H/G/0615.
- Materialized acceptance:
  - Y1-Y8: truth 10/10, filtered truth=0, report=21, internal_review=50,
    filtered audit-only=13, Branch S=13.
  - H1-H16: truth 10/10, filtered truth=0, H6 chr21 visible,
    report=6, internal_review=49, filtered audit-only=8, Branch S=42.
  - G1-G8: truth 10/10, filtered truth=0, report=15,
    internal_review=40, filtered audit-only=7, Branch S=13.
  - 2026-06-15: status `skipped_no_truth`, report=52,
    internal_review=76, filtered audit-only=23, Branch S=14.
- H6 chr21 remains `internal_review_event` with
  `LOWRES_2MB_3MB_SAME_DIRECTION_SUPPORT` and `REF_STABILITY_STABLE`.

Pending gate:

- Inspect lowres/ref-MAD evidence by truth event and by remaining report event.
- Use same-direction lowres support to improve confidence/explanation only.
- Do not use lowres absence as standalone demotion/filtering evidence.
- Continue Branch S/SCA validation with a broader locked truth panel.

## 2026-06-22 Branch S / SCA V2 Sex-Aware Review

Current handoff:
`docs/handoff/2026-06-22_0101_branch_s_sca_v2_handoff.md`.

Current report:
`docs/reports/branch_s_sca_v2_sex_aware_review_2026-06-22.md`.

This loop only updates Branch S/SCA report-layer review behavior. It does not
modify Branch A, the active reference, autosomal Branch B V2 report-layer
filtering, or the default `merge_gap_bp=0` behavior. It also does not use the
2Mb/3Mb low-resolution evidence outputs.

Current Branch S behavior:

- Branch A chrX/chrY candidates are anchors only.
- Strong Branch A sex-chromosome support must be corroborated by X/Y non-PAR
  median or robust bin evidence before it becomes strong SCA review.
- Normal XY samples with chrX Branch A support but neutral X non-PAR evidence
  are routed to `sca_filtered_or_sex_consistent_event`, not `X_GAIN` /
  `SCA_REVIEW_STRONG`.
- `XX + X_LOSS` with very strong Branch A support remains visible as review to
  preserve current locked X-loss truth evidence.
- Branch S still uses `review_development_only` / `review_reportable_with_limitations`
  and is not final SCA.

Remote materialized active gap2m Branch S results:

- Y1-Y8: 1 strong SCA review, 7 no-call; Y3 remains visible as `X_LOSS`,
  `SCA_REVIEW_STRONG`, `sca_report_review_event`.
- H1-H16: 2 strong SCA review, 1 weak SCA review, 13 no-call; H5 and H6 remain
  visible as `X_LOSS`, `SCA_REVIEW_STRONG`, `sca_report_review_event`.
- H7-H16: all normal XY context samples are `SCA_NO_CALL` and no longer produce
  false `X_GAIN` / `SCA_REVIEW_STRONG`.
- G1-G8: 2 weak SCA internal-review events and 6 no-call; context only for SCA
  because no locked SCA truth is currently assigned here.
- 2026-06-15: 2 weak SCA internal-review events and 3 no-call; burden/context
  only, no TP/FN/FP conclusion.

Current report-layer totals from `wisecondorx/cnv/report/cnv_summary.tsv`:

- Y1-Y8: report=21, internal_review=50, filtered audit-only=13,
  Branch-S-routed=13, sample burden flags=2.
- H1-H16: report=6, internal_review=49, filtered audit-only=8,
  Branch-S-routed=42, sample burden flags=1.
- G1-G8: report=15, internal_review=40, filtered audit-only=7,
  Branch-S-routed=13, sample burden flags=1.
- 2026-06-15: report=52, internal_review=76, filtered audit-only=23,
  Branch-S-routed=14, sample burden flags=5.

`Branch-S-routed` counts are sex-chromosome candidate routing burden, not SCA
positive counts. The SCA correction reduces false strong SCA interpretation
inside the Branch S section; it does not by itself reduce autosomal report
events.

Remote validation:

- `tests/unit/test_branch_s_shadow.py`,
  `tests/unit/test_cnv_report.py`,
  `tests/unit/test_branch_ab_phase12_workflow_contract.py`, and
  `tests/unit/test_current_context_index.py`: `38 passed in 0.92s`.
- `branch_s_review cnv_report` dry-run returned `RC=0` for all four active
  gap2m configs: Y, H, G, and 2026-06-15.
- Branch S/report outputs were forced and rematerialized for all four active
  gap2m configs after the report-layer reason review fix; `cnv_report`
  completed for Y=8, H=16, G=8, and 2026-06-15=5 samples.

Remaining limitation:

- SCA final promotion still needs a broader locked truth panel covering X gain,
  XXY, XYY, Y loss, mosaic SCA fraction series, and cross-batch clean XX/XY
  negatives.

## 2026-06-21 Current Context Index Override

Current context entrypoint:
`docs/CURRENT_CONTEXT_INDEX.md`.

This override applies to the whole R&D cycle, not only G1. The purpose is to
prevent future work from resuming from stale handoffs, legacy Branch B results,
or old N0/N1 assumptions.

Current fixed baseline:

- `active_reference_id=h_r0_shadow_ref_20260619`.
- Reference status is `fixed_shadow_baseline_not_production`.
- P2 Branch A no-FN validation passed under this reference for Y1-Y8 and H1-H6.
- P2 also exposed high Branch A candidate burden, so the next Branch A work is
  burden optimization under no-FN and H6 chr21 preservation constraints.

Remote mirror / Snakemake parser state:

- The 2026-06-21 remote Snakemake parse blocker was traced to CRLF line endings
  in workflow source files copied from the Windows worktree into the Linux
  mirror.
- Remote mirror workflow sources were normalized to LF:
  `Snakefile`, `rules/*.smk`, `pgta/**/*.py`, `scripts/**/*.py`,
  `tests/**/*.py`, and root YAML configs.
- `.gitattributes` now records the LF contract for `Snakefile`, `*.smk`,
  Python, and YAML files.
- The previously failing dry-run now parses successfully on `ssh fengxian`.

Current Branch B state:

- Branch B remains necessary for candidate-level refinement, including mosaic,
  LOH, UPD, CN-like amplitude, region-risk, sample-noise, background, and
  CNVpro-inspired consistency evidence.
- Legacy/current-code Branch B outputs are excluded from Branch B V2 performance
  evidence. `final_disposition`, `branch_b_keep_event`, legacy artifact labels,
  old `N0=0`, and N1-only matched-negative promotion must not drive V2
  decisions, benchmarks, or report release.
- The next Branch B gate is a V2-only truth benchmark and evidence/disposition
  path using Y1-Y8 and H1-H6 truth. Reference-cohort background can be used only
  as labeled limited context until formal N0/cross-fit evidence exists.

Current Branch S state:

- Branch S is not final SCA.
- Branch S must still be carried into report development as
  `review_reportable_with_limitations`, with visible sex-chromosome/SCA output
  mode and uncertainty.
- The immediate Branch S target is to reduce SCA/sex-chromosome FP burden in
  reference or negative-like samples while preserving sex-call concordance.

Current P6/report state:

- The already materialized 2026-06-15 P6 package remains historical
  `development_only_not_final_release` evidence.
- Future P6/report remains the delivery target after Branch A burden
  optimization and Branch B V2 evidence/disposition strengthening. It should not
  be permanently classified as "cannot produce a report"; it should be generated
  by workflow once the fixed A/B/S contracts are represented and limitations are
  explicit.

Current Branch B V2 burden stratification:

- Current report:
  `docs/reports/branch_b_v2_burden_stratification_2026-06-21.md`.
- This loop did not modify Branch A and did not promote `merge_gap_bp=2_000_000`
  to default. It stratified burden on top of the already materialized explicit
  gap2m Branch A candidates.
- New V2 burden fields include `v2_burden_reduction_tier`,
  `v2_burden_reduction_action`, `v2_burden_reduction_reason`, and
  `v2_burden_evidence_tags`.
- CNVpro/CNVseq borrowed concepts are explicitly tagged:
  `[CNVpro-inspired]` length tiers, `[CNVpro-confirmed]` acrocentric context,
  `[CNVseq-asset]` mask/mappability/sex-homology context, `[CNVpro-like]`
  GC/RC context, and `[Not used]` CNVcalling.R/cghFLasso.
- Materialized truth preservation remains:
  - Y1-Y8: 97 candidates; truth preserved 10/10; FN=0; hard-suppressed truth=0.
  - H1-H16: 105 candidates; truth preserved 10/10; FN=0; hard-suppressed
    truth=0; H6 chr21 retained.
  - 2026-06-15: 165 candidates; no locked truth; burden/context only.
- Current burden stratification does not yet reduce total review burden:
  - Y1-Y8: 84 background-unknown review, 13 Branch S review.
  - H1-H16: 63 background-unknown review, 42 Branch S review.
  - 2026-06-15: 151 background-unknown review, 14 Branch S review.
- Therefore this loop improves auditability and report/review routing context,
  but it is not a final Branch B V2 FP-reduction proof and not report-promotion
  evidence.

Current Branch B V2 report-contract integration:

- Current report:
  `docs/reports/branch_b_v2_report_contract_2026-06-21.md`.
- `cnv_report` is being updated to consume the V2 benchmark summary and
  sample-summary table as report-display evidence only.
- Required report fields are:
  `branch_b_v2_burden_status`,
  `branch_b_v2_background_unknown_review_count`,
  `branch_b_v2_branch_s_review_count`,
  `branch_b_v2_technical_risk_review_count`,
  `branch_b_v2_report_candidate_count`,
  `branch_b_v2_legacy_fields_used=false`, and
  `branch_b_v2_final_impact=development_review_only`.
- This integration does not change Branch A, does not add hard filters, does
  not use legacy/current-code Branch B decisions, and does not promote Branch B
  V2 or Branch S.
- Remote validation completed on `ssh fengxian`:
  - `tests/unit/test_cnv_report.py`,
    `tests/unit/test_branch_ab_phase12_workflow_contract.py`, and
    `tests/unit/test_current_context_index.py`: `28 passed in 0.93s`.
  - Snakemake dry-runs passed for `branch_b_v2_benchmark branch_s_review
    cnv_report` under all three active gap2m configs.
  - Materialized `cnv_report` for Y1-Y8, H1-H16, and 2026-06-15.
- Materialized report checks:
  - Y1-Y8 V2 truth preserved 10/10, FN=0, hard-suppressed truth=0.
  - H1-H16 V2 truth preserved 10/10, FN=0, hard-suppressed truth=0; H6 chr21
    remains visible.
  - 2026-06-15 V2 benchmark status is `skipped_no_truth`, so it remains
    burden/context only.

Current Branch B V2 report-layer filtering:

- Current report:
  `docs/reports/branch_b_v2_report_layer_filter_2026-06-21.md`.
- This loop is the first Branch B V2 step that performs actual report-layer
  filtering rather than display-only burden stratification.
- New V2 report-layer classes are `report_event`, `internal_review_event`,
  `filtered_event`, and `branch_s_event`.
- `filtered_event` rows are removed from report/internal-review main flow and
  retained only in `filtered_events.tsv/json` audit ledgers.
- The initial non-contract filter requires combination evidence:
  non-strong A signal, B-side non-support, GC/RC attenuation, and sensitive,
  short/focal, or low-clean/high-risk context. Single indicators remain
  insufficient for filtering.
- Strong A signals with B/GC disagreement are retained as internal review, not
  filtered. H6 chr21 remains protected as internal review.
- `cnv_report` now prefers V2 report-layer counts in
  `biological_candidate_conclusion` when V2 summaries are present. Legacy/current-code
  Branch B top events are labeled as legacy context and are not the main
  biological conclusion.
- Remote validation completed on `ssh fengxian`:
  - `tests/unit/test_branch_b_v2_classifier.py`,
    `tests/unit/test_branch_b_v2_benchmark.py`,
    `tests/unit/test_cnv_report.py`, and
    `tests/unit/test_branch_ab_phase12_workflow_contract.py`: `66 passed in
    1.15s`.
  - Snakemake dry-runs passed for `branch_b_v2_benchmark branch_s_review
    cnv_report` under all three active gap2m configs.
  - Forced materialization completed for Y1-Y8, H1-H16, and 2026-06-15.
- Materialized report-layer result:
  - Y1-Y8: candidates=97, truth preserved 10/10, FN=0,
    report-layer filtered truth=0, report=21, internal_review=50,
    filtered=13, Branch S=13.
  - H1-H16: candidates=105, truth preserved 10/10, FN=0,
    report-layer filtered truth=0, H6 chr21 retained, report=6,
    internal_review=49, filtered=8, Branch S=42.
  - 2026-06-15: candidates=165, no locked truth, report=52,
    internal_review=76, filtered=23, Branch S=14.
- This is still `development_review_only`; it is not production promotion
  because 2026-06-15 has no truth labels, Branch S is not final, and the active
  reference remains shadow.

G1-G8 current-scheme validation has also been materialized under the same
active reference and explicit gap2m Branch A overlay:

- Current report:
  `docs/reports/branch_b_v2_g1_g8_validation_2026-06-21.md`.
- Existing BAMs from
  `/data/project/CNV/PGT-A/g_reseq_qc_20260504/mapping/G*.sorted.bam` were
  reused; no BWA remapping was performed.
- WisecondorX convert/predict, Branch A, Branch B V2, Branch S, benchmark, and
  report were rerun in the remote mirror.
- Materialized result:
  - candidates=75;
  - locked truth events=10;
  - truth preserved 10/10;
  - FN=0;
  - report-layer filtered truth=0;
  - report=15;
  - internal_review=40;
  - filtered audit-only=7;
  - Branch S=13.
- G1-G8 supports the current sensitivity/preservation gate. It does not promote
  the shadow reference, Branch B V2, or Branch S to production-final status.
- G2 remains the main G-batch burden outlier with 28 candidates, 6 report
  events, 16 internal-review events, 3 filtered audit-only events, and 3
  Branch S events.

Current Branch B V2 pass2 correction:

- Current report:
  `docs/reports/branch_b_v2_pass2_correction_2026-06-21.md`.
- Superseded report:
  `docs/reports/branch_b_v2_report_event_audit_2026-06-21.md`.
- The previous pass2 demotion rule used same-sample report-event burden as part
  of candidate-level evidence. That boundary is now retracted.
- `sample report_event count >= 3` is retained only as a sample-level burden
  audit flag in benchmark summaries. It must not modify
  `v2_report_layer_class`, `v2_report_visibility`, `v2_filter_action`, or
  `v2_burden_reduction_action`.
- The current correction does not modify Branch A, does not rebuild the
  reference, does not promote `merge_gap_bp=2_000_000` to default, and does not
  promote Branch B V2 or Branch S to production-final.
- Full remote rematerialization for Y/H/G/2026-06-15 is complete. Corrected
  materialized result:
  - Y1-Y8: candidates=97, truth preserved 10/10, FN=0, report-layer filtered
    truth=0, report=21, internal_review=50, filtered audit-only=13,
    Branch S=13, sample burden flags=2.
  - H1-H16: candidates=105, truth preserved 10/10, FN=0, report-layer filtered
    truth=0, H6 chr21 retained as internal review, report=6,
    internal_review=49, filtered audit-only=8, Branch S=42,
    sample burden flags=1.
  - G1-G8: candidates=75, truth preserved 10/10, FN=0, report-layer filtered
    truth=0, report=15, internal_review=40, filtered audit-only=7,
    Branch S=13, sample burden flags=1. G2 truth remains visible.
  - 2026-06-15: candidates=165, no locked truth, context only, report=52,
    internal_review=76, filtered audit-only=23, Branch S=14,
    sample burden flags=5.
- Report summaries now display sample report-burden flags as audit-only
  information. Previous pass2 counts remain historical methodology-audit output
  only and must not be used as current V2 performance evidence.

Current Branch B V2 low-resolution/ref-MAD evidence interface:

- Current report:
  `docs/reports/branch_b_v2_lowres_ref_evidence_2026-06-21.md`.
- This loop adds a config-gated auxiliary evidence path for 2Mb/3Mb
  low-resolution same-direction support and reference-bin MAD stability.
- The active Branch A/WisecondorX/CBS chain remains the primary discovery
  layer. The low-resolution path does not replace WisecondorX predict/CBS, does
  not create B-only events, does not modify Branch A, and does not promote
  Branch B V2 or Branch S.
- The new optional workflow path emits low-resolution support fields and
  ref-MAD stability fields before `cnv_branch_b_v2_classifier`:
  `lowres_2mb_*`, `lowres_3mb_*`, `lowres_consensus_label`,
  `event_ref_mad_median`, `event_ref_mad_p90`,
  `high_ref_mad_bin_fraction`, and `ref_stability_context`.
- The V2 classifier exposes the evidence as
  `v2_lowres_context_label`, `v2_lowres_context_reason`, and
  `v2_ref_stability_context`, and adds `[lowres-shadow]` /
  `[ref-mad-shadow]` audit tags when configured.
- Low-resolution absence remains context only. It does not change
  `v2_report_layer_class`, `v2_report_visibility`, `v2_filter_action`, or
  `v2_burden_reduction_action`.
- `core.wisecondorx.cnv.lowres_evidence.enable=true` now requires
  `lowres_evidence.reference_npz`; if `reference_sample_ids` is supplied, its
  length must match the NPZ list.
- Default behavior remains unchanged when lowres evidence is disabled.
- This loop has not built `h_r0_shadow_ref_20260619_2mb` or
  `h_r0_shadow_ref_20260619_3mb`, and has not materialized low-res predict
  events for Y/H/G/2026-06-15. Those are the next long-task gate, not a
  completed result.

## 2026-06-21 Branch A Burden Optimization Phase 1

Current report:
`docs/reports/branch_a_burden_optimization_phase1_2026-06-21.md`.

This phase tightened the Branch A workflow contract without changing the
default Branch A behavior:

- Branch A is still WisecondorX predict/CBS-derived candidate discovery.
- Current WisecondorX predict parameters remain `zscore=5`, `alpha=0.001`,
  `maskrepeats=5`, and `minrefbins=150`.
- `core.wisecondorx.cnv.branch_a.merge_gap_bp` and
  `core.wisecondorx.cnv.branch_a.strong_z` are now passed through the workflow
  into `pgta.predict.branch_a`.
- Branch A output, validation, and log directories are now configurable so
  candidate parameter runs can be materialized without overwriting the default
  P2 evidence.
- Defaults are `merge_gap_bp=0` and `strong_z=10.0`, so current P2 behavior is
  preserved.
- No Branch B artifact, matched-negative, calibration, or legacy disposition
  field is used by Branch A.

Remote validation under default `merge_gap_bp=0`:

- Y1-Y8: truth detected 10/10, FN=0, Branch A candidates=131.
- H1-H16: truth detected 10/10, FN=0, H6 chr21 detected, Branch A
  candidates=221.
- 2026-06-15: no truth labels, Branch A candidates=201.

Ablation conclusions:

- Hard z tightening is unsafe as a Branch A inclusion filter: `z>=8` would miss
  H6 chr21 gain and `z>=10` would miss Y6 chr7 loss plus H6 chr21 gain.
- Existing-BED merge-gap ablation supports `merge_gap_bp=2_000_000` as the
  first candidate for a formal burden-reduction run.
- The explicit 2 Mb overlay has now been materialized remotely without
  overwriting the default P2 outputs:
  - Y1-Y8: truth detected 10/10, FN=0, Branch A candidates=97.
  - H1-H16: truth detected 10/10, FN=0, H6 chr21 detected, Branch A
    candidates=105.
  - 2026-06-15: no truth labels, Branch A candidates=165.
- `merge_gap_bp=2_000_000` is still not promoted to default. Before promotion,
  the fixed Branch B V2 / Branch S / report chain must be rerun under this
  explicit config and benchmarked.

## 2026-06-21 Branch B V2 Gap2m Benchmark

Current report:
`docs/reports/branch_b_v2_gap2m_benchmark_2026-06-21.md`.

This phase materialized the first Branch B V2-only benchmark on top of the
explicit Branch A `merge_gap_bp=2_000_000` overlay. It does not promote
Branch B V2, Branch S, the report package, or the shadow reference.

Workflow contract:

- New target: `branch_b_v2_benchmark`.
- Benchmark scope: `v2_classifier_rows_only`.
- Legacy/current-code Branch B decision fields are explicitly ignored:
  `final_disposition`, `branch_b_keep_event`, `branch_b_report_class`, and
  `branch_b_artifact_status`.
- `final_report_impact=none_shadow_only`.

Remote materialization after correcting the benchmark hard-suppression
semantics:

- Y1-Y8: candidates=97, truth preserved=10/10, FN=0,
  hard-suppressed truth=0.
- H1-H16: candidates=105, truth preserved=10/10, FN=0,
  hard-suppressed truth=0; H6 chr21 remains preserved with
  `top_a_abs_zscore=7.1135`.
- 2026-06-15: candidates=165, no locked truth table; burden/context only.

Important interpretation:

- V2 now proves it can avoid hard suppression of locked truth-overlap
  candidates under the gap2m overlay.
- It does not yet prove FP/review burden is low enough for final reporting.
- Sex-chromosome truth events such as Y3/H5/H6 chrX loss are preserved as
  `V2_NO_CALL_CONTRACT_RISK` / `V2_REVIEW_NO_HARD_SUPPRESSION`, not as final
  SCA calls. Branch S remains a separate review-reportable-with-limitations
  gate.

## 2026-06-21 Branch B V2 Sex-Route Refinement

Current report:
`docs/reports/branch_b_v2_sex_route_refinement_2026-06-21.md`.

This refinement separates sex-chromosome candidates from autosomal Branch B V2
positive-support review:

- `chrX`/`chrY` candidates now route to `V2_SEX_CHROMOSOME_REVIEW` with action
  `V2_ROUTE_BRANCH_S_REVIEW`.
- The original evidence tier, evidence gate, review priority, and
  `none_shadow_only` final-report impact are preserved.
- This does not change WisecondorX predict, Branch A candidate generation,
  mosaic logic, sex calling, final SCA promotion, or report release status.

Remote materialization under the explicit gap2m overlay produced:

- Y1-Y8: 97 candidates; class counts are 78 autosomal
  `V2_POSITIVE_SUPPORT_REVIEW`, 13 `V2_SEX_CHROMOSOME_REVIEW`, and 6
  `V2_NO_CALL_CONTRACT_RISK`; truth preserved 10/10, FN=0,
  hard-suppressed truth=0.
- H1-H16: 105 candidates; class counts are 57 autosomal
  `V2_POSITIVE_SUPPORT_REVIEW`, 42 `V2_SEX_CHROMOSOME_REVIEW`, and 6
  `V2_NO_CALL_CONTRACT_RISK`; truth preserved 10/10, FN=0,
  hard-suppressed truth=0.
- H7-H16 context subset: 57 candidate rows, including 31 sex-chromosome review
  rows and 24 autosomal positive-support review rows; no locked truth, so this
  is burden/context only.
- 2026-06-15: 165 candidate rows, including 14 sex-chromosome review rows and
  127 autosomal positive-support review rows; no locked truth, so this remains
  burden/context only.

Important interpretation:

- Y3/H5/H6 chrX truth events are preserved but now route to Branch S review
  instead of autosomal Branch B positive-support review.
- H6 chr21 gain remains autosomal `V2_POSITIVE_SUPPORT_REVIEW` with
  `top_a_abs_zscore=7.1135`.
- Existing benchmark summary field `v2_positive_support_candidate_count` is
  evidence-tier based and may include sex-chromosome positive-support evidence;
  class-level autosomal burden should be read from
  `v2_classifier/*.candidate_classification.tsv` or the sex-route report table.

## 2026-06-21 Branch B V2 Autosomal Burden Audit

Current report:
`docs/reports/branch_b_v2_autosomal_burden_audit_2026-06-21.md`.

This audit reviewed the remaining autosomal `V2_POSITIVE_SUPPORT_REVIEW` rows
after sex-route refinement.

Remaining autosomal positive-support rows:

- Y1-Y8 truth cohort: 78.
- H1-H6 truth subset: 33.
- H7-H16 context subset: 24.
- 2026-06-15 context: 127.

The dominant unresolved condition is:

```text
UNKNOWN_BACKGROUND + NO_NULL_SUPPORT
```

All audited autosomal positive rows currently have:

- `matched_negative_background_status=UNKNOWN_BACKGROUND`;
- `calibration_null_status=NO_NULL_SUPPORT`;
- `refmap_status=OK`;
- `sample_noise_status=OK`;
- `cnvpro_like_evidence_status=SHADOW_EVIDENCE_ONLY`.

A direction-support split was audited but not implemented as classifier logic.
It is not safe as a hard filter or universal positive-support downgrade because
several locked autosomal truth top candidates would be labeled as A-only /
weak Branch B direction support:

- Y2 chr14 gain;
- Y4 chr13 gain;
- H2 chr6 gain;
- H3 chr13 gain;
- H4 chr15 gain.

H6 chr21 remains Branch B-direction-supported, but that is not enough to make
the direction-support rule generally safe. Direction support can be exposed as
review evidence in a later loop, but it must not hard-suppress or final-demote
autosomal candidates without additional validation.

## 2026-06-21 Branch B V2 Direction-Support Review Label

Current report:
`docs/reports/branch_b_v2_direction_support_label_2026-06-21.md`.

This loop implemented the next safe step from the autosomal burden audit:
Branch B-side direction support is now emitted as review evidence only.

New Branch B V2 classifier output fields:

- `v2_direction_support_label`;
- `v2_direction_support_reason`.

The label contract is explicitly non-final:

- it does not change WisecondorX predict;
- it does not change Branch A candidate discovery;
- it does not change `v2_candidate_class`;
- it does not change `v2_classifier_action`;
- it does not change `v2_final_report_impact`;
- it must not hard-suppress or final-demote a candidate.

Remote materialization under the explicit gap2m overlay produced:

- Y1-Y8: 97 rows; `B_DIRECTION_SUPPORTED=66`,
  `A_ONLY_WEAK_B_DIRECTION=20`, `B_DIRECTION_CONFLICT=11`.
- H1-H16: 105 rows; `B_DIRECTION_SUPPORTED=68`,
  `A_ONLY_WEAK_B_DIRECTION=26`, `B_DIRECTION_CONFLICT=11`.
- 2026-06-15: 165 rows; `B_DIRECTION_SUPPORTED=97`,
  `A_ONLY_WEAK_B_DIRECTION=40`, `B_DIRECTION_CONFLICT=28`.

Autosomal `V2_POSITIVE_SUPPORT_REVIEW` subset:

- Y1-Y8: 78 rows; `B_DIRECTION_SUPPORTED=58`,
  `A_ONLY_WEAK_B_DIRECTION=20`.
- H1-H16: 57 rows; `B_DIRECTION_SUPPORTED=35`,
  `A_ONLY_WEAK_B_DIRECTION=22`.
- 2026-06-15: 127 rows; `B_DIRECTION_SUPPORTED=89`,
  `A_ONLY_WEAK_B_DIRECTION=38`.

Truth preservation after this implementation remains:

- Y1-Y8: truth preserved 10/10, FN=0, hard-suppressed truth=0.
- H1-H16: truth preserved 10/10, FN=0, hard-suppressed truth=0.
- H6 chr21 remains `V2_POSITIVE_SUPPORT_REVIEW` with
  `v2_direction_support_label=B_DIRECTION_SUPPORTED`,
  `v2_final_report_impact=none_shadow_only`, and
  `top_a_abs_zscore=7.113507302991461`.

Remote unit tests:

- `tests/unit/test_branch_b_v2_classifier.py`
- `tests/unit/test_branch_b_v2_benchmark.py`
- `tests/unit/test_branch_ab_phase12_workflow_contract.py`

Result:

```text
32 passed in 0.71s
```

Important interpretation:

- This is an evidence-label improvement, not a final classifier promotion.
- The main unresolved Branch B V2 condition remains
  `UNKNOWN_BACKGROUND + NO_NULL_SUPPORT`.
- The 2026-06-15 cohort remains burden/context only because it has no locked
  truth table.

## 2026-06-21 Branch B V2 Background Context Review Label

Current report:
`docs/reports/branch_b_v2_background_context_label_2026-06-21.md`.

This loop made the remaining background limitation explicit in Branch B V2
classifier outputs:

- `v2_background_context_label`;
- `v2_background_context_reason`;
- `background_context_label_counts` in classifier summaries.

The labels are review context only. They do not change:

- WisecondorX predict;
- Branch A candidate discovery;
- `v2_candidate_class`;
- `v2_classifier_action`;
- `v2_evidence_tier`;
- `v2_direction_support_label`;
- `v2_final_report_impact`;
- Branch S status;
- report-release status.

Remote materialization under the explicit gap2m overlay produced:

- Y1-Y8: 97 rows; `UNKNOWN_BACKGROUND_NO_NULL_SUPPORT=86`,
  `SHADOW_BACKGROUND_NO_NULL_SUPPORT=11`.
- H1-H16: 105 rows; `UNKNOWN_BACKGROUND_NO_NULL_SUPPORT=69`,
  `SHADOW_BACKGROUND_NO_NULL_SUPPORT=36`.
- 2026-06-15: 165 rows; `UNKNOWN_BACKGROUND_NO_NULL_SUPPORT=155`,
  `SHADOW_BACKGROUND_NO_NULL_SUPPORT=10`.

Truth preservation remains unchanged:

- Y1-Y8: truth preserved 10/10, FN=0, hard-suppressed truth=0.
- H1-H16: truth preserved 10/10, FN=0, hard-suppressed truth=0.
- H6 chr21 remains `V2_POSITIVE_SUPPORT_REVIEW`,
  `v2_background_context_label=UNKNOWN_BACKGROUND_NO_NULL_SUPPORT`,
  `v2_direction_support_label=B_DIRECTION_SUPPORTED`,
  `v2_final_report_impact=none_shadow_only`, and
  `a_abs_zscore=7.113507302991461`.

Important interpretation:

- `UNKNOWN_BACKGROUND_NO_NULL_SUPPORT` means background evidence is not
  informative; it must not be interpreted as benign or background-compatible.
- This is not FP reduction yet. It is a safer contract for the next
  truth-safe burden stratification step.
- The 2026-06-15 cohort remains burden/context only because it has no locked
  truth table.

## 2026-06-21 Branch B V2 Method Reset And Threshold Inventory

Current report:
`docs/reports/branch_b_v2_method_reset_threshold_inventory_2026-06-21.md`.

This loop resets the Branch B V2 interpretation from legacy/current-code
Branch B filtering toward Branch-A-anchored evidence stratification.

Key contract:

- Branch A remains the WisecondorX/CBS-derived primary discovery layer.
- Branch B V2 does not create B-only final/report events.
- GC/RC corrected signal is auxiliary context only; it must not define the
  primary gain/loss direction and must not hard-filter Branch A positives.
- The old direction-conflict wording is replaced by B-side signal-context
  wording:
  - `B_SIGNAL_SUPPORTED_A_DIRECTION`;
  - `A_ANCHORED_WEAK_B_SIGNAL`;
  - `B_SIGNAL_DISCORDANT_WITH_A_DIRECTION`;
  - `NO_POSITIVE_A_SIGNAL`.
- `B_SIGNAL_DISCORDANT_WITH_A_DIRECTION` means auxiliary B-side signal is
  discordant with Branch A direction. It does not mean Branch B called an
  opposite event.
- Current first-round V2 dispositions are:
  - `report_candidate`;
  - `review_candidate`;
  - `technical_risk_review`;
  - `background_unknown_review`;
  - `sca_branch_s_review`.

Legacy/current-code thresholds were inventoried and remain excluded from V2
final decision evidence unless a future truth-safe ablation promotes a specific
rule. This includes old `min_event_bins`, `min_abs_calibrated_z`,
`high_confidence_z`, A-branch artifact-protection thresholds, and CNVpro-like
length thresholds. The CNVpro-like size values are review/reportability tiers,
not truth/falsity filters.

Remote validation and materialization:

- Unit tests:
  `tests/unit/test_branch_b_v2_classifier.py`,
  `tests/unit/test_branch_b_v2_benchmark.py`,
  `tests/unit/test_branch_ab_phase12_workflow_contract.py`,
  `tests/unit/test_cnv_report.py`:
  `45 passed in 1.26s`.
- Snakemake dry-runs succeeded for
  `branch_b_v2_benchmark branch_s_review cnv_report` under all three active
  gap2m configs.
- Forced materialization completed for `cnv_branch_b_v2_classifier`,
  `cnv_branch_b_v2_benchmark`, and `branch_b_v2_benchmark` under all three
  active gap2m configs.

Materialized preservation result:

- Y1-Y8: candidates=97; truth preserved 10/10; FN=0;
  hard-suppressed truth=0.
- H1-H16: candidates=105; truth preserved 10/10; FN=0;
  hard-suppressed truth=0.
- 2026-06-15: candidates=165; no locked truth; status `skipped_no_truth`.

H6 chr21 remains retained:

```text
top_candidate_id=H6.A0003
top_a_abs_zscore=7.113507302991461
top_v2_disposition=background_unknown_review
top_v2_b_signal_context_label=B_SIGNAL_SUPPORTED_A_DIRECTION
top_v2_gc_rc_context_label=GC_RC_ATTENUATED
```

Materialized V2 disposition counts:

- Y1-Y8: `background_unknown_review=84`, `sca_branch_s_review=13`.
- H1-H16: `background_unknown_review=63`, `sca_branch_s_review=42`.
- 2026-06-15: `background_unknown_review=151`, `sca_branch_s_review=14`.

This loop improves interpretability and truth-safety. It still does not solve
FP/review burden and does not promote Branch B V2, Branch S, the gap2m overlay,
or the shadow reference to final release.

## 2026-06-21 Branch B V2 Truth-Safe Filter Contract

Current report:
`docs/reports/branch_b_v2_truth_safe_filter_2026-06-21.md`.

This loop adds an explicit filter-action contract on top of the Branch B V2
method reset. It keeps Branch A as the only discovery anchor and keeps Branch B
V2 as candidate-level evidence/disposition logic.

New classifier fields:

- `v2_filter_version`;
- `v2_filter_action`;
- `v2_filter_reason`;
- `v2_filter_scope`;
- `v2_filter_hard_suppression_allowed`.

New benchmark fields:

- `top_v2_filter_action` in truth-overlap metrics;
- `v2_filter_suppressed_count` per sample;
- `v2_filter_suppressed_candidate_count` in benchmark summary JSON.

Current filter contract:

- `suppress_workflow_contract_risk` is the only hard-suppression action.
- The current implemented hard-suppression reason is
  `same_chrom_ref_leakage_contract_risk`.
- GC/RC attenuation, B-side signal discordance, unknown background, length
  tier, low clean support, and high region-risk burden cannot hard-suppress a
  Branch A candidate.
- Low clean support/high-risk burden can only downgrade to
  `technical_risk_review`.
- Sex-chromosome candidates route to Branch S review; this does not finalize
  SCA.

Remote validation:

```text
tests/unit/test_branch_b_v2_classifier.py
tests/unit/test_branch_b_v2_benchmark.py
26 passed in 0.90s
```

Broader remote unit tests after syncing the full contract:

```text
55 passed in 1.37s
```

Remote Snakemake dry-runs succeeded for
`branch_b_v2_benchmark branch_s_review cnv_report` under all three active
gap2m configs. Forced materialization then refreshed
`cnv_branch_b_v2_classifier`, `cnv_branch_b_v2_benchmark`, and
`branch_b_v2_benchmark`.

Materialized result:

- Y1-Y8: candidates=97; truth preserved 10/10; FN=0;
  hard-suppressed truth=0; filter-suppressed candidates=0.
- H1-H16: candidates=105; truth preserved 10/10; FN=0;
  hard-suppressed truth=0; filter-suppressed candidates=0.
- 2026-06-15: candidates=165; no locked truth; status `skipped_no_truth`;
  filter-suppressed candidates=0.

Remote output files now include `v2_filter_action`, `top_v2_filter_action`,
and `v2_filter_suppressed_count`.

This is a materialized truth-safe filter contract validation. It does not yet
claim FP reduction, Branch B V2 final promotion, Branch S finalization, report
release, or production reference promotion.

## 2026-06-21 P1-P6 Credibility Audit Override

Current audit report:
`docs/reports/p1_p6_result_credibility_audit_2026-06-21.md`.

This audit reclassifies the current G1/P1-P6 outputs by evidentiary value:

- G1/P1 reference audit is valid current evidence for R0/R1/R2 reference-rebuild
  eligibility only. It does not promote `h_r0_shadow_ref_20260619` to
  production and does not define N0/N1/N2.
- P2 Branch A validation is the current performance-relevant evidence under the
  shadow reference: Y1-Y8 truth detected 10/10, H1-H16 truth detected 10/10,
  and 2026-06-15 has no truth table.
- P3 Branch B evidence ledger is a valid engineering artifact only. It provides
  one review-safe row per Branch A candidate with `report_impact=none_shadow_only`;
  it is not Branch B V2 performance validation.
- Current-code Branch B fields such as `final_disposition`,
  `branch_b_keep_event`, and artifact status are historical/shadow context only
  and must not be used as Branch B V2 evidence, N0/N1, benchmark truth, or
  report-release evidence.
- Current Branch B V2 classifier shadow output is not valid V2-only performance
  evidence because it still reads legacy `final_disposition` and produces
  legacy-derived shadow classes.
- P5 Branch S is a valid engineering/report-boundary artifact only. It remains
  `review_development_only` and is not final SCA.
- P6 2026-06-15 report package is development-only:
  `report_contract.status=development_only_not_final_release`.

Remote audit notes from 2026-06-21:

- Result-file inspection on `ssh fengxian` succeeded.
- The requested Snakemake dry-run for
  `branch_b_evidence branch_s_review cnv_report` failed at remote Snakefile
  parse time with an `IndentationError`. The first 40 lines of the remote
  Snakefile looked consistent with local tracked `Snakefile`, and local
  `git diff -- Snakefile` was empty.
- Therefore current materialized results can be used for credibility audit, but
  no new workflow-validation claim should be made until the remote mirror or
  Snakemake parse issue is fixed.

## 2026-06-20 Current State Override

Current source of truth for phase order:
`docs/reports/branch_ab_v2_gate_phase_plan_2026-06-20.md`.

Current constraints-only document:
`docs/reports/branch_ab_v2_rnd_constraints_2026-06-18.md`.

Older Branch B V2 / N1 shadow-background sections below remain historical
context only. They do not define the current report route, production reference
route, or N0/N1/N2 promotion state.

Current G1/P1 state:

- H7-H16 reference audit has been materialized remotely under
  `h_r0_shadow_ref_20260619`.
- Current H R labels are `R0=H9,H10,H11,H12,H15,H16`, `R1=H7,H8,H13,H14`,
  `R2=none`.
- The stable decision artifact is
  `docs/reports/h7_h16_reference_cohort_decision_2026-06-20.md` plus TSV.
- `R0/R1/R2` are reference-rebuild labels only and must not be reused as
  `N0/N1/N2`.
- The H-augmented reference remains shadow-only. P2 must rerun WisecondorX
  predict and Branch A under the same named reference/config before any Branch B
  or report-facing promotion.

Current P2 state:

- A formal `branch_a_validation` workflow target now exists for Branch A no-FN
  validation without consuming Branch B, Branch S, or report outputs.
- Under `h_r0_shadow_ref_20260619`, remote P2 materialization produced:
  - Y1-Y8: 10/10 truth detected, FN=0, Branch A candidates=131.
  - H1-H16: 10/10 truth detected, FN=0, Branch A candidates=221.
  - H6 chr21 status: detected.
  - 2026-06-15: no truth table, Branch A candidates=201, exploratory burden
    only.
- P2 supports Branch A sensitivity under the current shadow reference. It does
  not promote the H reference to production and does not make Branch B, Branch S,
  or 2026-06-15 reports final-reportable.
- Detailed report:
  `docs/reports/branch_a_validation_h_r0_shadow_2026-06-20.md`.

Current P3 state:

- A formal `branch_b_evidence` workflow target now exists for P3 Branch B
  candidate-level evidence materialization without requesting Branch B V2
  classifier, Branch S, matched-negative percentile, or report outputs.
- Under `h_r0_shadow_ref_20260619`, remote P3 evidence ledger materialization
  produced one row per Branch A candidate:
  - Y1-Y8: 8 ledger files, 131 candidate rows.
  - H1-H16: 16 ledger files, 221 candidate rows.
  - 2026-06-15: 5 ledger files, 201 candidate rows.
- The P3 schema includes explicit review-safe fields:
  `branch_b_direction_support`, `copy_number_like_amplitude`,
  `mosaic_proxy`, `loh_evidence`, `upd_evidence`, `background_source`,
  `background_status`, `region_risk_context`, `sample_noise_context`,
  `cnvpro_consistency_status`, `disposition`, `disposition_reason`, and
  `report_impact`.
- Current P3 output is evidence/refinement only. All rows are
  `REVIEW_REQUIRED` with `report_impact=none_shadow_only`; this does not
  promote Branch B, suppress Branch A positives, or make 2026-06-15 reports
  releasable.
- Detailed report:
  `docs/reports/branch_b_p3_evidence_contract_2026-06-20.md`.

Current P5 state:

- A formal `branch_s_review` workflow target now exists for Branch S/SCA
  report-boundary materialization without requesting Branch B V2 classifier,
  matched-negative percentile, or report outputs.
- Under `h_r0_shadow_ref_20260619`, remote P5 Branch S materialization produced
  report-boundary summaries:
  - Y1-Y8: 8 evidence tables, 8 score tables, 8 summary JSON files.
  - H1-H16: 16 evidence tables, 16 score tables, 16 summary JSON files.
  - 2026-06-15: 5 evidence tables, 5 score tables, 5 summary JSON files.
- The P5 summary schema includes explicit report-boundary fields:
  `expected_x_ploidy`, `expected_y_ploidy`, `x_nonpar_direction`,
  `x_par_context`, `y_nonpar_direction`, `y_par_or_homology_context`,
  `branch_a_x_support`, `branch_a_y_support`, `sca_candidate_state`,
  `sca_confidence_tier`, `sca_output_mode`, `sca_uncertainty_reason`, and
  `report_text_status`.
- Current Branch S remains review/development-only:
  `sca_output_mode=review_development_only` for all materialized summaries.
  This satisfies the G5-review boundary only; it does not satisfy G5-final
  promotion because locked SCA truth coverage remains incomplete.
- Detailed report:
  `docs/reports/branch_s_p5_report_boundary_2026-06-20.md`.

Current P6 state:

- The `cnv_report` workflow now carries P3/P5 review context through
  summary-only inputs:
  - `CNV_B_EVIDENCE_SUMMARY`;
  - `CNV_BRANCH_S_SUMMARY`.
- The report workflow still does not consume P3 raw evidence ledger TSV,
  Branch S raw evidence TSV, Branch S score TSV, matched-negative percentile
  TSV, or Branch B V2 classifier TSV.
- The 2026-06-15 P6 report package was materialized remotely under
  `h_r0_shadow_ref_20260619`:
  - `wisecondorx/cnv/report/cnv_summary.tsv`;
  - `wisecondorx/cnv/report/cnv_summary.json`;
  - `wisecondorx/cnv/report/cnv_summary.md`;
  - `wisecondorx/cnv/report/cnv_summary.html`.
- The report JSON has
  `report_contract.status=development_only_not_final_release`,
  `branch_b_evidence_summary_count=5`, and `branch_s_summary_count=5`.
- All 2026-06-15 report rows explicitly show
  `UNKNOWN_BACKGROUND` P3 background and
  `review_development_only` Branch S/SCA status.
- P6 is a development/report-boundary package only. It does not release the
  2026-06-15 samples as final reports and does not promote Branch B or Branch S.
- Detailed report:
  `docs/reports/p6_report_package_contract_2026-06-20.md`.

> 更新日期：2026-06-19

## 文档定位
本文件只记录当前已确认的项目事实、有效入口、有效验证结论与当前阻塞点，不承担计划文档职责。

## 当前阶段
当前项目处于：

`Branch B 收敛 + re-sequencing cohort 回流设计阶段`

补充说明：
- 项目算法主线仍然是 Branch B 结果收敛。
- 本轮代码结构审查、runtime tracking 修复、report 输入闭环修复和 modular dry-run，属于中间代码质控与工程稳定性工作。
- 这些工作与项目算法本身没有直接等价关系，不应被误写成“项目主阶段已切换”。
- 当前新增的生产约束是：FAIL 样本需要 re-sequencing 后重新进入 QC 与 reference cohort 候选路径，不能继续作为口头状态或手工补丁处理。

## 已确认的代码与仓库状态
1. 模块化代码已进入 `main`。当前 tracked `main` 已同步到 `origin/main`；本地仍存在若干未跟踪的历史脚本/样本表，不应在无明确需求时纳入提交。
2. 当前有效的主入口仍是根目录下：
   - `Snakefile`
   - `build_samples.yaml`
   - `config_qc.yaml`
   - `config_reference.yaml`
   - `predict_samples.yaml`
   - `config_predict.yaml`
3. 当前规则文件已按模块拆分在 `rules/` 下，实际落地的关键文件包括：
   - `common_preprocess.smk`
   - `pipeline_modes.smk`
   - `reference_layout.smk`
   - `reference_workflow.smk`
   - `reference_assets.smk`
   - `predict_layout.smk`
   - `predict_workflow.smk`
   - `runtime_layout.smk`
   - `runtime_tracking.smk`
   - `qc_workflow.smk`
   - `qc_framework.smk`
   - `target_assembly.smk`
   - `script_entrypoints.smk`
4. 当前 Python 实现已按包沉淀到 `pgta/` 下，核心子包为：
   - `pgta/core`
   - `pgta/qc`
   - `pgta/reference`
   - `pgta/predict`
   - `pgta/validation`
5. 顶层 `scripts/` 当前主要承担 CLI / 初始化 / 兼容入口，不再是全部分析逻辑的唯一承载位置。

## 本轮已确认的新事实
1. `rules/runtime_tracking.smk` 已修复为跟踪按 `ref_set` 生成的真实 benchmark 与输出路径，不再继续记录旧的非 `ref_set` reference benchmark 路径。
2. `rules/predict_workflow.smk` 中 `cnv_report_summary` 已改为只在 summary 真正属于本次 rule 输入时才传：
   - `--evaluation-summary`
   - `--ml-summary`
   - `--benchmark-summary`
3. 上述修复已提交并推送：
   - `e784187` `fix: align runtime tracking and report inputs`
4. 上述修复属于流程工程侧质控，不构成 Branch B 算法主线变更。
5. FAIL 样本 re-sequencing 后重新加入建模样本已有首轮正式 manifest / cohort 决策入口：
   - `build_reference.resequencing.manifest`
   - `build_reference.resequencing.allowed_statuses_for_reference`
   - 默认只有 `promoted` 状态可进入 reference cohort selection
   - `replacement_policy=replace_original` 会替换原始 `source_sample_id`
6. `build_reference.groups` 现在会在 Snakemake 解析阶段 fail fast，如果引用了不存在于 `samples` 的样本。
7. `rules/predict_workflow.smk` 中 Branch B 后处理 rule 重复定义已收敛为单一实现。
8. `pgta/predict/evaluation.py` 与 `pgta/predict/benchmark.py` 的 truth/event 匹配公共逻辑已收敛到 `pgta/predict/truth.py`。

## 当前已确认的 remote 验证结论
以下结论都来自 remote `fengxian`，不是本地替代验证：

1. 本轮 remote dry-run 使用的是：
   - 代码树：`/data/project/CNV/PGT-A/refactor_validation_20260419`
   - 配置/运行目录：`/data/project/CNV/PGT-A/refactor_validation_run_20260419`
2. `reference_qc reference` dry-run 已通过，DAG 已正确展开到按 `ref_set` 的 benchmark 路径：
   - `monitor/benchmarks/reference/select_cohort/pass_warn.tsv`
   - `reference_prefilter/pass_warn/XX.tsv`
   - `reference_prefilter/pass_warn/XY.tsv`
   - `tune_wisecondorx_reference_qc/pass_warn.tsv`
   - `build_reference/pass_warn/XX.tsv`
   - `build_reference/pass_warn/XY.tsv`
   - `build_gender_reference/pass_warn.tsv`
3. `cnv_qc cnv` dry-run 已通过；当前结果为 `Nothing to be done`。
4. `cnv_report` 在两种配置口径下都已做 remote dry-run：
   - 现有远端 predict 配置包含 `cnv_eval / cnv_benchmark / cnv_report`，因此 summary 输入存在，符合配置事实。
   - 额外构造的临时配置只保留 `mapping / metadata / cnv_qc / cnv / cnv_report` 时，`cnv_report_summary` 的输入列表中已不再出现 `evaluation_summary / benchmark_summary`，说明本轮修复已生效。
5. 本轮新增的 re-sequencing manifest 契约已通过 remote 单元测试与 Snakemake dry-run 验证：
   - promoted re-sequenced sample can replace source sample
   - candidate re-sequenced sample is not selected until promoted
   - manifest appears as `select_reference_cohort` input when configured
   - missing `build_reference.groups` sample raises a parse-time error
6. 本轮 remote 验证是单元测试与 dry-run，不是完整真实执行；不能把它写成“流程已完整跑通”。

## 当前仍有效的较早 remote 结果基线
以下较早结论仍来自 remote 结果文件，且本轮没有证据推翻它们：

1. 2026-04-23 的 predict 基线结果仍显示：
   - `truth_recall = 1.0`
   - `truth_detected_count = 10`
   - `branch_b_detection_rate = 1.0`
2. `sex-aware ranking / report suppression` 的结果口径仍保持有效。
3. 本轮没有重新做完整 predict 真实运行，因此这些结论本轮是“沿用既有 remote 结果”，不是“本轮重新实跑确认”。

## 当前有效路径与文件
1. 本地仓库路径：
   - `D:\Pipeline\PGT-A`
2. 当前 remote 验证代码树：
   - `/data/project/CNV/PGT-A/refactor_validation_20260419`
3. 当前 remote 验证运行目录：
   - `/data/project/CNV/PGT-A/refactor_validation_run_20260419`
4. 本轮直接相关文件：
   - `rules/runtime_tracking.smk`
   - `rules/predict_workflow.smk`
   - `tests/server_validation/03_snakemake_dryrun.sh`
   - `config_predict.yaml`
   - `README.md`

## 当前阻塞点与风险
1. remote `/data/project/CNV/PGT-A` 根目录不是 Git 工作树；remote 验证当前依赖已有代码副本目录，验证过程需要显式同步和恢复文件。
2. 当前 remote 验证入口仍分散在：
   - 代码树目录
   - 运行结果目录
   - 临时配置副本
   这增加了“代码版本与运行目录脱钩”的风险。
3. 本轮只完成了 dry-run 级别验证，尚未完成“基于当前 `main` 的完整 remote 真实运行验证”。
4. 当前本地默认 `config_predict.yaml` 仅请求：
   - `mapping`
   - `metadata`
   - `cnv_qc`
   - `cnv`
   而 remote 验证配置还显式请求了：
   - `cnv_eval`
   - `cnv_benchmark`
   - `cnv_report`
   后续需要明确哪套配置口径作为持续验证默认口径。
5. re-sequencing 样本回流已有版本化 manifest 入口，但真实 FAIL 样本重测数据尚未在本轮完成非 dry-run 验证。
6. 默认 reference 发布路径会复制单一 `DEFAULT_REF_SET` 到通用输出；后续比较 re-sequencing cohort 时，需要区分 cohort 比较产物与发布产物。

## 当前不应误判为已完成的事项
1. 不能因为 dry-run 已通过，就宣称模块化 workflow 已完成最终 remote 验证。
2. 不能把 remote 临时同步补丁目录的 dry-run，等同于“remote Git checkout 已统一”。
3. 不能把较早的 Branch B full-result 基线，写成“本轮重新实跑完成”。
4. 不能把 FAIL 样本 re-sequencing 需求写成“已重新纳入建模”，除非已有真实 manifest、QC 结果、cohort 决策和 remote 非 dry-run 验证证据。
5. 不能通过新增临时脚本或平行 workflow 来处理 re-sequencing 回流；应优先收敛到正式 config、manifest、rule 或 `pgta/` 模块内。

## 2026-06-18 已确认研发事实

1. `docs/branch_b_independent_detector_methodology_2026-04-22.md` 中“Branch B 独立 detector”方向已不再作为当前研发主线。
2. 当前主线以 `docs/reports/branch_ab_v2_rnd_constraints_2026-06-18.md` 为准：
   - WisecondorX-first 保持不变；
   - Branch A 是高灵敏度候选发现层；
   - Branch B V2 应收缩为 candidate-level evidence/classification shadow layer；
   - legacy Branch B 只作为对照，不继续增加补丁式 hard filter。
3. 2026-06-17 的 refmap/calibration 风险仍有效：
   - 当前 0615 run 中 event-level refmap count 证据缺失；
   - `calibration_null_eligible` 全为 0 的情况下，Branch B calibration 会 fallback 到宽背景；
   - 这会使强 Branch A 候选被弱 Branch B calibrated evidence 过度过滤。
4. 远端 `results_h_20260608_mask_only/wisecondorx/cnv/report/cnv_summary.tsv` 显示：
   - H7-H16 均为 QC PASS/XY；
   - H8/H13/H14 有 Branch B kept review event；
   - H7/H9/H10/H11/H12/H15/H16 虽然 Branch B kept=0，但仍存在强 A-sensitive signal。
5. 因此 H7-H16 当前不能直接作为 N0 high-confidence negative；但这只限制 Branch B matched-negative calibration，不等于禁止 reference-rebuild exploration。
6. 当前可写入计划的 reference 候选状态需要拆分为两套标签：
   - `R0/R1/R2`：reference-rebuild eligibility；
   - `N0/N1/N2`：Branch B empirical-null / matched-negative calibration。
7. 旧 32-ref 下的 H7-H16 Branch A signal 可能代表 reference mismatch、batch shift、高重复区域不稳定或真实个体异常，不能单独作为 R0/R1 排除标准。
8. 当前建议先建立 H7-H16 的 `R0/R1/R2` 表，并构建 named shadow reference variant，而不是直接 promote 到 production reference。
9. 若构建新的 reference，但 binsize/mask/preprocess contract 不变，已有 predict NPZ 可复用；但 WisecondorX predict、Branch A candidate、Branch B、evaluation、benchmark、report 必须重跑。
10. 若 binsize、mask 或 preprocessing 发生变化，reference 和 predict 样本 NPZ 都必须重建。
11. Branch A/B V2 Phase 1 evidence ledger 已落地为 shadow-only 输出：
   - `wisecondorx/cnv/postprocess/evidence_ledger/{sample}.candidate_evidence.tsv`
   - 不改变当前 final report。
12. Branch A/B V2 Phase 2 negative bank 已落地：
   - seed 文件：`docs/reports/branch_ab_v2_negative_bank_seed_2026-06-18.tsv`
   - 当前 `N0=0`、`N1=6`、`N2=4`。
13. Branch A/B V2 Phase 3 matched-negative percentile shadow path 已落地：
   - `wisecondorx/cnv/postprocess/matched_negative/{sample}.candidate_evidence.tsv`
   - 当前因为无 N0，Y1-Y8 结果全部为 `UNKNOWN_BACKGROUND` / `REVIEW_NO_CALL`
   - 该输出不进入 `cnv_report`，也不作为 hard artifact filter。
14. Branch A/B V2 negative-bank readiness contract 已落地：
   - `negative_bank_summary.json` 现在输出 `matched_negative_ready`、`matched_negative_blocking_reason`、`n0_sample_ids`、`n1_shadow_reference_candidate_count`、`n2_holdout_count`；
   - 当前远端结果为 `matched_negative_ready=false`；
   - 当前阻塞原因是 `no_n0_locked_clean_negative_samples`；
   - 详细报告：`docs/reports/branch_ab_v2_negative_bank_readiness_2026-06-19.md`。

## 2026-06-18 仍待验证事项

1. Phase 3 目前只有安全 fallback 结果；当前 summary 已机器可读地标记为 `matched_negative_ready=false`。这是旧 32-ref / 当前 seed labels 下的状态，不是新 reference rebuild 后的最终 Phase 3 结论。
2. `matched_negative_ready=false` 只阻塞 Branch B calibration，不阻塞 reference-rebuild shadow experiment。
3. H7-H16 需要重新按 `R0/R1/R2` reference-rebuild eligibility 分层；该分层不能只依赖旧 32-ref 下的 Branch A signal。
4. 若 H7-H16 或其子集进入 named shadow reference，必须在新 ref 后重跑 WisecondorX predict、Branch A、Branch B evidence、negative-bank labeling、evaluation、benchmark 和 report，并生成新的 negative-bank version。
5. H7-H16 的 named shadow reference 尚未完成正式 Snakemake workflow 验证。
6. Branch B V2 final-report promotion 仍需 post-rebuild FN=0、H6 chr21 跟踪、FP/review burden 评估、unknown background 安全行为和 rollback path；当前只能作为 shadow evidence。
7. Branch S 仍只是研发方向，尚未替换或改变当前 sex calling / SCA report 逻辑；它需要 locked SCA truth 和 clean XX/XY negatives，不能因 reference rebuild 完成而自动 promotion。
8. 任何基于 H7-H16 的阈值、artifact rule 或 matched-negative null，都必须先经过 ablation，不能直接写入 production 规则。

## 2026-06-19 Branch S / SCA 验证状态

1. Branch S 当前仍为 shadow-only：
   - 不替换 WisecondorX predict；
   - 不替换 Branch A candidate；
   - 不替换 legacy Branch B final events；
   - 不替换当前 sex calling；
   - 不进入 `cnv_report`。
2. 远端当前已 materialized 的 Branch S 输出覆盖：
   - Y1-Y8：`results_build_ref_v2_mask_only/wisecondorx/cnv/postprocess/branch_s`，8 个 evidence / 8 个 score / 8 个 summary；
   - H1-H16：`results_h_20260608_mask_only/wisecondorx/cnv/postprocess/branch_s`，16 个 evidence / 16 个 score / 16 个 summary；
   - 2026-06-15 五样本：`results_20260615_mask_only/wisecondorx/cnv/postprocess/branch_s`，5 个 evidence / 5 个 score / 5 个 summary；
   - summary 均保持 `final_report_impact=none_shadow_only`、`replaces_final_report=false`。
3. Branch S 已完成正式结果目录物化，但还不能用于 SCA 性能 promotion；它只能作为 shadow evidence 和 review 设计依据。
4. 当前 SCA truth 覆盖不足：
   - 已有 Y3 `chrX loss / 45,XO`；
   - H5/H6 有 `chrX loss` truth，且 Branch S 输出已 materialize，当前 score direction 与 Branch A chrX-loss evidence 一致；
   - 缺少 X gain、XXY、XXX、XYY、Y loss、mosaic SCA fraction series 和跨批次 clean XX/XY negative。
5. 新增研发文档：
   - `docs/reports/branch_s_sca_validation_plan_2026-06-19.md`
   - `docs/reports/branch_s_materialization_2026-06-19.md`
6. 当前判断：
   - Branch S 可以用于 SCA review 设计；
   - 不能生成 `PASS` / `CONFIRMED` SCA 报告结论；
   - 不能作为 2026-06-15 五个样本报告放行或过滤依据。
## 2026-06-19 Branch A/B Phase 1B asset decision matrix

1. Phase 1B has been consolidated as an asset-decision document:
   - `docs/reports/branch_ab_v2_phase1b_asset_decision_matrix_2026-06-19.md`
2. Current production-path asset evidence was refreshed from remote:
   - `/data/project/CNV/PGT-A/refactor_validation_20260419/results_build_ref_v2_mask_only/reference/assets`
3. Current asset contract remains unchanged:
   - WisecondorX-first;
   - Branch A remains WisecondorX predict/CBS candidate discovery;
   - Branch B V2 asset evidence remains shadow/candidate-level unless locked ablation proves no recall regression;
   - Branch S remains shadow-only.
4. Current Phase 1B decision:
   - `annotation_only` remains the default for PAR, sex-homology, low-map, repeat/segdup, CNVseq unique-bin, and mappability-style evidence.
   - `reference_build_mask` is allowed only as a named shadow reference experiment.
   - `full_remap_reference` using `hg19.mapability.fa.gz`, `maskpar`, or another altered FASTA is not justified now.
5. Current remote mask asset counts:
   - atomic 50 kb bins: `61,927`;
   - analysis 100 kb bins: `30,970`;
   - QC 1 Mb bins: `3,113`;
   - combined mask labels: `pass=81,406`, `hard=7,892`, `soft=6,710`, `dynamic=2`.
6. Current hard-mask BED contains `4,881` atomic BED3 intervals and is still a WisecondorX predict blacklist asset.
7. Current analysis-level PAR and sex-homology annotations remain annotation-only:
   - PAR bins: `31`;
   - sex-homology bins: `31`;
   - XTR bins: `0`.
8. Repeat/segmental-duplication/low-mappability burden must not be promoted into a hard event-level filter without locked ablation, because previous evidence showed known positives and FP-like events can share similar soft-mask burden.
