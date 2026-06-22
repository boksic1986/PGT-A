---
status: active_current_index
decision_use: active_context_entrypoint
last_updated: 2026-06-22
active_reference_id: h_r0_shadow_ref_20260619
---

# Current Context Index

This file is the first context-restoration entrypoint for the current PGT-A
CNV-seq R&D cycle. It does not replace result files, configs, or workflow code.
It pins the current decision hierarchy so future work does not resume from stale
handoffs or legacy Branch B outputs.

## Required Read Order

1. `docs/CURRENT_CONTEXT_INDEX.md`
2. `docs/handoff/2026-06-22_1618_branch_a_ref_z_plot_sex_chrom_cn_trend_handoff.md`
3. `AGENTS.md`
4. `skills/conversation_handoff/SKILL.md`
5. `skills/pgta_reference_modeling_analysis/SKILL.md`
6. Current result files and configs for the active task

## Active Inputs

active_handoff: docs/handoff/2026-06-22_1618_branch_a_ref_z_plot_sex_chrom_cn_trend_handoff.md
previous_handoff: docs/handoff/2026-06-22_1535_sex_aware_cn_z_plot_branch_s_overlay_handoff.md
active_reference_id: h_r0_shadow_ref_20260619
reference_status: fixed_shadow_baseline_not_production
remote_snakemake_parse_status: repaired_lf_normalized_2026-06-21
branch_a_status: burden_phase1_gap2m_materialized_default_unchanged
branch_b_status: v2_report_visibility_materialized_development_only
branch_s_status: sex_aware_branch_s_overlay_not_final
report_status: branch_a_ref_z_and_sex_aware_cn_plots_materialized_development_only

## Current Evidence Files

- `docs/reports/p1_p6_result_credibility_audit_2026-06-21.md`
- `docs/reports/branch_a_ref_z_plot_sex_chrom_cn_trend_2026-06-22.md`
- `docs/handoff/2026-06-22_1618_branch_a_ref_z_plot_sex_chrom_cn_trend_handoff.md`
- `docs/reports/sex_aware_cn_z_plot_branch_s_overlay_2026-06-22.md`
- `docs/handoff/2026-06-22_1535_sex_aware_cn_z_plot_branch_s_overlay_handoff.md`
- `docs/reports/copy_number_centering_and_sca_fix_2026-06-22.md`
- `docs/handoff/2026-06-22_1445_cn_centering_branch_s_fix_handoff.md`
- `docs/reports/copy_number_ratio_cnv_plot_2026-06-22.md`
- `docs/handoff/2026-06-22_1425_copy_number_ratio_cnv_plot_handoff.md`
- `docs/reports/0615_high_confidence_report_candidates_2026-06-22.md`
- `docs/handoff/2026-06-22_1335_copy_number_cnv_plot_cn_threshold_scatter_handoff.md`
- `docs/handoff/2026-06-22_1305_copy_number_cnv_plot_bin_z_scatter_handoff.md`
- `docs/handoff/2026-06-22_1234_copy_number_cnv_plot_centromere_scatter_handoff.md`
- `docs/handoff/2026-06-22_1155_copy_number_cnv_plot_scatter_background_handoff.md`
- `docs/handoff/2026-06-22_1130_copy_number_cnv_plot_background_fix_handoff.md`
- `docs/handoff/2026-06-22_1100_copy_number_cnv_plot_v2_handoff.md`
- `docs/handoff/2026-06-22_0907_0615_high_confidence_report_handoff.md`
- `docs/handoff/2026-06-22_1025_copy_number_cnv_plot_handoff.md`
- `docs/handoff/2026-06-22_0941_cnv_plot_wisecondor_style_handoff.md`
- `docs/reports/report_main_convergence_cnv_plot_2026-06-22.md`
- `docs/handoff/2026-06-22_0437_report_main_cnv_plot_handoff.md`
- `docs/reports/h7_h16_reference_cohort_decision_2026-06-20.md`
- `docs/reports/branch_a_validation_h_r0_shadow_2026-06-20.md`
- `docs/reports/branch_a_burden_optimization_phase1_2026-06-21.md`
- `docs/reports/branch_b_v2_gap2m_benchmark_2026-06-21.md`
- `docs/reports/branch_b_v2_sex_route_refinement_2026-06-21.md`
- `docs/reports/branch_b_v2_autosomal_burden_audit_2026-06-21.md`
- `docs/reports/branch_b_v2_direction_support_label_2026-06-21.md`
- `docs/reports/branch_b_v2_background_context_label_2026-06-21.md`
- `docs/reports/branch_b_v2_method_reset_threshold_inventory_2026-06-21.md`
- `docs/reports/branch_b_v2_truth_safe_filter_2026-06-21.md`
- `docs/reports/branch_b_v2_burden_stratification_2026-06-21.md`
- `docs/reports/branch_b_v2_report_contract_2026-06-21.md`
- `docs/reports/branch_b_v2_report_layer_filter_2026-06-21.md`
- `docs/reports/branch_b_v2_g1_g8_validation_2026-06-21.md`
- `docs/reports/branch_b_v2_pass2_correction_2026-06-21.md`
- `docs/reports/branch_b_v2_lowres_ref_evidence_2026-06-21.md`
- `docs/reports/branch_b_s_lowres_integration_2026-06-22.md`
- `docs/handoff/2026-06-22_0930_lowres_branch_bs_integration_handoff.md`
- `docs/reports/branch_s_sca_v2_sex_aware_review_2026-06-22.md`
- `docs/handoff/2026-06-22_0101_branch_s_sca_v2_handoff.md`
- `docs/reports/branch_b_v2_report_event_audit_2026-06-21.md`
- `docs/handoff/2026-06-21_2253_branch_b_v2_lowres_ref_evidence_handoff.md`
- `docs/handoff/2026-06-21_2123_branch_b_v2_pass2_correction_handoff.md`
- `docs/handoff/2026-06-21_2049_branch_b_v2_report_burden_pass2_handoff.md`
- `docs/handoff/2026-06-21_1810_g1_g8_current_scheme_validation_handoff.md`
- `docs/handoff/2026-06-21_1730_branch_b_v2_report_layer_filter_handoff.md`
- `docs/reports/branch_b_v2_reference_background_and_sca_design_2026-06-20.md`
- `docs/reports/branch_s_p5_report_boundary_2026-06-20.md`
- `docs/reports/p6_report_package_contract_2026-06-20.md`
- `docs/reports/branch_ab_v2_gate_phase_plan_2026-06-20.md`
- `docs/reports/branch_ab_v2_rnd_constraints_2026-06-18.md`

## Current Module State

### Branch A Ref-Z Plot And Sex-Chrom CN Trend

The active plot contract is now documented in
`docs/reports/branch_a_ref_z_plot_sex_chrom_cn_trend_2026-06-22.md`.

This loop supersedes the z-plot limitation from
`docs/handoff/2026-06-22_1535_sex_aware_cn_z_plot_branch_s_overlay_handoff.md`.
The previous combined plot used Branch B residual `calibrated_z` as the z-axis.
The active z plot now uses per-bin `branch_a_ref_z`:

`branch_a_ref_z = (normalized_signal - ref_bin_stability.ref_median) / max(1.4826 * ref_mad, min_ref_scale)`

The z trend line for report events uses direction quantiles, not median:

- gain/dup: `Q75(branch_a_ref_z)`;
- loss/del: `Q25(branch_a_ref_z)`.

This preserves normal points inside tolerated 2Mb merge gaps while making the
event-level Branch A direction visible. `plot_bins.tsv` now carries
`branch_a_ref_z`, `residual_calibrated_z`, `z_source`, `ref_z_scale`, and
`ref_z_scale_source`.

Branch S `sca_report_review_event` intervals still overlay chrX/chrY. The CN
plot now additionally draws sex-aware Branch S CN trend lines only when the
sex-chromosome bins are interpretable and the median CN deviates from expected
CN by the visual threshold. G8 chrY remains
`sex_chrom_ref_ratio_not_interpretable`.

Remote validation:

- pytest: `58 passed in 2.06s`;
- G1-G8 dry-run and 0615 dry-run parsed successfully without mapping/reference
  rebuild;
- G1-G8 materialization completed: `20 of 20 steps (100%) done`;
- 0615 materialization completed: `14 of 14 steps (100%) done`.

Local synchronized plot paths:

- `D:\Pipeline\PGT-A\reports\g1_g8_cnv_plots\`
- `D:\Pipeline\PGT-A\reports\0615_cnv_plots\`

This remains a visualization/report-support update only. It does not change
Branch A calling, Branch B V2 filtering, Branch S classifier, reference, or
report-event classification.

### 0615 High-Confidence Review Candidates

The latest completed loop is
`docs/reports/0615_high_confidence_report_candidates_2026-06-22.md`.

This is a read-only review of the current 2026-06-15 materialized report output.
It does not change workflow code, thresholds, reference, Branch A, Branch B V2,
Branch S, or remote results.

Current 0615 input check:

- `report_events.tsv`: 71 autosomal report rows.
- plot-bin TSVs: 5/5 samples present.
- Branch S remains development-only context.

Current conservative high-confidence autosomal review candidates:

- `JZ26125845-60-60`: 10 rows.
- `JZ26125843-56-56`, `JZ26125844-59-59`,
  `JZ26125846-61-61`, and `JZ26125847-62-62`: 0 rows under the conservative
  high-confidence screen.

The current 5/5 batch-shared region `chr4:67.50-101.25Mb gain` is excluded
from high-confidence interpretation. This cohort still has no locked truth, so
no TP/FP/FN conclusion is allowed.

### Report Main And CNV Plot

The active CN plot and Branch S correction report is now
`docs/reports/copy_number_centering_and_sca_fix_2026-06-22.md`.

The previous `calibrated_z_mosaic30_cn_proxy` contract is superseded.
`calibrated_z` remains the only signal for the z plot, but it is no longer used
to estimate copy number. The CN plot now derives bin-level CN from the
sample/reference depth ratio and then recenters each sample by the autosomal
non-gap median:

`raw_log2R = log2((sample_cpm + 0.001) / (ref_cpm + 0.001))`

`centered_log2R = raw_log2R - median(raw_log2R over non-gap autosomal bins)`

`CN = 2 * 2^centered_log2R`

where `sample_cpm = expm1(normalized_signal)` and
`ref_cpm = expm1(ref_bin_stability.ref_median)`. `plot_bins_cn.tsv` carries
`raw_log2r`, centered `log2r`, `copy_number_centering_log2_shift`, and
`copy_number_source=normalized_signal_ref_median_log2r_autosome_centered`.
CN values are not clipped in `plot_bins_cn.tsv`; out-of-range points are clipped
only for SVG display and marked in the SVG. `plot_event_support.tsv` records per-report-event
`valid_bin_count`, `cn_support_bin_count`, `cn_same_direction_fraction`,
`median_bin_cn`, `median_log2r`, `median_calibrated_z`, and
`cn_direction_consistency_status`.

G1-G8 and 2026-06-15 have been re-materialized with this centered ratio-derived
CN plot. All materialized samples have autosomal median CN = 2.000 after
centering.
0615 remains no-truth review/context only. G1-G8 remains the first truth-backed
visual check: locked truth is 10/10 preserved, FN=0, and G2 chr8 remains
internal-review visible rather than filtered.

Branch S correction in the same loop:

- segment-level corroboration uses median/robust-median only, not mean;
- G3 and G5 X-loss are now `sca_report_review_event`;
- G2/G6/G8 XY X-gain-like Branch A signals are no-call/sex-consistent audit,
  not strong SCA review.

Local synchronized plot paths:

- `D:\Pipeline\PGT-A\reports\g1_g8_cnv_plots\`
- `D:\Pipeline\PGT-A\reports\0615_cnv_plots\`

### Superseded CN Proxy Plot

The latest completed loop is
`docs/reports/report_main_convergence_cnv_plot_2026-06-22.md`.

`cnv_report` now uses Branch B V2 benchmark report events as the final autosomal
main-table event source when V2 benchmark is available. The final autosomal
report event layer is split into `report_strong_event` and
`report_weak_event`; internal review, filtered audit rows, and Branch S rows
are kept out of the autosomal main table. The sample report summary is expanded
from the V2 sample summary, so zero-autosomal-report-event samples remain
visible.

The CNV plot contract is fixed to `calibrated_z` only. Missing or non-finite
bin values are skipped, and missing `calibrated_z` is an error. Each sample has
`{sample}.final_cnv.svg` and `{sample}.plot_bins.tsv`; plot states are limited
to `dup`, `del`, and `neutral`. The current SVG style uses yellow `dup` bins,
blue `del` bins, grey neutral bins, and red horizontal `report-z-trend` lines
only across final autosomal report event intervals. It does not draw genome-wide
or chromosome-wide smooth polylines.

Each 0615 sample also has a copy-number proxy plot:
`{sample}.final_cnv_cn.svg` plus `{sample}.plot_bins_cn.tsv`. The CN scatter is
drawn for every non-centromere 1Mb bin and is colored only by the bin-level
visual CN proxy, not by report-event intervals. The proxy is computed from each
bin's `calibrated_z` as `CN_proxy = clip(2 + calibrated_z * 0.05, 0, 4)`.
This maps the 30% mosaic review band to `CN=1.7..2.3`: bins below 1.7 are
colored deletion red, bins above 2.3 are colored duplication blue, and bins in
the interval are grey. Final report events still draw separate horizontal
event-level CN trend lines. The CN plot is a report visualization supplement,
not a bin-level CN caller and not a filtering input.

Materialized report-layer acceptance:

- Y1-Y8: truth 10/10, FN=0, report events 40, plots 8/8.
- H1-H16: truth 10/10, FN=0, H6 chr21 `report_weak_event`, report events
  23, plots 16/16.
- G1-G8: truth 10/10, FN=0, G2 truth not filtered, report events 26, plots
  8/8.
- 2026-06-15: no locked truth, burden/context only, report events 71, z plots
  5/5, CN plots 5/5.

This remains `development_only_not_final_release`.

### Reference

`h_r0_shadow_ref_20260619` is the current fixed shadow baseline for R&D. It
contains 38 selected reference samples, including H9, H10, H11, H12, H15, and
H16 as H-batch additions. It is not promoted as production-final.

Reference rebuild eligibility labels remain `R0/R1/R2`. They must not be reused
as `N0/N1/N2` background labels.

### Branch A

Branch A remains WisecondorX predict/CBS-derived candidate discovery. P2 no-FN
validation passed on Y1-Y8 and H1-H6 under the active reference.

Branch A burden optimization Phase 1 is now completed as config plumbing,
ablation evidence, and isolated `merge_gap_bp=2_000_000` materialization. The
workflow exposes `core.wisecondorx.cnv.branch_a` settings for `merge_gap_bp`,
`strong_z`, `output_dir`, `validation_dir`, and `log_dir`; defaults remain
`merge_gap_bp=0` and `strong_z=10.0`, preserving current P2 behavior.

`strong_z=10` is a strong/sensitive tier marker only. It is not a hard inclusion
or exclusion filter. Signals below 10 can still represent weak positives or
mosaic-positive candidates.

Remote materialization with default `merge_gap_bp=0` still gives FN=0 on
Y1-Y8/H1-H6 and keeps H6 chr21 detected. The explicit 2 Mb overlay has also
been materialized without overwriting default P2 outputs:

- Y1-Y8: candidates 131 -> 97, truth detected 10/10, FN=0.
- H1-H16: candidates 221 -> 105, truth detected 10/10, FN=0, H6 chr21
  retained.
- 2026-06-15: candidates 201 -> 165, no truth labels.

`merge_gap_bp=2_000_000` is not promoted as a default until the fixed A/B/S
chain is rerun and benchmarked under that explicit config.

### Branch B

Branch B is still needed for candidate-level refinement, including mosaic, LOH,
UPD, CN-like amplitude, region-risk, sample-noise, background, and
CNVpro-inspired consistency evidence. It must refine Branch A candidates only
and must not create B-only final events.

Current legacy/current-code Branch B output is excluded from current Branch B V2
performance evidence. Do not use `final_disposition`, `branch_b_keep_event`,
legacy/current-code Branch B artifact labels, old N0=0, or N1-only
matched-negative promotion as V2 evidence, benchmark truth, or report-release
criteria.

Branch B V2-only benchmark has now been materialized on the explicit Branch A
`merge_gap_bp=2_000_000` overlay. The benchmark consumes V2 classifier rows
only and explicitly ignores legacy/current-code Branch B decision fields.

Current V2 benchmark preservation result:

- Y1-Y8: 97 Branch A gap2m candidates; truth preserved 10/10; FN=0;
  hard-suppressed truth=0.
- H1-H16: 105 Branch A gap2m candidates; truth preserved 10/10; FN=0;
  hard-suppressed truth=0; H6 chr21 remains preserved.
- 2026-06-15: 165 Branch A gap2m candidates; no locked truth, context only.

This proves truth-overlap preservation and no hard suppression under the current
V2 shadow contract. It does not prove FP/review burden is solved and does not
promote Branch B V2 to final-report logic.

The first Branch B V2 burden refinement after the benchmark is also
materialized: `chrX`/`chrY` candidates now route to
`V2_SEX_CHROMOSOME_REVIEW` / `V2_ROUTE_BRANCH_S_REVIEW`, preserving evidence
tier and `none_shadow_only` final-report impact. This separates sex-chromosome
review burden from autosomal Branch B positive-support review without changing
sex calling or final SCA status.

Current class-level burden after sex routing:

- Y1-Y8: 97 candidates; 78 autosomal `V2_POSITIVE_SUPPORT_REVIEW`, 13
  `V2_SEX_CHROMOSOME_REVIEW`, 6 `V2_NO_CALL_CONTRACT_RISK`; truth preserved
  10/10, FN=0.
- H1-H16: 105 candidates; 57 autosomal `V2_POSITIVE_SUPPORT_REVIEW`, 42
  `V2_SEX_CHROMOSOME_REVIEW`, 6 `V2_NO_CALL_CONTRACT_RISK`; truth preserved
  10/10, FN=0.
- H7-H16 context subset: 57 candidate rows; 31 sex-chromosome review rows, 24
  autosomal positive-support review rows, and 2 no-call contract-risk rows.
- 2026-06-15: 165 candidate rows; 14 sex-chromosome review rows, 127 autosomal
  positive-support review rows, and 24 no-call contract-risk rows.

Y3/H5/H6 chrX truth events are now routed to Branch S review and remain
preserved. H6 chr21 remains autosomal positive-support review with
`top_a_abs_zscore=7.1135`.

The remaining autosomal Branch B V2 burden has been audited in
`docs/reports/branch_b_v2_autosomal_burden_audit_2026-06-21.md`.
The key result is that a tempting direction-support split is not safe as a hard
filter or universal downgrade: multiple locked autosomal truth top candidates
would be labeled `A_ONLY_NO_B_DIRECTION_SUPPORT` by that audit, including
Y2 chr14 gain, Y4 chr13 gain, H2 chr6 gain, H3 chr13 gain, and H4 chr15 gain.
H6 chr21 remains direction-supported. Therefore Branch B V2 can expose
direction-support as review evidence, but it must not hard-suppress or
final-demote candidates solely for weak Branch B-side direction support.

The review-label-only direction-support contract is now implemented and
materialized in `docs/reports/branch_b_v2_direction_support_label_2026-06-21.md`.
The classifier emits `v2_direction_support_label` and
`v2_direction_support_reason` without changing `v2_candidate_class`,
`v2_classifier_action`, or `v2_final_report_impact`.

Materialized direction-label counts:

- Y1-Y8: 97 rows; `B_DIRECTION_SUPPORTED=66`,
  `A_ONLY_WEAK_B_DIRECTION=20`, `B_DIRECTION_CONFLICT=11`.
- H1-H16: 105 rows; `B_DIRECTION_SUPPORTED=68`,
  `A_ONLY_WEAK_B_DIRECTION=26`, `B_DIRECTION_CONFLICT=11`.
- 2026-06-15: 165 rows; `B_DIRECTION_SUPPORTED=97`,
  `A_ONLY_WEAK_B_DIRECTION=40`, `B_DIRECTION_CONFLICT=28`.

Truth preservation remains unchanged after materialization:

- Y1-Y8: truth preserved 10/10, FN=0, hard-suppressed truth=0.
- H1-H16: truth preserved 10/10, FN=0, hard-suppressed truth=0.
- H6 chr21 remains `V2_POSITIVE_SUPPORT_REVIEW`,
  `B_DIRECTION_SUPPORTED`, `none_shadow_only`.

The unresolved background limitation is now explicit in workflow outputs as
review context. The classifier emits `v2_background_context_label` and
`v2_background_context_reason`, with `background_context_label_counts` in the
summary JSON. This does not change class/action/final impact.

Materialized background-context counts:

- Y1-Y8: 97 rows; `UNKNOWN_BACKGROUND_NO_NULL_SUPPORT=86`,
  `SHADOW_BACKGROUND_NO_NULL_SUPPORT=11`.
- H1-H16: 105 rows; `UNKNOWN_BACKGROUND_NO_NULL_SUPPORT=69`,
  `SHADOW_BACKGROUND_NO_NULL_SUPPORT=36`.
- 2026-06-15: 165 rows; `UNKNOWN_BACKGROUND_NO_NULL_SUPPORT=155`,
  `SHADOW_BACKGROUND_NO_NULL_SUPPORT=10`.

H6 chr21 remains `V2_POSITIVE_SUPPORT_REVIEW`,
`UNKNOWN_BACKGROUND_NO_NULL_SUPPORT`, `B_DIRECTION_SUPPORTED`,
`none_shadow_only`.

The Branch B V2 method reset loop has now materialized a safer interpretation
of those labels. Branch B V2 is Branch-A-anchored evidence stratification, not
legacy/current-code Branch B filtering. The classifier now emits:

- `v2_signal_strength_tier`;
- `v2_length_tier`;
- `v2_clean_support_label`;
- `v2_gc_rc_context_label`;
- `v2_b_signal_context_label`;
- `v2_b_signal_context_reason`;
- `v2_disposition`.

The old direction-conflict wording is replaced by B-side signal-context
wording. `B_SIGNAL_DISCORDANT_WITH_A_DIRECTION` means auxiliary B-side signal is
discordant with Branch A direction; it does not mean Branch B called an
opposite CNV event. GC/RC context remains auxiliary and cannot hard-suppress a
Branch A positive.

Materialized method-reset preservation result:

- Y1-Y8: candidates 97; truth preserved 10/10; FN=0; hard-suppressed truth=0.
- H1-H16: candidates 105; truth preserved 10/10; FN=0; hard-suppressed
  truth=0; H6 chr21 retained.
- 2026-06-15: candidates 165; no locked truth; status `skipped_no_truth`.

Materialized disposition counts:

- Y1-Y8: `background_unknown_review=84`, `sca_branch_s_review=13`.
- H1-H16: `background_unknown_review=63`, `sca_branch_s_review=42`.
- 2026-06-15: `background_unknown_review=151`, `sca_branch_s_review=14`.

This improves interpretability and truth-safety. It does not prove FP/review
burden is solved and does not promote Branch B V2 to final-report logic.

The next Branch B V2 truth-safe filter contract is now materialized on the
active gap2m benchmark path. The classifier emits explicit filter-action
fields:

- `v2_filter_version`;
- `v2_filter_action`;
- `v2_filter_reason`;
- `v2_filter_scope`;
- `v2_filter_hard_suppression_allowed`.

The benchmark now records filter actions in truth-overlap output and counts
filter-level hard suppression. Current hard suppression is allowed only for
workflow/reference contract risk, currently
`same_chrom_ref_leakage_contract_risk`. GC/RC attenuation, B-side signal
discordance, unknown background, length tier, low clean support, and high
region-risk burden cannot hard-suppress a Branch A candidate. Low clean support
can only downgrade to technical review.

Remote validation:

```text
tests/unit/test_branch_b_v2_classifier.py
tests/unit/test_branch_b_v2_benchmark.py
26 passed in 0.90s
```

Broader remote unit tests:

```text
55 passed in 1.37s
```

Snakemake dry-runs passed for
`branch_b_v2_benchmark branch_s_review cnv_report` under all three active
gap2m configs. Forced materialization refreshed the V2 classifier and
benchmark outputs.

Materialized result after the filter-action contract:

- Y1-Y8: candidates 97; truth preserved 10/10; FN=0; hard-suppressed truth=0;
  filter-suppressed candidates=0.
- H1-H16: candidates 105; truth preserved 10/10; FN=0; hard-suppressed
  truth=0; filter-suppressed candidates=0.
- 2026-06-15: candidates 165; no locked truth; status `skipped_no_truth`;
  filter-suppressed candidates=0.

Remote output files include `v2_filter_action`, `top_v2_filter_action`, and
`v2_filter_suppressed_count`.

The next Branch B V2 burden stratification loop is now materialized in
`docs/reports/branch_b_v2_burden_stratification_2026-06-21.md`. It adds
`v2_burden_reduction_tier`, `v2_burden_reduction_action`,
`v2_burden_reduction_reason`, and `v2_burden_evidence_tags`.

CNVseq/CNVpro-derived concepts are explicitly marked as evidence tags only:

- `[CNVpro-inspired]` length tiers for report/review routing;
- `[CNVpro-confirmed]` acrocentric qter context for chr13/14/15/21/22 review;
- `[CNVseq-asset]` mask/mappability/PAR/sex-homology annotation context;
- `[CNVpro-like]` GC/RC context;
- `[Not used]` CNVcalling.R/cghFLasso as primary caller replacements.

Materialized burden-stratification result:

- Y1-Y8: candidates 97; truth preserved 10/10; FN=0; hard-suppressed truth=0;
  84 `background_unknown_review`; 13 `branch_s_review`.
- H1-H16: candidates 105; truth preserved 10/10; FN=0; hard-suppressed
  truth=0; H6 chr21 retained; 63 `background_unknown_review`; 42
  `branch_s_review`.
- 2026-06-15: candidates 165; no locked truth; context only; 151
  `background_unknown_review`; 14 `branch_s_review`.

This improves auditability and future report/review routing, but it still does
not prove FP/review burden reduction. It must not be treated as Branch B V2
final-report promotion.

The Branch B V2 report-contract integration loop is now represented in
`docs/reports/branch_b_v2_report_contract_2026-06-21.md`. `cnv_report` can
consume the V2 benchmark summary and sample-summary table and expose:

- `branch_b_v2_burden_status`;
- `branch_b_v2_background_unknown_review_count`;
- `branch_b_v2_branch_s_review_count`;
- `branch_b_v2_technical_risk_review_count`;
- `branch_b_v2_report_candidate_count`;
- `branch_b_v2_legacy_fields_used=false`;
- `branch_b_v2_final_impact=development_review_only`.

These fields are report-display context only. They do not modify Branch A,
do not add hard filters, and do not promote Branch B V2 or Branch S to final
report logic.

The earlier G1-G8 positive cohort has also been evaluated under the same active
reference and explicit gap2m Branch A overlay. Existing G-batch BAM files were
reused; WisecondorX convert/predict, Branch A, Branch B V2, Branch S,
benchmark, and report were rerun in the remote mirror. Current G1-G8 result:

- candidates=75;
- locked truth events=10;
- truth preserved=10/10;
- FN=0;
- truth report-layer filtered=0;
- report events=15;
- internal review events=40;
- filtered audit-only events=7;
- Branch S events=13.

This supports the current sensitivity/preservation gate but does not promote
the shadow reference, Branch B V2, or Branch S to production-final status. G2
remains the main G-batch review-burden outlier and should be reviewed before
designing a second demotion pass.

The second Branch B V2 report-burden pass has been retracted as a current
candidate-level decision rule. It is documented only as superseded methodology
audit in `docs/reports/branch_b_v2_report_event_audit_2026-06-21.md`.

Current correction report:
`docs/reports/branch_b_v2_pass2_correction_2026-06-21.md`.

Corrected rule boundary:

```text
sample report_event count >= 3
```

is a sample-level burden audit flag only. It must not change a candidate's
`v2_report_layer_class`, `v2_report_visibility`, `v2_filter_action`, or
`v2_burden_reduction_action`.

Benchmark output still includes `report_events.tsv` and `report_events.json`.
The corrected classifier/benchmark/report outputs have been rematerialized for
Y, H, G, and 2026-06-15. Corrected materialized counts:

- Y1-Y8: candidates=97, truth preserved 10/10, FN=0, truth filtered=0,
  report=21, internal_review=50, filtered=13, Branch S=13,
  sample burden flags=2.
- H1-H16: candidates=105, truth preserved 10/10, FN=0, truth filtered=0,
  report=6, internal_review=49, filtered=8, Branch S=42,
  sample burden flags=1; H6 chr21 remains visible.
- G1-G8: candidates=75, truth preserved 10/10, FN=0, truth filtered=0,
  report=15, internal_review=40, filtered=7, Branch S=13,
  sample burden flags=1; G2 truth remains visible.
- 2026-06-15: candidates=165, no locked truth, context only,
  report=52, internal_review=76, filtered=23, Branch S=14,
  sample burden flags=5.

Further rules must inspect candidate-level evidence, not sample report count.

The 2Mb/3Mb low-resolution and ref-MAD evidence interface is now implemented
as an optional Branch B V2 auxiliary context layer. It is disabled by default.
When enabled, it computes:

- per-candidate 2Mb/3Mb same-direction support labels;
- per-bin reference median/MAD/CV/MAD-z stability context from configured
  reference NPZs;
- per-event ref-MAD summaries;
- V2 classifier context fields `v2_lowres_context_label`,
  `v2_lowres_context_reason`, and `v2_ref_stability_context`.

This interface does not change Branch A, does not replace WisecondorX/CBS, does
not create B-only events, and does not let low-resolution absence change
report/filter class by itself. It has not yet built or materialized
`h_r0_shadow_ref_20260619_2mb` or `h_r0_shadow_ref_20260619_3mb`; those remain
the next long-task gate.

### Branch S

Branch S is not final SCA, but it must not be omitted from report development.
Current status is `review_reportable_with_limitations`: report packages should
carry a visible SCA/sex-chromosome section with output mode and uncertainty,
while SCA confirmation/suppression remains blocked until locked SCA truth
coverage is available.

The immediate Branch S development target is to control SCA false positives in
reference/negative-like samples and preserve sex-call concordance. Final SCA
promotion still needs dedicated truth coverage.

### P6 / Report

The 2026-06-15 P6 package already materialized is historical
development-only evidence. That historical package is not final-release
evidence.

Future P6/report remains the delivery target after Branch A burden optimization
and Branch B V2 evidence/disposition strengthening. The report must be workflow
generated, must carry the same reference/config contract, and must show Branch
A, Branch B, Branch S, background source, and limitations explicitly.

## Stale Routes To Exclude

- legacy/current-code Branch B as Branch B V2 performance evidence
- `final_disposition` as V2 decision evidence
- `branch_b_keep_event` as V2 decision evidence
- old N0=0 as a post-rebuild blocker
- N1-only matched-negative promotion
- treating P3 evidence ledger as V2 performance validation
- treating current V2 shadow classifier output as V2-only performance evidence
- treating Branch S as final SCA
- treating the old 0615 development report package as final release evidence

## Next Gates

1. Keep this context index current and commit it with each completed R&D loop.
2. Keep workflow source line endings LF on the remote mirror before making new
   Snakemake validation claims. The 2026-06-21 parser blocker was repaired by
   normalizing `Snakefile`, `.smk`, Python, and YAML workflow sources to LF.
3. Use the materialized explicit `merge_gap_bp=2_000_000` overlay for the next
   Branch B V2-only benchmark unless a blocking regression appears; keep
   default `merge_gap_bp=0` as rollback/control.
4. Continue refining Branch B V2 evidence/disposition for remaining autosomal
   report/internal-review burden while preserving the materialized
   no-FN/no-truth-filtered benchmark.
5. Direction support is now materialized as a review label only; do not turn it
   into a hard filter, final benign/artifact call, or universal demotion.
6. Background context is now materialized as a review label only; do not treat
   `UNKNOWN_BACKGROUND_NO_NULL_SUPPORT` as benign or background-compatible.
7. B-side signal context and length tiers are review/disposition context only;
   do not convert them into hard filters without a locked-truth ablation.
8. The V2 filter-action fields are now materialized and truth-preserving under
   locked Y/H truth. Do not treat this as FP/review-burden reduction yet. Only
   `suppress_workflow_contract_risk` can count as a hard suppression action in
   the current contract.
9. The V2 burden-stratification fields are now materialized and truth-preserving
   under locked Y/H truth, with explicit CNVpro/CNVseq evidence tags. Treat this
   as audit/report routing context, not as proof that FP/review burden has been
   reduced. Do not promote length, clean support, GC/RC, acrocentric, sex
   homology, or B-side signal tags into hard filters without locked-truth
   ablation and explicit H6 chr21 preservation.
10. Upgrade Branch S toward review-reportable output with controlled negative/ref
   FP burden; final SCA promotion remains a separate truth gate.
11. Generate the next P6/report package only after the fixed A/B/S contracts are
   represented in workflow outputs.
12. The next Branch B V2 pass must start from remaining `report_events.tsv`
   rows, especially G2 and 0615 sample `JZ26125846-61-61`; do not add broad
   demotion/filter rules without Y/H/G locked-truth ablation.
