---
status: active_current_index
decision_use: active_context_entrypoint
last_updated: 2026-06-21
active_reference_id: h_r0_shadow_ref_20260619
---

# Current Context Index

This file is the first context-restoration entrypoint for the current PGT-A
CNV-seq R&D cycle. It does not replace result files, configs, or workflow code.
It pins the current decision hierarchy so future work does not resume from stale
handoffs or legacy Branch B outputs.

## Required Read Order

1. `docs/CURRENT_CONTEXT_INDEX.md`
2. `docs/handoff/2026-06-21_1415_branch_b_v2_method_reset_handoff.md`
3. `AGENTS.md`
4. `skills/conversation_handoff/SKILL.md`
5. `skills/pgta_reference_modeling_analysis/SKILL.md`
6. Current result files and configs for the active task

## Active Inputs

active_handoff: docs/handoff/2026-06-21_1415_branch_b_v2_method_reset_handoff.md
active_reference_id: h_r0_shadow_ref_20260619
reference_status: fixed_shadow_baseline_not_production
remote_snakemake_parse_status: repaired_lf_normalized_2026-06-21
branch_a_status: burden_phase1_gap2m_materialized_default_unchanged
branch_b_status: v2_gap2m_method_reset_disposition_materialized_truth_preserved
branch_s_status: review_reportable_with_limitations
report_status: final_delivery_target_after_a_b_strengthening

## Current Evidence Files

- `docs/reports/p1_p6_result_credibility_audit_2026-06-21.md`
- `docs/reports/h7_h16_reference_cohort_decision_2026-06-20.md`
- `docs/reports/branch_a_validation_h_r0_shadow_2026-06-20.md`
- `docs/reports/branch_a_burden_optimization_phase1_2026-06-21.md`
- `docs/reports/branch_b_v2_gap2m_benchmark_2026-06-21.md`
- `docs/reports/branch_b_v2_sex_route_refinement_2026-06-21.md`
- `docs/reports/branch_b_v2_autosomal_burden_audit_2026-06-21.md`
- `docs/reports/branch_b_v2_direction_support_label_2026-06-21.md`
- `docs/reports/branch_b_v2_background_context_label_2026-06-21.md`
- `docs/reports/branch_b_v2_method_reset_threshold_inventory_2026-06-21.md`
- `docs/reports/branch_b_v2_reference_background_and_sca_design_2026-06-20.md`
- `docs/reports/branch_s_p5_report_boundary_2026-06-20.md`
- `docs/reports/p6_report_package_contract_2026-06-20.md`
- `docs/reports/branch_ab_v2_gate_phase_plan_2026-06-20.md`
- `docs/reports/branch_ab_v2_rnd_constraints_2026-06-18.md`

## Current Module State

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
   FP and review burden while preserving the materialized
   no-FN/no-hard-suppression benchmark.
5. Direction support is now materialized as a review label only; do not turn it
   into a hard filter, final benign/artifact call, or universal demotion.
6. Background context is now materialized as a review label only; do not treat
   `UNKNOWN_BACKGROUND_NO_NULL_SUPPORT` as benign or background-compatible.
7. B-side signal context and length tiers are review/disposition context only;
   do not convert them into hard filters without a locked-truth ablation.
8. Upgrade Branch S toward review-reportable output with controlled negative/ref
   FP burden; final SCA promotion remains a separate truth gate.
9. Generate the next P6/report package only after the fixed A/B/S contracts are
   represented in workflow outputs.
