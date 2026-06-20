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
2. `docs/handoff/2026-06-21_0127_p1_p6_credibility_audit_handoff.md`
3. `AGENTS.md`
4. `skills/conversation_handoff/SKILL.md`
5. `skills/pgta_reference_modeling_analysis/SKILL.md`
6. Current result files and configs for the active task

## Active Inputs

active_handoff: docs/handoff/2026-06-21_0127_p1_p6_credibility_audit_handoff.md
active_reference_id: h_r0_shadow_ref_20260619
reference_status: fixed_shadow_baseline_not_production
remote_snakemake_parse_status: repaired_lf_normalized_2026-06-21
branch_a_status: p2_no_fn_passed_burden_optimization_next
branch_b_status: v2_evidence_design_legacy_excluded
branch_s_status: review_reportable_with_limitations
report_status: final_delivery_target_after_a_b_strengthening

## Current Evidence Files

- `docs/reports/p1_p6_result_credibility_audit_2026-06-21.md`
- `docs/reports/h7_h16_reference_cohort_decision_2026-06-20.md`
- `docs/reports/branch_a_validation_h_r0_shadow_2026-06-20.md`
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
validation passed on Y1-Y8 and H1-H6 under the active reference. The next Branch
A task is burden optimization while preserving known-positive recall, especially
the weak H6 chr21 event.

`strong_z=10` is a strong/sensitive tier marker only. It is not a hard inclusion
or exclusion filter. Signals below 10 can still represent weak positives or
mosaic-positive candidates.

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

The next Branch B gate is to implement V2-only evidence/disposition evaluation
against Y1-Y8 and H1-H6 truth, using reference-cohort background only as a
limited and explicitly labeled context source until formal N0/cross-fit evidence
exists.

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
3. Evaluate Branch A burden under the fixed reference while preserving FN=0 on
   Y1-Y8/H1-H6 and preserving H6 chr21.
4. Implement Branch B V2-only truth benchmark and evidence/disposition output,
   excluding legacy/current-code Branch B fields from decisions.
5. Upgrade Branch S toward review-reportable output with controlled negative/ref
   FP burden; final SCA promotion remains a separate truth gate.
6. Generate the next P6/report package only after the fixed A/B/S contracts are
   represented in workflow outputs.
