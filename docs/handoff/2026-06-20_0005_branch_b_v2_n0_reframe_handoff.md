# 2026-06-20 00:05 Branch B V2 N0 Reframe Handoff

## Goal

Correct the Branch B V2 R&D documentation before code changes. The key
correction is that legacy Branch B and the legacy negative-bank seed cannot be
used as the final `N0` benchmark for Branch B V2. After this documentation PR is
reviewed and merged, implementation should continue from the corrected V2
validation contract.

## Context Read

- `docs/handoff/2026-06-19_2145_h_r0_branch_b_refresh_handoff.md`
- `AGENTS.md`
- `skills/conversation_handoff/SKILL.md`
- `skills/pgta_reference_modeling_analysis/SKILL.md`
- `docs/reports/branch_ab_v2_rnd_constraints_2026-06-18.md`
- `docs/reports/branch_ab_v2_negative_bank_readiness_2026-06-19.md`
- `docs/reports/branch_ab_v2_negative_bank_seed_2026-06-18.tsv`
- `pgta/predict/branch_b/negative_bank.py`

## Documentation Changes

### Branch A/B V2 R&D constraints

File:
`docs/reports/branch_ab_v2_rnd_constraints_2026-06-18.md`

Updated contract:

- `branch_ab_v2_negative_bank_seed_2026-06-18.tsv` is a `legacy safety seed`,
  not a Branch B V2 benchmark truth table.
- Legacy/current Branch B kept or review status can explain why empirical null
  is disabled, but it cannot be the final authority for V2 `N0`.
- Branch B V2 validation order is:
  1. positive recall gate on Y1-Y8 and H1-H6;
  2. candidate-level V2 evidence gate;
  3. independent negative-background gate;
  4. V2 negative benchmark gate;
  5. final-report promotion gate.
- V2 positive validation must be at least no worse than the frozen comparator
  and should preserve weak truth-overlap candidates as review instead of hard
  suppression.
- Phase 2 now creates provisional N0/N1/N2 labels.
- Phase 3 may interpret matched-negative percentiles only with validated N0,
  held-out, or cross-fit background. Otherwise it must emit
  `UNKNOWN_BACKGROUND`.
- Phase 4 remains `REVIEW_REQUIRED` / `none_shadow_only` when no validated
  background exists.

### Negative-bank readiness audit

File:
`docs/reports/branch_ab_v2_negative_bank_readiness_2026-06-19.md`

Updated contract:

- The audit is a legacy safety-seed audit, not a V2 benchmark conclusion.
- H7-H16 `N0=0` only means the legacy seed has no locked-clean negative samples.
- After a post-rebuild reference run, N0/N1/N2 must be recomputed using the new
  reference ID, Branch A candidates, Branch B V2 evidence, and independent
  clean-negative review.
- Samples used to build the same reference are not ideal empirical-null controls
  unless a leave-one-out or cross-fit design is used.

## Current Effective State

- H R0 shadow ref was completed before this handoff.
- H9, H10, H11, H12, H15, and H16 were included in the shadow reference.
- Post-refresh Y1-Y8 and H1-H16 truth recall were both 1.0 in the recorded
  validation context.
- Current matched-negative status remains:
  - `N0=0`
  - `matched_negative_ready=false`
  - `matched_negative_blocking_reason=no_n0_locked_clean_negative_samples`
- This blocker applies only to the configured legacy seed and current
  matched-negative run. It is not the final Branch B V2 benchmark.

## Next Implementation Gates

After PR review and merge:

1. Validate Branch B V2 on Y1-Y8 and H1-H6 positives first.
2. Recompute H7-H16 stratification using post-H-R0 V2 evidence.
3. Design a held-out or cross-fit negative background.
4. Evaluate matched-negative percentile only after a valid background exists.
5. Produce the 2026-06-15 five-sample report only after the V2/report contract is
   clear for the selected workflow state.

## Do Not

- Do not use old Branch B kept=0 as the only standard for V2 `N0`.
- Do not treat old-reference `N0=0` as the final post-rebuild conclusion.
- Do not start Branch B V2 code changes before this documentation PR is
  reviewed.
- Do not claim validation from local execution; real workflow validation must run
  on `ssh fengxian`.
