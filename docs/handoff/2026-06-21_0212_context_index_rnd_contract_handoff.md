# 2026-06-21 02:12 Context Index And R&D Contract Handoff

## 1. Goal

Update the whole current PGT-A CNV-seq R&D context recovery mechanism after the
user clarified that:

- Branch S is not final SCA, but report development must still carry an SCA /
  sex-chromosome section with explicit uncertainty.
- P6/report should remain the delivery target after strengthening Branch A and
  Branch B, rather than being permanently classified as development-only.
- The correction applies to the whole R&D cycle, not only G1.

No algorithm thresholds, CNV calling logic, mosaic logic, sex calling logic, or
result schema were changed.

## 2. Context Read

- `docs/handoff/2026-06-21_0127_p1_p6_credibility_audit_handoff.md`
- `AGENTS.md`
- `skills/conversation_handoff/SKILL.md`
- `skills/pgta_reference_modeling_analysis/SKILL.md`
- `CURRENT_STATE.md`
- `PLANS.md`
- `docs/reports/p1_p6_result_credibility_audit_2026-06-21.md`
- `docs/reports/h7_h16_reference_cohort_decision_2026-06-20.md`
- `docs/reports/branch_a_validation_h_r0_shadow_2026-06-20.md`
- `docs/reports/branch_ab_v2_gate_phase_plan_2026-06-20.md`
- `docs/reports/branch_ab_v2_rnd_constraints_2026-06-18.md`
- `docs/reports/branch_b_v2_reference_background_and_sca_design_2026-06-20.md`
- `docs/reports/branch_s_p5_report_boundary_2026-06-20.md`
- `docs/reports/p6_report_package_contract_2026-06-20.md`

Active handoff used:

```text
docs/handoff/2026-06-21_0127_p1_p6_credibility_audit_handoff.md
```

## 3. Completed

Added:

```text
docs/CURRENT_CONTEXT_INDEX.md
tests/unit/test_current_context_index.py
docs/handoff/2026-06-21_0212_context_index_rnd_contract_handoff.md
```

Updated:

```text
CURRENT_STATE.md
PLANS.md
docs/reports/p1_p6_result_credibility_audit_2026-06-21.md
docs/reports/h7_h16_reference_cohort_decision_2026-06-20.md
docs/reports/branch_a_validation_h_r0_shadow_2026-06-20.md
docs/reports/branch_ab_v2_gate_phase_plan_2026-06-20.md
docs/reports/branch_ab_v2_rnd_constraints_2026-06-18.md
docs/reports/branch_b_v2_reference_background_and_sca_design_2026-06-20.md
docs/reports/branch_s_p5_report_boundary_2026-06-20.md
docs/reports/p6_report_package_contract_2026-06-20.md
```

## 4. Current Conclusions

- `docs/CURRENT_CONTEXT_INDEX.md` is now the first context-restoration
  entrypoint for this R&D cycle.
- `h_r0_shadow_ref_20260619` remains the fixed shadow baseline, not production
  final.
- P2 Branch A no-FN evidence remains current sensitivity evidence, but Branch A
  candidate burden optimization is still required under no-FN and H6 chr21
  preservation constraints.
- Legacy/current-code Branch B fields remain excluded from Branch B V2
  performance evidence and report-release decisions.
- Branch B remains necessary for candidate refinement, including mosaic, LOH,
  UPD, CN-like amplitude, region-risk, sample-noise, background, and
  CNVpro-inspired evidence.
- Branch S is not final SCA, but it should be carried into report development as
  `review_reportable_with_limitations`.
- The already materialized 0615 P6 package is historical
  `development_only_not_final_release` evidence. Future P6/report remains the
  delivery target after A/B/S contracts are strengthened.

## 5. Verification

Local static check:

```text
executable: git
command: git diff --check
result: no whitespace errors; only LF/CRLF warnings from the Windows worktree
```

Local dependency note:

```text
executable: python
command: python -m pytest tests/unit/test_current_context_index.py -q
result: not run; local Python lacks pytest and local tests are not workflow validation
```

Local direct static check:

```text
executable: python
command: direct invocation of tests/unit/test_current_context_index.py test functions
result: direct static context-index tests passed
```

Remote validation:

```text
executable: ssh fengxian
command: cd /data/project/CNV/PGT-A/refactor_validation_20260419 && PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 /biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest tests/unit/test_current_context_index.py -q
result: 3 passed in 0.15s
```

Remote sync note:

- Formal docs and the new test were copied into
  `/data/project/CNV/PGT-A/refactor_validation_20260419`.
- A first `scp` command briefly placed four root/test files under
  `docs/reports/`; those four mistaken remote copies were removed immediately:
  `CURRENT_CONTEXT_INDEX.md`, `CURRENT_STATE.md`, `PLANS.md`,
  `test_current_context_index.py`.

## 6. Not Changed

- No CNV calling logic.
- No mosaic logic.
- No sex calling logic.
- No result schema.
- No Snakemake rule behavior.
- No local result package, shell helper, Excel sample sheet, image, or
  `validate_npz.py` file was promoted into the formal change set.

## 7. Remaining Risks

- The previous remote Snakemake parse-time issue is not resolved by this
  documentation/test change. Do not make new Snakemake validation claims until
  the remote mirror is re-synced or repaired and a dry-run passes.
- Branch B V2 performance is still not validated. The next implementation gate
  must exclude legacy/current-code Branch B fields.
- Branch S final SCA remains blocked by insufficient locked SCA truth coverage.

## 8. Recommended Next Step

Proceed in this order:

1. Commit this context-index and contract cleanup.
2. Re-sync the remote mirror from committed main.
3. Fix/verify the remote Snakefile parse issue before new workflow dry-run
   claims.
4. Start Branch A burden optimization under fixed reference with no-FN and H6
   chr21 gates.
5. Implement Branch B V2-only benchmark/evidence disposition using Y1-Y8 and
   H1-H6 truth.
6. Carry Branch S into report development as review-reportable with explicit
   limitations.

## 9. Core File Sync

- `AGENTS.md`: not updated; no repository hard-rule change.
- `REPO_MAP.md`: not updated; no stable structure/entrypoint change beyond the
  new context index.
- `PLANS.md`: updated; active plan now covers the whole R&D cycle and report
  delivery target.
- `CURRENT_STATE.md`: updated; current state now reflects Branch S
  review-reportable limitations and P6/report as future delivery target.
