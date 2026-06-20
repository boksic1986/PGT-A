# H7-H16 Reference Cohort Decision

Date: 2026-06-20

## Scope

This is the G1/P1 decision artifact derived from the H7-H16 reference-candidate
audit. It records which H samples can be used in the next named shadow
reference contract before Branch A no-FN validation.

This is not a production-reference promotion and not a formal N0 background
definition.

## Controlling Documents

- `docs/reports/branch_ab_v2_gate_phase_plan_2026-06-20.md`
- `docs/reports/branch_ab_v2_rnd_constraints_2026-06-18.md`
- `docs/reports/h7_h16_reference_candidate_audit_2026-06-20.md`

## Primary Decision

Use the existing R0-only H shadow reference as the primary P2 Branch A
validation contract:

```text
reference_id: h_r0_shadow_ref_20260619
reference_config: config_reference_h_r0_shadow_20260619.yaml
predict_audit_config: config_predict_h_reference_audit_20260620.yaml
reference_root: /data/project/CNV/PGT-A/refactor_validation_20260419/results_h_r0_shadow_ref_20260619/reference
selected_samples: /data/project/CNV/PGT-A/refactor_validation_20260419/results_h_r0_shadow_ref_20260619/reference/cohorts/pass_warn_h_r0/selected_samples.txt
```

Current selected reference sample count is 38. The H additions are:

```text
H9,H10,H11,H12,H15,H16
```

## R Labels

| label | samples | decision |
|---|---|---|
| R0 | H9, H10, H11, H12, H15, H16 | keep in primary H R0 shadow reference for P2 Branch A validation |
| R1 | H7, H8, H13, H14 | exclude from primary; use only in R0+R1 ablation or shadow review |
| R2 | none | no H7-H16 sample is hard-excluded by this audit |

The machine-readable decision table is:

```text
docs/reports/h7_h16_reference_cohort_decision_2026-06-20.tsv
```

## Interpretation

The R0 samples are reference-rebuild candidates only. They are not clean N0
controls, and they must not be used as formal matched-negative empirical null
samples without a cross-fit or leave-one-out design.

The R1 samples remain useful for ablation because they are QC PASS/XY samples
from the same H batch, but their altered-genome fraction or sample-specific
Branch A burden is too high for silent inclusion in the primary shadow
reference.

The decision intentionally avoids rejecting H7-H16 solely because the older
32-reference analysis produced Branch A signals. The labels use the
post-rebuild/ref-specific audit output and current QC/gender/reference context.

## P2 Entry Contract

P2 should now proceed with:

```text
WisecondorX predict
-> Branch A candidates
-> Branch A no-FN validation
```

using the same `h_r0_shadow_ref_20260619` reference/config contract.

Required P2 sample groups:

- Y1-Y8 known positives;
- H1-H6 known positives, including explicit H6 chr21 tracking;
- H7-H16 reference-candidate burden review;
- 2026-06-15 exploratory samples for burden review only.

Branch B, Branch S, evaluation, benchmark, and report outputs generated before
P2/G2 passes remain diagnostic snapshots only. They must not be used as the
report line or promotion evidence.

## Verification State

Remote audit outputs exist at:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419/results_h_20260608_h_r0_shadow_ref_20260619/wisecondorx/cnv/reference_audit/reference_candidate_audit.tsv
/data/project/CNV/PGT-A/refactor_validation_20260419/results_h_20260608_h_r0_shadow_ref_20260619/wisecondorx/cnv/reference_audit/reference_candidate_audit.summary.json
```

Remote selected sample file was checked and contains 38 samples, with H samples
limited to `H9,H10,H11,H12,H15,H16`.

## Current Limitations

- `h_r0_shadow_ref_20260619` remains a shadow reference, not final production.
- R0/R1/R2 labels do not define N0/N1/N2.
- Formal N0/cross-fit background is still a later promotion gate.
- SCA/Branch S remains review/development-only until its own validation gate.
- 2026-06-15 reports must not be released before fixed Branch A, Branch B, and
  Branch S contracts pass their gates.
