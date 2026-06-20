# H7-H16 Reference Candidate Audit

Date: 2026-06-20

## Scope

This report records the first G1/P1 reference-candidate audit run for H7-H16
under the current H-augmented shadow reference.

This is not a production-reference promotion decision. It is an audit artifact
for the next reference reset step.

## Context

Controlling plan:

```text
docs/reports/branch_ab_v2_gate_phase_plan_2026-06-20.md
```

Current constraints:

```text
docs/reports/branch_ab_v2_rnd_constraints_2026-06-18.md
```

Remote code mirror:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419
```

Reference ID:

```text
h_r0_shadow_ref_20260619
```

Workflow config:

```text
config_predict_h_reference_audit_20260620.yaml
```

Remote output:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419/results_h_20260608_h_r0_shadow_ref_20260619/wisecondorx/cnv/reference_audit/reference_candidate_audit.tsv
/data/project/CNV/PGT-A/refactor_validation_20260419/results_h_20260608_h_r0_shadow_ref_20260619/wisecondorx/cnv/reference_audit/reference_candidate_audit.summary.json
```

## Validation Commands

Unit tests:

```bash
ssh fengxian "cd /data/project/CNV/PGT-A/refactor_validation_20260419 && PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 /biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest tests/unit/test_reference_candidate_audit.py tests/unit/test_branch_ab_phase12_workflow_contract.py -q"
```

Result:

```text
12 passed in 0.63s
```

Snakemake dry-run:

```bash
ssh fengxian "cd /data/project/CNV/PGT-A/refactor_validation_20260419 && /biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile /data/project/CNV/PGT-A/refactor_validation_20260419/config_predict_h_reference_audit_20260620.yaml --cores 1 -n reference_audit"
```

Result:

```text
2 jobs: cnv_reference_candidate_audit, reference_audit
No fastp_bwa / mapping / convert / predict rerun required.
```

Materialization:

```bash
ssh fengxian "cd /data/project/CNV/PGT-A/refactor_validation_20260419 && /biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile /data/project/CNV/PGT-A/refactor_validation_20260419/config_predict_h_reference_audit_20260620.yaml --cores 1 reference_audit"
```

Result:

```text
2 of 2 steps done.
```

## Audit Summary

The audit uses post-rebuild/ref-specific Branch A candidate tables, CNV QC,
gender calls, current H R0 shadow reference membership, bin annotations, and the
existing Branch B evidence ledger only as optional region-risk context.

It does not use old Branch B kept counts for R-label assignment.

```json
{
  "formal_n0_used": false,
  "legacy_branch_b_kept_counts_used_for_r_label": false,
  "label_counts": {
    "R0": 6,
    "R1": 4
  },
  "sample_count": 10
}
```

## Sample-Level Result

| sample | QC | sex | membership | A candidates | strong | sensitive | broad | altered genome fraction | risk burden | shared batch signal | sample-specific signal | R label | reason | action |
|---|---|---|---|---:|---:|---:|---:|---:|---:|---|---|---|---|---|
| H7 | PASS | XY | candidate_only | 14 | 12 | 2 | 5 | 0.073167 | 0.167748 | yes | no | R1 | altered_genome_fraction>=0.05 | shadow_or_ablation_only |
| H8 | PASS | XY | candidate_only | 15 | 12 | 3 | 7 | 0.080919 | 0.171014 | yes | no | R1 | altered_genome_fraction>=0.05 | shadow_or_ablation_only |
| H9 | PASS | XY | included_in_evaluated_reference | 12 | 12 | 0 | 4 | 0.042398 | 0.171490 | yes | no | R0 | branch_a_signal_shared_across_batch_not_sample_specific | candidate_reference_with_batch_context |
| H10 | PASS | XY | included_in_evaluated_reference | 12 | 12 | 0 | 4 | 0.042398 | 0.170551 | yes | no | R0 | branch_a_signal_shared_across_batch_not_sample_specific | candidate_reference_with_batch_context |
| H11 | PASS | XY | included_in_evaluated_reference | 12 | 12 | 0 | 4 | 0.042398 | 0.176364 | yes | no | R0 | branch_a_signal_shared_across_batch_not_sample_specific | candidate_reference_with_batch_context |
| H12 | PASS | XY | included_in_evaluated_reference | 13 | 13 | 0 | 5 | 0.041671 | 0.174772 | yes | no | R0 | branch_a_signal_shared_across_batch_not_sample_specific | candidate_reference_with_batch_context |
| H13 | PASS | XY | candidate_only | 23 | 14 | 9 | 13 | 0.198664 | 0.251563 | yes | yes | R1 | sample_specific_branch_a_signal;altered_genome_fraction>=0.05;high_repeat_or_acrocentric_burden>=0.2 | shadow_or_ablation_only |
| H14 | PASS | XY | candidate_only | 22 | 12 | 10 | 14 | 0.217804 | 0.175860 | yes | no | R1 | altered_genome_fraction>=0.05 | shadow_or_ablation_only |
| H15 | PASS | XY | included_in_evaluated_reference | 12 | 12 | 0 | 4 | 0.042398 | 0.171404 | yes | no | R0 | branch_a_signal_shared_across_batch_not_sample_specific | candidate_reference_with_batch_context |
| H16 | PASS | XY | included_in_evaluated_reference | 12 | 12 | 0 | 4 | 0.042398 | 0.166904 | yes | no | R0 | branch_a_signal_shared_across_batch_not_sample_specific | candidate_reference_with_batch_context |

## Interpretation

Current audit labels:

- `R0`: H9, H10, H11, H12, H15, H16.
- `R1`: H7, H8, H13, H14.
- `R2`: none in this audit.

The R0 samples are not clean N0 controls. They are reference-rebuild candidates
with shared batch/context signals. They can support a named shadow reference
variant, but they do not prove a formal matched-negative background.

The R1 samples should remain shadow/ablation-only until reference reset and
Branch A no-FN validation are complete.

## Current G1 Status

The audit portion of G1 is complete. The cohort decision is recorded separately
in:

```text
docs/reports/h7_h16_reference_cohort_decision_2026-06-20.md
docs/reports/h7_h16_reference_cohort_decision_2026-06-20.tsv
```

This permits moving into P2 Branch A no-FN validation under the
`h_r0_shadow_ref_20260619` shadow reference contract. It does not promote the
reference to production.

Completed:

- formal audit module and workflow target exist;
- H7-H16 audit is materialized remotely;
- audit labels do not use old Branch B kept counts;
- Snakemake dry-run confirms no mapping/BAM regeneration is required for this
  audit target;
- R0/R1/R2 cohort decision is now written to a stable decision artifact.

Still required before any downstream promotion:

- rerun WisecondorX predict and Branch A under the chosen reference/config for
  P2;
- validate Y1-Y8 and H1-H6 no-FN behavior before moving to Branch B evidence
  refinement.
