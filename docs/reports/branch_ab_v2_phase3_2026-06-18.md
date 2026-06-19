# Branch A/B V2 Phase 3 Implementation Report

Date: 2026-06-18

## Scope

This report implements Phase 3 from `docs/reports/branch_ab_v2_rnd_constraints_2026-06-18.md`.

Phase 3 adds matched-negative empirical percentile as a shadow feature for Branch A/B V2 evidence. It does not change:

- WisecondorX predict
- Branch A candidate inclusion
- legacy Branch B final events
- `cnv_report`
- artifact hard filtering

## Phase 3 Contract

Input evidence:

- per-sample Phase 1 evidence ledger
- Phase 2 negative-bank labels
- optional external background evidence ledgers

Output paths:

```text
wisecondorx/cnv/postprocess/matched_negative/{sample}.candidate_evidence.tsv
wisecondorx/cnv/postprocess/matched_negative/{sample}.summary.json
```

Implementation:

```text
pgta/predict/branch_b/matched_negative.py
scripts/predict.py matched_negative_percentile
rule cnv_branch_b_matched_negative
```

New shadow columns include:

```text
matched_negative_version
matched_negative_feature
matched_negative_query_abs
matched_negative_background_status
matched_negative_scope
matched_negative_n
matched_negative_abs_percentile
matched_negative_abs_median
matched_negative_abs_p95
matched_negative_action
```

## Background Fallback

The fallback order is:

1. same region in eligible N0 samples with enough background rows
2. same chromosome and similar length
3. autosome and similar length
4. `UNKNOWN_BACKGROUND`

`UNKNOWN_BACKGROUND` always maps to:

```text
matched_negative_action=REVIEW_NO_CALL
```

It must not become a hard artifact rule.

## Current Seed Result

The current negative-bank seed has:

```text
N0: 0
N1: 6
N2: 4
matched_negative_eligible: 0
```

Therefore the current expected Phase 3 result is:

```text
all Y1-Y8 candidate rows -> UNKNOWN_BACKGROUND
all Y1-Y8 candidate rows -> REVIEW_NO_CALL
```

This is the correct safe behavior because H7-H16 are not promoted to N0.

This statement applies only to Branch B matched-negative empirical-null use.
It should not be used to exclude H7-H16 from reference-rebuild exploration.
Reference rebuild eligibility is handled by separate `R0/R1/R2` labels, and
old-reference Branch A signal can represent reference mismatch or batch shift
rather than a true individual positive.

## Remote Validation

RED checks:

```bash
ssh fengxian "cd /data/project/CNV/PGT-A/refactor_validation_20260419 && PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 /biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest tests/unit/test_branch_b_matched_negative.py -q"
```

Initial result:

```text
ModuleNotFoundError: No module named 'pgta.predict.branch_b.matched_negative'
```

Workflow contract RED check:

```bash
ssh fengxian "cd /data/project/CNV/PGT-A/refactor_validation_20260419 && PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 /biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest tests/unit/test_branch_ab_phase12_workflow_contract.py::test_phase3_matched_negative_is_shadow_only_and_report_safe -q"
```

Initial result:

```text
failed because dispatcher/workflow did not expose matched_negative_percentile
```

GREEN check:

```bash
ssh fengxian "cd /data/project/CNV/PGT-A/refactor_validation_20260419 && PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 /biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest tests/unit/test_branch_b_matched_negative.py tests/unit/test_branch_ab_phase12_workflow_contract.py -q"
```

Result:

```text
7 passed in 0.57s
```

Dry-run:

```bash
ssh fengxian "cd /data/project/CNV/PGT-A/refactor_validation_20260419 && /biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile /data/project/CNV/PGT-A/refactor_validation_20260419/config_predict_build_ref_v2_mask_only.yaml --cores 1 -n cnv"
```

Result:

```text
cnv_branch_b_matched_negative: 8 jobs
```

Execution:

```bash
ssh fengxian "cd /data/project/CNV/PGT-A/refactor_validation_20260419 && /biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile /data/project/CNV/PGT-A/refactor_validation_20260419/config_predict_build_ref_v2_mask_only.yaml --cores 4 cnv"
```

Result:

```text
8 matched-negative jobs completed.
Complete log: .snakemake/log/2026-06-18T164037.257699.snakemake.log
```

Output summary:

| sample | rows | background status | action |
|---|---:|---|---|
| Y1 | 17 | UNKNOWN_BACKGROUND=17 | REVIEW_NO_CALL=17 |
| Y2 | 22 | UNKNOWN_BACKGROUND=22 | REVIEW_NO_CALL=22 |
| Y3 | 15 | UNKNOWN_BACKGROUND=15 | REVIEW_NO_CALL=15 |
| Y4 | 25 | UNKNOWN_BACKGROUND=25 | REVIEW_NO_CALL=25 |
| Y5 | 9 | UNKNOWN_BACKGROUND=9 | REVIEW_NO_CALL=9 |
| Y6 | 15 | UNKNOWN_BACKGROUND=15 | REVIEW_NO_CALL=15 |
| Y7 | 19 | UNKNOWN_BACKGROUND=19 | REVIEW_NO_CALL=19 |
| Y8 | 21 | UNKNOWN_BACKGROUND=21 | REVIEW_NO_CALL=21 |

## Interpretation

Phase 3 is now implemented as an auditable shadow feature path.

The current run does not produce usable matched-negative percentiles because there are no N0 samples. This is intentional. The code preserves the safe behavior required by the plan: unknown background cannot be used as clean support and cannot suppress candidates as hard artifacts.

## Remaining Gate

To make Phase 3 informative instead of `UNKNOWN_BACKGROUND`, at least one
locked clean N0 cohort must be defined and processed through the same
evidence-ledger path. H7-H16 should not be used for matched-negative empirical
null unless a post-rebuild or holdout/ablation run promotes specific samples to
N0.

This does not block a reference-rebuild experiment. A shadow rebuild may include
documented `R0/R1` candidates, then rerun WisecondorX predict, Branch A, Branch
B, evaluation, benchmark, and report to decide whether the expanded reference
improves FP burden without increasing known-positive FN.

If H7-H16 are added to a named shadow reference, the current `N0=0` result is no
longer the operative Phase 3 decision. The post-rebuild workflow must create a
new negative-bank version tied to the new reference ID and recompute
`N0/N1/N2` from the regenerated WisecondorX predict, Branch A candidates, and
Branch B evidence. The report must keep separate whether each sample was:

- included in the reference;
- held out as negative evidence;
- used only for shadow review.

Phase 3 remains blocked only for percentile-based filtering until this
post-rebuild negative-bank decision exists. It does not block R0/R1 shadow
reference construction.
