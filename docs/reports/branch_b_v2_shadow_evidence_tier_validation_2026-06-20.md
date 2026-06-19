# Branch B V2 Shadow Evidence Tier Validation

Date: 2026-06-20

## Scope

This report records the first implementation pass for Branch B V2 shadow evidence tiers.

Branch B V2 remains shadow-only:

- `v2_final_report_impact = none_shadow_only`
- no Branch B V2 standalone final candidate
- no V2 hard filter or promotion into final report
- no matched-negative-derived hard decision when the background is missing

The implementation only adds candidate-level evidence labels to help review Branch A candidates while the R0/R1/R2 and N0/N1/N2 validation design is still being developed.

## Context Documents Read

- `docs/handoff/2026-06-20_0005_branch_b_v2_n0_reframe_handoff.md`
- `AGENTS.md`
- `skills/conversation_handoff/SKILL.md`
- `skills/pgta_reference_modeling_analysis/SKILL.md`
- `REPO_MAP.md`
- `PLANS.md`
- `CURRENT_STATE.md`
- `docs/reports/branch_ab_v2_rnd_constraints_2026-06-18.md`
- `docs/reports/h_r0_shadow_post_rebuild_predict_2026-06-19.md`
- `docs/reports/h7_h16_reference_rebuild_eligibility_2026-06-19.tsv`
- `docs/reports/branch_ab_v2_phase4_unknown_background_fix_2026-06-19.md`

## Code Changes

Changed files:

- `pgta/predict/branch_b/v2_classifier.py`
- `tests/unit/test_branch_b_v2_classifier.py`

New classifier output columns:

- `v2_evidence_tier`
- `v2_evidence_gate`
- `v2_review_priority`

Summary JSON now includes:

- `evidence_tier_counts`
- `review_priority_counts`

## Evidence Tier Semantics

The current tiering is intentionally conservative.

`UNKNOWN_BACKGROUND_POSITIVE_SUPPORT` means:

- Branch A has positive evidence, but
- matched-negative background is not validated, so
- the event remains review-only.

`MATCHED_NEGATIVE_OUTLIER_POSITIVE_SUPPORT` means:

- matched-negative background exists,
- the candidate is an outlier relative to that background,
- Branch A direction/support is positive,
- but the event still remains shadow-only until promotion criteria are satisfied.

`UNKNOWN_BACKGROUND_DIRECTION_CONFLICT` means:

- candidate state and amplitude direction conflict under the available evidence,
- no final hard decision is made by V2.

`UNKNOWN_BACKGROUND_REF_CONTRACT_RISK` means:

- a reference-contract risk such as same-chromosome leakage was detected,
- this is a no-call contract risk in V2 shadow output.

## Remote Unit Test Evidence

Executable:

`/biosoftware/miniconda/envs/snakemake_env/bin/python`

Command:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 /biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest tests/unit/test_branch_b_v2_classifier.py -q'
```

Result:

```text
8 passed in 0.51s
```

Regression command:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 /biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest tests/unit/test_branch_a_candidates.py tests/unit/test_branch_b_evidence_ledger.py tests/unit/test_branch_b_matched_negative.py tests/unit/test_branch_b_v2_classifier.py tests/unit/test_negative_bank.py tests/unit/test_branch_ab_phase12_workflow_contract.py -q'
```

Result:

```text
30 passed in 0.74s
```

## Remote Dry-Run Evidence

Executable:

`/biosoftware/miniconda/envs/snakemake_env/bin/snakemake`

Command:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && for cfg in config_predict_y_h_r0_shadow_20260619.yaml config_predict_h_20260608_h_r0_shadow_20260619.yaml config_predict_20260615_h_r0_shadow_20260619.yaml; do /biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile "$cfg" --cores 1 --rerun-triggers mtime --forcerun cnv_branch_b_v2_classifier -n cnv; done'
```

Result:

- Y config planned 8 `cnv_branch_b_v2_classifier` jobs.
- H 2026-06-08 config planned 16 `cnv_branch_b_v2_classifier` jobs.
- 2026-06-15 config planned 5 `cnv_branch_b_v2_classifier` jobs.
- No BAM regeneration was planned.
- No WisecondorX convert or predict rerun was planned.

## Remote Classifier Refresh Evidence

Command:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && for cfg in config_predict_y_h_r0_shadow_20260619.yaml config_predict_h_20260608_h_r0_shadow_20260619.yaml config_predict_20260615_h_r0_shadow_20260619.yaml; do /biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile "$cfg" --cores 16 --rerun-triggers mtime --forcerun cnv_branch_b_v2_classifier cnv; done'
```

Result:

- completed successfully
- latest visible Snakemake log for 2026-06-15: `.snakemake/log/2026-06-20T002031.606922.snakemake.log`
- 2026-06-15 report files exist:
  - `results_20260615_h_r0_shadow_ref_20260619/wisecondorx/cnv/report/cnv_summary.tsv`
  - `results_20260615_h_r0_shadow_ref_20260619/wisecondorx/cnv/report/cnv_summary.json`
  - `results_20260615_h_r0_shadow_ref_20260619/wisecondorx/cnv/report/cnv_summary.md`
  - `results_20260615_h_r0_shadow_ref_20260619/wisecondorx/cnv/report/cnv_summary.html`

Log audit:

| config | classifier jobs | map reads | convert | predict |
|---|---:|---:|---:|---:|
| Y H-R0 shadow | 8 | 0 | 0 | 0 |
| H 2026-06-08 H-R0 shadow | 16 | 0 | 0 | 0 |
| 2026-06-15 H-R0 shadow | 5 | 0 | 0 | 0 |

## V2 Shadow Output Summary

| cohort | files | candidates | v2 final impact | high priority | medium priority |
|---|---:|---:|---|---:|---:|
| Y1-Y8 | 8 | 131 | none_shadow_only | 119 | 12 |
| H1-H16 | 16 | 221 | none_shadow_only | 211 | 10 |
| 2026-06-15 five samples | 5 | 201 | none_shadow_only | 170 | 31 |

Evidence tier counts:

| cohort | UNKNOWN_BACKGROUND_POSITIVE_SUPPORT | UNKNOWN_BACKGROUND_DIRECTION_CONFLICT |
|---|---:|---:|
| Y1-Y8 | 119 | 12 |
| H1-H16 | 211 | 10 |
| 2026-06-15 five samples | 170 | 31 |

This result confirms that V2 currently provides evidence prioritization only. It does not reduce or add final report events.

## 2026-06-15 Current Report Summary

Source:

`/data/project/CNV/PGT-A/refactor_validation_20260419/results_20260615_h_r0_shadow_ref_20260619/wisecondorx/cnv/report/cnv_summary.tsv`

| sample | QC | sex | Branch A candidates | Branch B kept | reportable | top report event |
|---|---|---|---:|---:|---:|---|
| JZ26125843-56-56 | PASS | XX | 28 | 0 | 0 | none |
| JZ26125844-59-59 | PASS | XX | 38 | 3 | 3 | chr4:52500001-186000000 gain |
| JZ26125845-60-60 | PASS | XY | 57 | 26 | 26 | chr14:60750001-106500000 gain |
| JZ26125846-61-61 | PASS | XY | 44 | 3 | 3 | chr12:37500001-124500000 gain |
| JZ26125847-62-62 | PASS | XY | 34 | 1 | 1 | chr4:121500001-186000000 loss |

Interpretation:

- 56 has Branch A signals but no Branch B reportable event under the current report logic.
- 59, 61, and 62 have low-count review/reportable events.
- 60 remains high-burden and needs separate review before it can be treated like the other four.
- V2 evidence-tier output did not change this report.

## Constraints Still Active

- H7-H16 must not be defined as N0 only by the old Branch B behavior.
- Adding H7-H16 or any inlier samples into a reference means WisecondorX predict, Branch A candidates, Branch B, evaluation, benchmark, and report must all be rerun.
- Branch A remains the WisecondorX predict/CBS high-sensitivity candidate layer.
- Branch B remains a downstream candidate refinement layer.
- V2 cannot hard-filter or promote final calls until matched-negative/cross-fit validation is available.

## Next Gate

The next implementation step should not be a hard filter. It should be:

1. build cross-fit or held-out negative background using the R0/R1/R2 reference plan;
2. rerun WisecondorX predict and Branch A for each reference candidate set;
3. evaluate positive recall on Y1-Y8 and H1-H6;
4. evaluate negative review burden on held-out H7-H16 and later true-negative batches;
5. only then decide whether V2 tiers can affect report disposition.

