# Branch A/B V2 Phase 5: Branch S Shadow Evidence

## Scope

Phase 5 adds a Branch S shadow path for sex-chromosome evidence only.

It does not replace:

- WisecondorX predict.
- Branch A candidate generation.
- Legacy Branch B final events.
- Current sex routing.
- `cnv_report`.

## Implementation

New module:

- `pgta/predict/branch_s.py`

Workflow outputs:

- `wisecondorx/cnv/postprocess/branch_s/{sample}.sex_chrom_evidence.tsv`
- `wisecondorx/cnv/postprocess/branch_s/{sample}.sca_state_scores.tsv`
- `wisecondorx/cnv/postprocess/branch_s/{sample}.summary.json`

Branch S consumes:

- Branch B calibrated bins.
- Branch A candidates.
- Gender TSV when sex routing is enabled.

Branch S emits:

- X non-PAR evidence.
- X PAR evidence when available.
- Y non-PAR evidence.
- Y PAR evidence when available.
- Shadow SCA state scores for `X_GAIN`, `X_LOSS`, `Y_GAIN`, and `Y_LOSS`.

All rows keep:

- `branch_s_action=SHADOW_ONLY`
- `final_report_impact=none_shadow_only`

Summary metadata keeps:

- `replaces_current_sex_calling=false`
- `replaces_final_report=false`

## Validation

Environment:

- remote: `fengxian`
- code path: `/data/project/CNV/PGT-A/refactor_validation_20260419`
- python: `/biosoftware/miniconda/envs/snakemake_env/bin/python`
- snakemake: `/biosoftware/miniconda/envs/snakemake_env/bin/snakemake`

RED:

```text
tests/unit/test_branch_s_shadow.py failed with:
ModuleNotFoundError: No module named 'pgta.predict.branch_s'
```

GREEN subset:

```bash
ssh fengxian "cd /data/project/CNV/PGT-A/refactor_validation_20260419 && PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 /biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest tests/unit/test_branch_s_shadow.py tests/unit/test_branch_ab_phase12_workflow_contract.py::test_phase5_branch_s_is_shadow_only_and_report_safe -q"
```

Result:

```text
3 passed in 0.57s
```

Regression:

```bash
ssh fengxian "cd /data/project/CNV/PGT-A/refactor_validation_20260419 && PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 /biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest tests/unit/test_branch_ab_phase12_workflow_contract.py tests/unit/test_branch_s_shadow.py tests/unit/test_branch_b_v2_classifier.py tests/unit/test_branch_b_matched_negative.py tests/unit/test_branch_b_evidence_ledger.py tests/unit/test_negative_bank.py tests/unit/test_branch_b_calling.py tests/unit/test_branch_b_artifact_rules.py tests/unit/test_branch_b_correction.py -q"
```

Result:

```text
81 passed in 1.01s
```

Dry-run after implementation:

```bash
ssh fengxian "cd /data/project/CNV/PGT-A/refactor_validation_20260419 && /biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile /data/project/CNV/PGT-A/refactor_validation_20260419/config_predict_build_ref_v2_mask_only.yaml --cores 1 -n cnv"
```

Result before execution:

```text
Job stats:
job                    count
-------------------  -------
cnv                        1
cnv_branch_s_shadow        8
total                      9
```

Execution:

```bash
ssh fengxian "cd /data/project/CNV/PGT-A/refactor_validation_20260419 && /biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile /data/project/CNV/PGT-A/refactor_validation_20260419/config_predict_build_ref_v2_mask_only.yaml --cores 4 cnv"
```

Result:

```text
8 cnv_branch_s_shadow jobs completed.
Complete log: .snakemake/log/2026-06-18T221547.235491.snakemake.log
```

Final dry-run:

```text
Nothing to be done (all requested files are present and up to date).
```

Output check:

```text
Y1 True True 3 4 none_shadow_only False
Y2 True True 3 4 none_shadow_only False
Y3 True True 3 4 none_shadow_only False
Y4 True True 3 4 none_shadow_only False
Y5 True True 3 4 none_shadow_only False
Y6 True True 3 4 none_shadow_only False
Y7 True True 3 4 none_shadow_only False
Y8 True True 3 4 none_shadow_only False
```

## Current Conclusion

Phase 5 is implemented as a shadow-only Branch S evidence path.

It is safe as a workflow extension because the outputs are part of `cnv` for materialization but are not consumed by `cnv_report`, sex routing, or legacy Branch B final-event logic.

## Next Step

Use the Branch S tables for SCA review design and truth-set planning. Do not promote SCA state scores to reportable calls until there is a separately locked SCA validation set.

Follow-up validation plan:

- `docs/reports/branch_s_sca_validation_plan_2026-06-19.md`

Current limitation:

- Branch S output has been materialized for Y1-Y8 only.
- H1-H6, H7-H16, and the 2026-06-15 exploratory samples still need formal Branch S materialization before SCA performance can be evaluated.
