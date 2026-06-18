# Branch A/B V2 Phase 4 Implementation Report

Date: 2026-06-18

## Scope

This report implements Phase 4 from `docs/reports/branch_ab_v2_rnd_constraints_2026-06-18.md`.

Phase 4 adds a Branch B V2 candidate-level classifier in shadow mode. It does not change:

- WisecondorX predict
- Branch A candidate inclusion
- legacy Branch B final events
- `cnv_report`
- artifact hard filtering

## Phase 4 Contract

Input evidence:

- Phase 3 matched-negative candidate evidence when `negative_bank.samples_tsv` is configured
- otherwise Phase 1 evidence ledger

Output paths:

```text
wisecondorx/cnv/postprocess/v2_classifier/{sample}.candidate_classification.tsv
wisecondorx/cnv/postprocess/v2_classifier/{sample}.summary.json
```

Implementation:

```text
pgta/predict/branch_b/v2_classifier.py
scripts/predict.py branch_b_v2_classifier
rule cnv_branch_b_v2_classifier
```

New shadow columns:

```text
v2_classifier_version
v2_candidate_class
v2_classifier_action
v2_classifier_reason
v2_final_report_impact
```

The classifier keeps one output row per Branch A candidate and preserves the legacy Branch B comparator fields, including `final_disposition`, `branch_b_report_class`, `branch_b_artifact_status`, and `branch_b_keep_event`.

## Safety Behavior

`UNKNOWN_BACKGROUND` from Phase 3 forces:

```text
v2_candidate_class=REVIEW_REQUIRED
v2_classifier_action=SHADOW_REVIEW_ONLY
v2_classifier_reason=unknown_matched_negative_background
```

It is not converted to `LIKELY_ARTIFACT` and is not consumed by `cnv_report`.

## Remote Validation

Remote environment:

```text
host: fengxian
path: /data/project/CNV/PGT-A/refactor_validation_20260419
python: /biosoftware/miniconda/envs/snakemake_env/bin/python
snakemake: /biosoftware/miniconda/envs/snakemake_env/bin/snakemake
config: /data/project/CNV/PGT-A/refactor_validation_20260419/config_predict_build_ref_v2_mask_only.yaml
```

RED checks:

```bash
ssh fengxian "cd /data/project/CNV/PGT-A/refactor_validation_20260419 && PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 /biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest tests/unit/test_branch_b_v2_classifier.py tests/unit/test_branch_ab_phase12_workflow_contract.py::test_phase4_v2_classifier_is_shadow_only_and_report_safe -q"
```

Initial result:

```text
ModuleNotFoundError: No module named 'pgta.predict.branch_b.v2_classifier'
```

Workflow contract RED check:

```bash
ssh fengxian "cd /data/project/CNV/PGT-A/refactor_validation_20260419 && PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 /biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest tests/unit/test_branch_ab_phase12_workflow_contract.py::test_phase4_v2_classifier_is_shadow_only_and_report_safe -q"
```

Initial result:

```text
failed because dispatcher/workflow did not expose branch_b_v2_classifier
```

GREEN check:

```bash
ssh fengxian "cd /data/project/CNV/PGT-A/refactor_validation_20260419 && PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 /biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest tests/unit/test_branch_b_v2_classifier.py tests/unit/test_branch_ab_phase12_workflow_contract.py -q"
```

Result:

```text
9 passed in 0.68s
```

Regression check:

```bash
ssh fengxian "cd /data/project/CNV/PGT-A/refactor_validation_20260419 && PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 /biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest tests/unit/test_branch_ab_phase12_workflow_contract.py tests/unit/test_branch_b_v2_classifier.py tests/unit/test_branch_b_matched_negative.py tests/unit/test_branch_b_evidence_ledger.py tests/unit/test_negative_bank.py tests/unit/test_branch_b_calling.py tests/unit/test_branch_b_artifact_rules.py tests/unit/test_branch_b_correction.py -q"
```

Result:

```text
78 passed in 0.90s
```

Dry-run:

```bash
ssh fengxian "cd /data/project/CNV/PGT-A/refactor_validation_20260419 && /biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile /data/project/CNV/PGT-A/refactor_validation_20260419/config_predict_build_ref_v2_mask_only.yaml --cores 1 -n cnv"
```

Result before execution:

```text
cnv_branch_b_v2_classifier: 8 jobs
cnv: 1 job
```

Execution:

```bash
ssh fengxian "cd /data/project/CNV/PGT-A/refactor_validation_20260419 && /biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile /data/project/CNV/PGT-A/refactor_validation_20260419/config_predict_build_ref_v2_mask_only.yaml --cores 4 cnv"
```

Result:

```text
8 cnv_branch_b_v2_classifier jobs completed.
Complete log: .snakemake/log/2026-06-18T215036.590993.snakemake.log
```

Final dry-run:

```text
Nothing to be done (all requested files are present and up to date).
```

## Current Output Summary

Because the current negative-bank seed has no N0 samples, Phase 3 provides only `UNKNOWN_BACKGROUND`. Phase 4 therefore keeps all current Y1-Y8 candidates in shadow review:

| sample | candidate rows | v2 class counts | final report impact |
|---|---:|---|---|
| Y1 | 17 | REVIEW_REQUIRED=17 | none_shadow_only |
| Y2 | 22 | REVIEW_REQUIRED=22 | none_shadow_only |
| Y3 | 15 | REVIEW_REQUIRED=15 | none_shadow_only |
| Y4 | 25 | REVIEW_REQUIRED=25 | none_shadow_only |
| Y5 | 9 | REVIEW_REQUIRED=9 | none_shadow_only |
| Y6 | 15 | REVIEW_REQUIRED=15 | none_shadow_only |
| Y7 | 19 | REVIEW_REQUIRED=19 | none_shadow_only |
| Y8 | 21 | REVIEW_REQUIRED=21 | none_shadow_only |

## Interpretation

Phase 4 is now implemented as a shadow classifier path. The output is intentionally conservative because the current project state has no locked clean N0 cohort. It gives a stable place to attach future candidate-level classification evidence without adding more hard filters to legacy Branch B.

## Remaining Gate

Do not promote the V2 classifier to final reporting until ablation has been run against Y1-Y8, H1-H6, H7-H16, and the 2026-06-15 shadow-report set, with explicit confirmation that known-positive recall does not regress.
