# 2026-06-20 01:05 Branch B V2 N1 Shadow Background Handoff

## Source Context

Latest handoff read before implementation:

`docs/handoff/2026-06-20_0035_branch_b_v2_shadow_evidence_handoff.md`

Other required context:

- `AGENTS.md`
- `skills/conversation_handoff/SKILL.md`
- `skills/pgta_reference_modeling_analysis/SKILL.md`
- `docs/reports/branch_ab_v2_rnd_constraints_2026-06-18.md`
- `docs/reports/branch_ab_v2_phase4_unknown_background_fix_2026-06-19.md`
- `docs/reports/h7_h16_reference_rebuild_eligibility_2026-06-19.tsv`

## Branch

Local branch:

`codex/branch-b-v2-shadow-background`

## Completed

Implemented N1 shadow background for Branch B V2 matched-negative evidence.

Key behavior:

- N0 remains the only formal empirical null.
- Configured N1 labels can provide `SHADOW_BACKGROUND`.
- Shadow background uses `SHADOW_CONTEXT_ONLY`.
- V2 classifier emits shadow evidence tiers but keeps `none_shadow_only`.

Changed files:

- `pgta/predict/branch_b/matched_negative.py`
- `pgta/predict/branch_b/v2_classifier.py`
- `rules/predict_layout.smk`
- `rules/predict_workflow.smk`
- `config_predict_y_h_r0_shadow_20260619.yaml`
- `config_predict_h_20260608_h_r0_shadow_20260619.yaml`
- `config_predict_20260615_h_r0_shadow_20260619.yaml`
- `tests/unit/test_branch_b_matched_negative.py`
- `tests/unit/test_branch_b_v2_classifier.py`
- `tests/unit/test_branch_ab_phase12_workflow_contract.py`
- `docs/reports/branch_b_v2_n1_shadow_background_validation_2026-06-20.md`
- `docs/handoff/2026-06-20_0105_branch_b_v2_n1_shadow_background_handoff.md`

## Remote Validation

Remote path:

`/data/project/CNV/PGT-A/refactor_validation_20260419`

Executable:

`/biosoftware/miniconda/envs/snakemake_env/bin/python`

Test result:

```text
33 passed in 0.69s
```

The final regression includes the pandas NA fallback edge case for V2 classifier:

- blank `matched_negative_background_status` plus `matched_negative_source=UNKNOWN_BACKGROUND` remains review-only;
- `SHADOW_BACKGROUND` takes precedence over stale `matched_negative_source=UNKNOWN_BACKGROUND`.

Snakemake executable:

`/biosoftware/miniconda/envs/snakemake_env/bin/snakemake`

Dry-run and actual refresh were run on:

- `config_predict_y_h_r0_shadow_20260619.yaml`
- `config_predict_h_20260608_h_r0_shadow_20260619.yaml`
- `config_predict_20260615_h_r0_shadow_20260619.yaml`

No BAM, convert, or WisecondorX predict jobs were rerun in the final classifier refresh.

## Materialized Result Summary

| cohort | candidates | shadow background | unknown background |
|---|---:|---:|---:|
| Y1-Y8 | 131 | 38 | 93 |
| H1-H16 | 221 | 145 | 76 |
| 2026-06-15 five samples | 201 | 32 | 169 |

2026-06-15 current report:

| sample | QC | sex | Branch A candidates | Branch B kept | reportable |
|---|---|---|---:|---:|---:|
| JZ26125843-56-56 | PASS | XX | 28 | 0 | 0 |
| JZ26125844-59-59 | PASS | XX | 38 | 3 | 3 |
| JZ26125845-60-60 | PASS | XY | 57 | 26 | 26 |
| JZ26125846-61-61 | PASS | XY | 44 | 3 | 3 |
| JZ26125847-62-62 | PASS | XY | 34 | 1 | 1 |

## Important Caveat

The 2026-06-15 report remains driven by the current Branch A+B report chain.

N1 shadow background does not change final report events. It only supplies review evidence.

## Next Step

The next meaningful gate is not another V2 hard filter. It is a validation design:

1. build held-out N0 or cross-fit negative background;
2. rerun WisecondorX predict and Branch A after every reference change;
3. verify Y1-Y8 and H1-H6 recall;
4. measure held-out negative review burden;
5. only then decide whether any V2 evidence tier may affect final report disposition.
