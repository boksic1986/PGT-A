# Branch B/S Lowres Evidence Integration

Date: 2026-06-22

## Scope

This loop connects the completed 2Mb/3Mb shadow references to the Branch B/Branch S
development workflow as auxiliary evidence only.

Unchanged boundaries:

- 1Mb Branch A remains the primary WisecondorX/CBS discovery layer.
- `merge_gap_bp=2_000_000` remains an explicit overlay, not the default.
- `h_r0_shadow_ref_20260619` remains the active 1Mb shadow reference, not a
  production reference.
- 2Mb/3Mb absence is not a standalone demotion/filter.
- Branch S remains `review_development_only`, not final SCA.

## Code Changes

- `rules/predict_layout.smk`
  - `lowres_evidence.reference_npz_glob` is now accepted and expanded into the
    existing ref-MAD input list.
  - Explicit `reference_npz` remains supported.
  - Lowres-enabled configs no longer fail parse when they provide
    `reference_npz_glob`.
- `rules/predict_workflow.smk`
  - Branch S receives optional `--lowres-2mb-events` and
    `--lowres-3mb-events` only when lowres evidence is enabled.
- `pgta/predict/branch_s.py`
  - Added segment-level X/Y non-PAR evidence for Branch A sex-chromosome
    candidates.
  - Whole-chromosome X/Y non-PAR median remains the global SCA trend evidence.
  - Local segment non-PAR median/mean can preserve small sex-chromosome CNV
    review when global non-PAR median is neutral.
  - PAR evidence remains secondary context and does not independently create
    XO/XXY/XYY calls.
  - Added lowres sex-chromosome context fields:
    `sex_chrom_lowres_2mb_context`,
    `sex_chrom_lowres_2mb_same_direction_count`,
    `sex_chrom_lowres_3mb_context`,
    `sex_chrom_lowres_3mb_same_direction_count`,
    `sex_chrom_lowres_final_impact=context_only_not_filter`.

## Configs Added To The Branch

- `config_predict_y_h_r0_shadow_2mb_branch_a_gap2m_20260622.yaml`
- `config_predict_y_h_r0_shadow_3mb_branch_a_gap2m_20260622.yaml`
- `config_predict_h_20260608_h_r0_shadow_2mb_branch_a_gap2m_20260622.yaml`
- `config_predict_h_20260608_h_r0_shadow_3mb_branch_a_gap2m_20260622.yaml`
- `config_predict_g1_g8_h_r0_shadow_2mb_branch_a_gap2m_20260622.yaml`
- `config_predict_g1_g8_h_r0_shadow_3mb_branch_a_gap2m_20260622.yaml`
- `config_predict_20260615_h_r0_shadow_2mb_branch_a_gap2m_20260622.yaml`
- `config_predict_20260615_h_r0_shadow_3mb_branch_a_gap2m_20260622.yaml`
- `config_predict_y_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`
- `config_predict_h_20260608_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`
- `config_predict_g1_g8_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`
- `config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`

## Remote State

2Mb/3Mb shadow references are complete:

- `results_h_r0_shadow_ref_20260619_2mb/reference/...`
- `results_h_r0_shadow_ref_20260619_3mb/reference/...`
- `common_best_binsize.txt`: `2000000` and `3000000`.

Lowres predict materialization was started in the remote mirror:

- PID file:
  `/data/project/CNV/PGT-A/refactor_validation_20260419/logs/driver/lowres_predict_2mb_3mb_20260622.pid`
- PID at launch: `4690`
- log:
  `/data/project/CNV/PGT-A/refactor_validation_20260419/logs/driver/lowres_predict_2mb_3mb_20260622.log`
- command file:
  `/data/project/CNV/PGT-A/refactor_validation_20260419/logs/driver/lowres_predict_2mb_3mb_20260622.command.sh`
- runtime database:
  `monitor/runtime.db` was absent, so no historical duration estimate is
  available.

Final observed progress during this implementation:

- Y 2Mb/3Mb lowres predict completed: 8/8 predict done, 8/8 A-branch events,
  8/8 V2 classifications, 8/8 Branch S summaries at each binsize.
- H 2Mb/3Mb lowres predict completed: 16/16 at each output level and binsize.
- G 2Mb/3Mb lowres predict completed: 8/8 at each output level and binsize.
- 2026-06-15 2Mb/3Mb lowres predict completed: 5/5 at each output level and
  binsize.

## Remote Verification Completed

Unit tests on `ssh fengxian`:

```text
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_s_shadow.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py \
  tests/unit/test_lowres_evidence.py \
  tests/unit/test_ref_stability.py \
  tests/unit/test_branch_b_v2_classifier.py -q
```

Result:

```text
68 passed in 1.53s
```

Dry-run status:

- 2Mb/3Mb predict configs parsed and planned successfully for Y/H/G/2026-06-15.
- Job stats used existing BAM inputs and did not include BWA/mapping rules.
- Lowres-enabled main-chain dry-runs passed for all four active configs:
  Y, H, G, and 2026-06-15.
- Lowres-enabled main-chain materialization completed for all four active
  configs.

## Materialized Acceptance

| cohort | samples | candidates | truth | preserved | filtered truth | report | internal review | filtered audit | Branch S | status |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| Y1-Y8 | 8 | 97 | 10 | 10 | 0 | 21 | 50 | 13 | 13 | ready |
| H1-H16 | 16 | 105 | 10 | 10 | 0 | 6 | 49 | 8 | 42 | ready |
| G1-G8 | 8 | 75 | 10 | 10 | 0 | 15 | 40 | 7 | 13 | ready |
| 2026-06-15 | 5 | 165 | 0 | 0 | 0 | 52 | 76 | 23 | 14 | skipped_no_truth |

H6 chr21 remains visible:

```text
sample=H6
chrom=chr21
start=15000001
end=42000000
state=gain
v2_report_layer_class=internal_review_event
v2_filter_action=keep_background_unknown_review
v2_lowres_context_label=LOWRES_2MB_3MB_SAME_DIRECTION_SUPPORT
v2_ref_stability_context=REF_STABILITY_STABLE
```

Example Branch S lowres context is present in sample summaries and remains
`context_only_not_filter`.

## Interpretation

Branch B lowres/ref-MAD evidence remains an auxiliary candidate-level context
layer. Same-direction lowres support can increase confidence and improve report
explanation. Lowres absence is not a standalone reason to demote or filter,
especially for:

- H6 chr21-like weak positives;
- 3-4Mb or shorter boundary-diluted events;
- small sex-chromosome CNV candidates;
- mosaic-sensitive events.

Branch S now separates global dosage evidence from local segment evidence:

- Global X/Y non-PAR median supports whole-chromosome SCA trends such as XO,
  XXY, or XYY.
- Local segment X/Y non-PAR median/mean protects small X/Y CNV candidates from
  being hidden by neutral whole-chromosome medians.
- PAR/XY homology context can explain whole sex-chromosome dosage consistency,
  but cannot independently create a final SCA call.

## Pending Gate

The next gate is interpretation, not another workflow contract fix:

1. Inspect lowres support labels per truth event and per remaining report event.
2. Use lowres same-direction support to improve report explanation.
3. Do not use lowres absence as standalone filtering evidence.
4. Review high ref-MAD regions before interpreting lowres absence.
5. Continue Branch S as review/development-only until a broader SCA truth panel
   is available.
