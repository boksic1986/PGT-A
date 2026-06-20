# Branch A Burden Optimization Phase 1

Date: 2026-06-21

Status: `active_current_evidence`

Decision use: records the first Branch A burden optimization pass under
`h_r0_shadow_ref_20260619`. It supports config-gated Branch A merge-gap
evaluation. It does not promote a new production report route.

## Scope

This phase keeps Branch A as WisecondorX predict/CBS-derived candidate
discovery.

Unchanged:

- WisecondorX predict remains the primary Branch A caller.
- Current predict parameters remain `zscore=5`, `alpha=0.001`,
  `maskrepeats=5`, `minrefbins=150`.
- `strong_z=10` remains a support tier label only, not an inclusion hard
  filter.
- No Branch B artifact, calibration, matched-negative, or legacy disposition
  field is introduced into Branch A.
- No CNV calling, mosaic, sex calling, or final report schema is changed.

## Current Baseline

Current P2 Branch A validation after the workflow parameter plumbing was
rematerialized remotely with default `merge_gap_bp=0`.

| cohort | samples | truth events | detected | FN | Branch A candidates | key note |
|---|---:|---:|---:|---:|---:|---|
| Y1-Y8 | 8 | 10 | 10 | 0 | 131 | no-FN retained |
| H1-H16 | 16 | 10 | 10 | 0 | 221 | H6 chr21 retained |
| 2026-06-15 | 5 | 0 | 0 | 0 | 201 | burden context only; no truth |

Per-sample Branch A candidate summary JSON now records:

```text
merge_gap_bp=0
strong_z=10.0
```

## Z Threshold Ablation

This is a post-hoc read of existing Branch A candidates. It is not a WisecondorX
predict rerun.

| cohort | threshold | candidates kept | truth detected | FN | missed |
|---|---:|---:|---:|---:|---|
| Y1-Y8 | 5 | 131 | 10/10 | 0 | none |
| Y1-Y8 | 6 | 112 | 10/10 | 0 | none |
| Y1-Y8 | 7 | 98 | 10/10 | 0 | none |
| Y1-Y8 | 8 | 88 | 10/10 | 0 | none |
| Y1-Y8 | 10 | 72 | 9/10 | 1 | Y6 chr7 loss |
| H1-H16 | 5 | 221 | 10/10 | 0 | none |
| H1-H16 | 6 | 194 | 10/10 | 0 | none |
| H1-H16 | 7 | 183 | 10/10 | 0 | none |
| H1-H16 | 8 | 177 | 9/10 | 1 | H6 chr21 gain |
| H1-H16 | 10 | 170 | 9/10 | 1 | H6 chr21 gain |
| 2026-06-15 | 5 | 201 | 0/0 | 0 | no truth |
| 2026-06-15 | 8 | 138 | 0/0 | 0 | no truth |
| 2026-06-15 | 10 | 120 | 0/0 | 0 | no truth |

Interpretation:

- Raising Branch A hard z above 7 is unsafe for the current truth set.
- `z>=8` would drop H6 chr21 gain.
- `z>=10` would drop Y6 chr7 loss and H6 chr21 gain.
- `strong_z=10` must remain a tier label for strong evidence, not a hard
  Branch A inclusion threshold.

## Merge-Gap Ablation

This ablation reassembled existing WisecondorX aberration BED files with
different Branch A same-sample/same-chromosome/same-direction merge gaps. It
did not rerun WisecondorX predict.

| cohort | merge gap | candidates | truth detected | FN | H6 chr21 |
|---|---:|---:|---:|---:|---|
| Y1-Y8 | 0 bp | 131 | 10/10 | 0 | NA |
| Y1-Y8 | 1 Mb | 131 | 10/10 | 0 | NA |
| Y1-Y8 | 2 Mb | 97 | 10/10 | 0 | NA |
| Y1-Y8 | 4 Mb | 87 | 10/10 | 0 | NA |
| H1-H16 | 0 bp | 221 | 10/10 | 0 | detected |
| H1-H16 | 1 Mb | 220 | 10/10 | 0 | detected |
| H1-H16 | 2 Mb | 105 | 10/10 | 0 | detected |
| H1-H16 | 4 Mb | 77 | 10/10 | 0 | detected |
| 2026-06-15 | 0 bp | 201 | 0/0 | 0 | NA |
| 2026-06-15 | 1 Mb | 201 | 0/0 | 0 | NA |
| 2026-06-15 | 2 Mb | 165 | 0/0 | 0 | NA |
| 2026-06-15 | 4 Mb | 143 | 0/0 | 0 | NA |

Interpretation:

- `merge_gap_bp=2_000_000` is the first candidate setting worth formal
  materialization because it lowers burden while preserving all current
  Y1-Y8/H1-H6 truth events in the existing-BED ablation.
- `merge_gap_bp=4_000_000` lowers burden more, but it has greater risk of
  merging biologically separate same-direction events. It should remain
  exploratory until more truth and negative evidence exists.
- Merge-gap optimization changes candidate granularity only. It does not
  suppress WisecondorX signals or apply Branch B rules.

## Implementation

Added config plumbing:

```yaml
core:
  wisecondorx:
    cnv:
      branch_a:
        merge_gap_bp: 0
        strong_z: 10.0
```

Workflow behavior:

- `a_branch_candidate_assembly` now passes `--merge-gap-bp` and `--strong-z`
  to the existing `pgta.predict.branch_a` CLI.
- Defaults preserve the previous behavior.
- Branch A candidate summaries record the active `merge_gap_bp` and `strong_z`.

## Remote Validation

Remote unit tests:

```text
executable: /biosoftware/miniconda/envs/snakemake_env/bin/python
command: python -m pytest tests/unit/test_branch_a_candidates.py tests/unit/test_branch_a_validation.py tests/unit/test_branch_ab_phase12_workflow_contract.py tests/unit/test_workflow_line_endings_contract.py -q
result: 18 passed in 0.98s
```

Remote dry-runs:

```text
executable: /biosoftware/miniconda/envs/snakemake_env/bin/snakemake
command: snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile config_predict_y_h_r0_shadow_branch_a_validation_20260620.yaml --cores 1 -n branch_a_validation
result: DAG built; Branch A assembly and validation planned because code/params changed
```

```text
executable: /biosoftware/miniconda/envs/snakemake_env/bin/snakemake
command: snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile config_predict_y_h_r0_shadow_branch_b_evidence_20260620.yaml --cores 1 -n branch_b_evidence branch_s_review cnv_report
result: DAG built; downstream jobs planned from updated Branch A candidates
```

Remote materialization with default `merge_gap_bp=0` completed for:

```text
config_predict_y_h_r0_shadow_branch_a_validation_20260620.yaml
config_predict_h_20260608_h_r0_shadow_branch_a_validation_20260620.yaml
config_predict_20260615_h_r0_shadow_branch_a_validation_20260620.yaml
```

Result:

```text
Y1-Y8: FN=0, truth_detected=10/10, Branch A candidates=131
H1-H16: FN=0, truth_detected=10/10, H6 chr21=detected, Branch A candidates=221
2026-06-15: no truth table, Branch A candidates=201
```

## Current Decision

This phase does not promote a new Branch A default. It makes the already
implemented Branch A merge-gap behavior workflow-configurable and records
evidence that `merge_gap_bp=2_000_000` is the next formal materialization
candidate.

Before promoting `merge_gap_bp=2_000_000`, rerun the fixed A/B/S chain under the
same reference/config contract:

```text
WisecondorX predict output
-> Branch A candidates with merge_gap_bp=2_000_000
-> Branch A validation
-> Branch B V2-only evidence/disposition
-> Branch S review
-> evaluation / benchmark / report package
```

## Risks

- The current truth set is small. Existing-BED ablation can show that known
  truth events are retained, but it cannot prove general no-FN behavior.
- Merge gaps can merge distinct same-direction events separated by small
  intervals; this is acceptable for burden reduction only if downstream Branch B
  and report text still expose uncertainty and candidate provenance.
- 2026-06-15 has no truth labels and remains burden/context only.
- Z threshold tightening is not safe as a Branch A hard filter under the current
  truth set.
