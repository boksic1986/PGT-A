# 2026-06-21 02:59 Branch A Burden Phase 1 Handoff

## Scope

This handoff closes the first Branch A burden optimization loop after the remote
mirror / Snakemake parse repair. It does not promote a new production reference,
Branch A default, Branch B V2 decision rule, Branch S final rule, or P6 release
package.

## Context Entry Point

Read first:

1. `docs/CURRENT_CONTEXT_INDEX.md`
2. `docs/reports/branch_a_burden_optimization_phase1_2026-06-21.md`
3. `docs/handoff/2026-06-21_0127_p1_p6_credibility_audit_handoff.md`
4. `AGENTS.md`
5. `skills/conversation_handoff/SKILL.md`
6. `skills/pgta_reference_modeling_analysis/SKILL.md`

## Remote Mirror / Parse Status

Remote mirror:
`/data/project/CNV/PGT-A/refactor_validation_20260419`

The previous Snakemake parse blocker was repaired before this Branch A loop.
Remote workflow sources were LF-normalized and parse validation passed.

Remote checks:

```text
executable: /biosoftware/miniconda/envs/snakemake_env/bin/python
command: python -m pytest tests/unit/test_workflow_line_endings_contract.py tests/unit/test_current_context_index.py -q
result: 5 passed in 0.13s
```

```text
executable: /biosoftware/miniconda/envs/snakemake_env/bin/snakemake
command: snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile config_predict_y_h_r0_shadow_branch_b_evidence_20260620.yaml --cores 1 -n branch_b_evidence branch_s_review cnv_report
result: DAG built successfully
```

```text
executable: /biosoftware/miniconda/envs/snakemake_env/bin/snakemake
command: snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile config_predict_y_h_r0_shadow_branch_a_validation_20260620.yaml --cores 1 -n branch_a_validation
result: Nothing to be done before Branch A code changes; parse succeeded
```

## Code Changes

Branch:
`codex/branch-a-burden-merge-gap`

Changed tracked files:

- `pgta/core/init_project.py`
- `rules/predict_layout.smk`
- `rules/predict_workflow.smk`
- `tests/unit/test_branch_a_candidates.py`
- `tests/unit/test_branch_ab_phase12_workflow_contract.py`
- `tests/unit/test_current_context_index.py`
- `docs/CURRENT_CONTEXT_INDEX.md`
- `CURRENT_STATE.md`
- `PLANS.md`
- `docs/reports/branch_a_burden_optimization_phase1_2026-06-21.md`

Untracked temporary shell scripts, local report packages, images, Excel files,
and `validate_npz.py` are not part of this handoff and should not be committed
unless explicitly promoted into the formal workflow.

## Branch A Contract Change

Branch A remains WisecondorX predict/CBS-derived candidate discovery.

New workflow-configurable settings:

```yaml
core:
  wisecondorx:
    cnv:
      branch_a:
        merge_gap_bp: 0
        strong_z: 10.0
```

`a_branch_candidate_assembly` now passes `--merge-gap-bp` and `--strong-z` to
`pgta.predict.branch_a`.

Defaults preserve previous behavior:

```text
merge_gap_bp=0
strong_z=10.0
```

No Branch B artifact, matched-negative, calibration, or legacy disposition
field enters Branch A.

## Branch A Evidence

Remote materialization with default `merge_gap_bp=0` completed:

```text
Y1-Y8: truth_detected=10/10, FN=0, Branch A candidates=131
H1-H16: truth_detected=10/10, FN=0, H6 chr21=detected, Branch A candidates=221
2026-06-15: no truth table, Branch A candidates=201
```

Z-threshold ablation:

```text
z>=8 would miss H6 chr21 gain.
z>=10 would miss Y6 chr7 loss and H6 chr21 gain.
```

Therefore `strong_z=10` remains a strong/sensitive tier label only, not a hard
inclusion filter.

Existing-BED merge-gap ablation:

```text
merge_gap_bp=2_000_000 preserved all current Y1-Y8/H1-H6 truth events and reduced candidate burden.
merge_gap_bp=4_000_000 reduced more burden but is riskier for merging distinct same-direction events.
```

`merge_gap_bp=2_000_000` is the first formal candidate for the next run, but it
is not promoted as default in this handoff.

## Remote Validation

Synced changed code/tests to remote and LF-normalized them.

Remote unit tests:

```text
executable: /biosoftware/miniconda/envs/snakemake_env/bin/python
command: PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 python -m pytest tests/unit/test_branch_a_candidates.py tests/unit/test_branch_a_validation.py tests/unit/test_branch_ab_phase12_workflow_contract.py tests/unit/test_workflow_line_endings_contract.py -q
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

Remote materialization:

```text
executable: /biosoftware/miniconda/envs/snakemake_env/bin/snakemake
command: snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile config_predict_y_h_r0_shadow_branch_a_validation_20260620.yaml --cores 4 branch_a_validation
result: completed successfully
```

```text
executable: /biosoftware/miniconda/envs/snakemake_env/bin/snakemake
command: snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile config_predict_h_20260608_h_r0_shadow_branch_a_validation_20260620.yaml --cores 4 branch_a_validation
result: completed successfully
```

```text
executable: /biosoftware/miniconda/envs/snakemake_env/bin/snakemake
command: snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile config_predict_20260615_h_r0_shadow_branch_a_validation_20260620.yaml --cores 4 branch_a_validation
result: completed successfully
```

## Current Decision

Branch A Phase 1 is complete as an engineering and ablation loop.

It does not change:

- active reference ID;
- WisecondorX predict parameters;
- default Branch A inclusion behavior;
- Branch B V2 status;
- Branch S status;
- report-release status.

## Next Gate

Choose one before Branch B V2 truth benchmarking:

1. Keep default Branch A `merge_gap_bp=0` and proceed to V2-only Branch B
   evidence/disposition benchmark.
2. Materialize an explicit Branch A `merge_gap_bp=2_000_000` config and rerun:
   WisecondorX predict outputs already present -> Branch A candidates -> Branch
   A validation -> Branch B V2-only evidence/disposition -> Branch S review ->
   evaluation / benchmark / report package.

Do not hard-filter Branch A by `z>=8` or `z>=10` under the current evidence.
