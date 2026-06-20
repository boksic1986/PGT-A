# Handoff: Branch A Gap2m Materialization

Date: 2026-06-21 03:18 Asia/Shanghai

Status: `active_current_handoff`

## Context Source

Use `docs/CURRENT_CONTEXT_INDEX.md` as the first restoration entrypoint.

This handoff supersedes the prior Branch A burden Phase 1 handoff for current
execution state. Older handoffs remain historical context only when they
conflict with this file, current configs, or current remote results.

## Scope Completed

This loop handled the remote mirror / Snakemake parse validation issue and
materialized the first explicit Branch A burden-reduction candidate:

```text
merge_gap_bp=2_000_000
strong_z=10.0
reference_id=h_r0_shadow_ref_20260619
```

Branch A remains WisecondorX predict/CBS-derived candidate discovery. No
Branch B, Branch S, mosaic, sex-calling, or final-report logic was changed.

## Code / Config Changes

Workflow plumbing:

- `rules/predict_layout.smk`
  - Added config-controlled `CNV_BRANCH_A_OUTPUT_DIR`.
  - Added config-controlled `CNV_BRANCH_A_VALIDATION_DIR`.
  - Added config-controlled `CNV_BRANCH_A_LOG_DIR`.
  - Branch A candidate and validation outputs now use these variables.
- `rules/predict_workflow.smk`
  - Branch A assembly and validation logs now use `CNV_BRANCH_A_LOG_DIR`.
  - Reference candidate audit now reads candidates from
    `CNV_BRANCH_A_OUTPUT_DIR`.
- `pgta/core/init_project.py`
  - Added default Branch A output, validation, and log dir template fields.
- `tests/unit/test_branch_ab_phase12_workflow_contract.py`
  - Added contract coverage for isolated Branch A output/log dirs.

New overlay configs:

- `config_predict_y_h_r0_shadow_branch_a_gap2m_validation_20260621.yaml`
- `config_predict_h_20260608_h_r0_shadow_branch_a_gap2m_validation_20260621.yaml`
- `config_predict_20260615_h_r0_shadow_branch_a_gap2m_validation_20260621.yaml`

Overlay output contract:

```text
branch_a.output_dir=wisecondorx/cnv/a_branch_gap2m
branch_a.validation_dir=wisecondorx/cnv/branch_a_validation_gap2m
branch_a.log_dir=logs/cnv_branch_a_gap2m
```

This avoids overwriting the default P2 `merge_gap_bp=0` outputs.

## Remote Mirror / Parse Validation

Remote mirror:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419
```

Synced files:

```text
rules/predict_layout.smk
rules/predict_workflow.smk
pgta/core/init_project.py
tests/unit/test_branch_ab_phase12_workflow_contract.py
config_predict_y_h_r0_shadow_branch_a_gap2m_validation_20260621.yaml
config_predict_h_20260608_h_r0_shadow_branch_a_gap2m_validation_20260621.yaml
config_predict_20260615_h_r0_shadow_branch_a_gap2m_validation_20260621.yaml
```

Remote LF normalization was applied to those files and confirmed with
`crlf=False`.

Remote tests:

```text
executable: /biosoftware/miniconda/envs/snakemake_env/bin/python
command: PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 python -m pytest tests/unit/test_branch_a_candidates.py tests/unit/test_branch_a_validation.py tests/unit/test_branch_ab_phase12_workflow_contract.py tests/unit/test_workflow_line_endings_contract.py tests/unit/test_current_context_index.py -q
result: 22 passed in 1.05s
```

Remote dry-runs:

```text
executable: /biosoftware/miniconda/envs/snakemake_env/bin/snakemake
command: snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile <gap2m-config> --cores 1 -n branch_a_validation
result: parsed successfully for Y1-Y8, H1-H16, and 2026-06-15 overlays
```

Dry-runs planned only:

```text
a_branch_candidate_assembly
cnv_branch_a_validation
branch_a_validation
```

No WisecondorX predict jobs were requested by the gap2m overlays.

## Remote Materialization Results

Actual Snakemake materialization completed on `ssh fengxian` for:

```text
config_predict_y_h_r0_shadow_branch_a_gap2m_validation_20260621.yaml
config_predict_h_20260608_h_r0_shadow_branch_a_gap2m_validation_20260621.yaml
config_predict_20260615_h_r0_shadow_branch_a_gap2m_validation_20260621.yaml
```

Results:

| cohort | default candidates | gap2m candidates | truth events | detected | FN | note |
|---|---:|---:|---:|---:|---:|---|
| Y1-Y8 | 131 | 97 | 10 | 10 | 0 | no-FN retained |
| H1-H16 | 221 | 105 | 10 | 10 | 0 | H6 chr21 retained |
| 2026-06-15 | 201 | 165 | 0 | 0 | 0 | no truth; burden context only |

All summary JSON files record:

```text
branch_b_used_for_pass_fail=false
branch_s_used_for_pass_fail=false
report_outputs_used_for_pass_fail=false
reference_id=h_r0_shadow_ref_20260619
```

## Current Decision

`merge_gap_bp=2_000_000` is now a materialized Branch A candidate setting under
the active shadow reference. It reduced candidate burden while preserving the
current truth set.

It is not promoted to default yet.

Default Branch A remains:

```text
merge_gap_bp=0
strong_z=10.0
```

`strong_z=10` remains a strong/sensitive tier label only. It is not an inclusion
hard filter; weaker positives or mosaic candidates can remain below this value.

## Next Gate

Use the explicit gap2m overlay for the next Branch B V2-only benchmark unless a
blocking regression appears. Keep default `merge_gap_bp=0` as rollback/control.

Next implementation should:

1. Consume Branch A gap2m candidates.
2. Implement V2-only Branch B evidence/disposition benchmark.
3. Exclude legacy/current-code Branch B fields from decisions:
   `final_disposition`, `branch_b_keep_event`, legacy artifact labels,
   old `N0=0`, and N1-only matched-negative promotion.
4. Validate only against Y1-Y8/H1-H6 truth for TP/FN.
5. Keep H7-H16 and 2026-06-15 as burden/context unless locked truth labels are
   added.
6. Carry Branch S as review-reportable-with-limitations, not final SCA.
7. Generate a new P6/report package only after fixed A/B/S contracts are
   represented and limitations are visible.

## Do Not Do

- Do not treat old legacy/current-code Branch B result fields as V2 evidence.
- Do not treat the old 2026-06-15 P6 package as final report evidence.
- Do not promote `merge_gap_bp=2_000_000` to default before downstream A/B/S
  benchmark and report-gate validation.
- Do not hard-filter Branch A by `z>=10`.
- Do not rerun BAM or WisecondorX predict for the gap2m Branch A validation;
  this stage intentionally reuses existing predict BED outputs under the same
  reference/config contract.
