# Handoff: Branch B V2 Gap2m Benchmark

Date: 2026-06-21 09:37 Asia/Shanghai

Status: `active_current_handoff`

## 1. Goal

Continue the current PGT-A / CNV-seq R&D mainline by materializing the Branch B
V2-only benchmark on top of the explicit Branch A `merge_gap_bp=2_000_000`
overlay.

This loop focused on preservation/no-hard-suppression against locked Y/H truth.
It did not attempt final Branch B promotion, final Branch S/SCA promotion, or a
new P6/report package.

## 2. Context Restored

Read before execution:

- `docs/CURRENT_CONTEXT_INDEX.md`
- `docs/handoff/2026-06-21_0318_branch_a_gap2m_materialization_handoff.md`
- `AGENTS.md`
- `skills/conversation_handoff/SKILL.md`
- `skills/pgta_reference_modeling_analysis/SKILL.md`
- `CURRENT_STATE.md`
- `PLANS.md`
- active Branch B V2 benchmark configs and workflow files

Active reference:

```text
h_r0_shadow_ref_20260619
```

Active Branch A overlay:

```text
merge_gap_bp=2_000_000
strong_z=10.0
branch_a.output_dir=wisecondorx/cnv/a_branch_gap2m
postprocess.output_dir=wisecondorx/cnv/postprocess_gap2m
```

## 3. Completed Changes

Implemented Branch B V2 benchmark workflow contract:

- Added `pgta/predict/branch_b/v2_benchmark.py`.
- Added dispatcher entrypoint action `branch_b_v2_benchmark`.
- Added workflow target/rule `branch_b_v2_benchmark`.
- Added config-controlled postprocess output directory so gap2m V2 outputs do
  not overwrite default outputs.
- Added three gap2m V2 benchmark configs:
  - `config_predict_y_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`
  - `config_predict_h_20260608_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`
  - `config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`
- Updated Branch B V2 classifier summary to record ignored legacy decision
  fields.
- Added unit tests for classifier and benchmark contracts.

Important bug fixed during materialization:

- Initial benchmark counted `V2_REVIEW_NO_HARD_SUPPRESSION` as hard suppression
  because it used broad substring matching on `SUPPRESS`.
- The hard-suppression test now requires explicit hard-action or artifact
  tokens.
- `NO_HARD_SUPPRESSION` is now treated as review-preserved.

## 4. Remote Materialization

Remote mirror:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419
```

Background jobs launched first:

| cohort | PID | log | target |
|---|---:|---|---|
| Y1-Y8 | 98857 | `logs/branch_b_v2_gap2m_y_20260621.log` | `results_y_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m/v2_benchmark/summary.json` |
| H1-H16 | 98858 | `logs/branch_b_v2_gap2m_h_20260621.log` | `results_h_20260608_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m/v2_benchmark/summary.json` |
| 2026-06-15 | 98859 | `logs/branch_b_v2_gap2m_0615_20260621.log` | `results_20260615_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m/v2_benchmark/summary.json` |

`monitor/runtime.db` did not exist in the remote mirror before launch, so there
was no historical runtime estimate.

After the hard-suppression semantics fix, the aggregate benchmark rule was
forced for all three configs:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile <branch_b_v2_gap2m_config> --cores 4 --forcerun cnv_branch_b_v2_benchmark branch_b_v2_benchmark
```

Result:

```text
completed for Y1-Y8, H1-H16, and 2026-06-15
```

Post-commit remote mirror sync:

- The committed files were archived from the local worktree and extracted into
  `/data/project/CNV/PGT-A/refactor_validation_20260419`.
- The first post-sync dry-run failed because `Snakefile` was again present with
  CRLF line endings on the Linux mirror.
- Remote Python was used to rewrite the workflow/Python/YAML files touched by
  this loop to LF.
- After LF normalization, dry-runs for all three configs parsed successfully and
  returned `Nothing to be done`.

## 5. Result Summary

| cohort | samples | candidates | truth events | truth preserved | FN | hard-suppressed truth | positive-support candidates | status |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| Y1-Y8 | 8 | 97 | 10 | 10 | 0 | 0 | 88 | ready |
| H1-H16 | 16 | 105 | 10 | 10 | 0 | 0 | 67 | ready |
| 2026-06-15 | 5 | 165 | 0 | 0 | 0 | 0 | 138 | skipped_no_truth |

Key notes:

- H6 chr21 remains preserved with `top_a_abs_zscore=7.1135`.
- Y3/H5/H6 chrX truth events are preserved as review/no-hard-suppression, not
  as final SCA positive calls.
- 2026-06-15 has no locked truth table and remains burden/context only.

## 6. Verification

Remote unit tests:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 /biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest tests/unit/test_branch_b_v2_classifier.py tests/unit/test_branch_b_v2_benchmark.py tests/unit/test_branch_ab_phase12_workflow_contract.py -q'
```

Result:

```text
27 passed in 0.77s
```

Post-sync remote unit tests:

```text
27 passed in 0.91s
```

Post-sync dry-runs:

```text
Y1-Y8: Nothing to be done
H1-H16: Nothing to be done
2026-06-15: Nothing to be done
```

Local syntax check:

```text
python -m py_compile pgta\predict\branch_b\v2_benchmark.py
```

Result:

```text
passed
```

## 7. Current Conclusion

Branch B V2 gap2m benchmark now satisfies the preservation gate:

- locked Y/H truth preserved;
- FN=0 by V2 preservation;
- truth-overlap candidates are not hard-suppressed;
- legacy/current-code Branch B decision fields are excluded from V2 benchmark
  decisions.

This does not prove FP/review burden is solved. Branch B V2 remains
`none_shadow_only` and is not a final report decision layer.

## 8. Next Gate

Next work should refine Branch B V2 evidence/disposition for FP/review burden
while preserving the current no-FN/no-hard-suppression benchmark.

Recommended immediate checks:

1. Quantify positive-support burden on H7-H16 and 2026-06-15 context samples.
2. Separate autosomal CNV review evidence from Branch S/SCA review evidence.
3. Decide which V2 fields can become report review labels without hard
   suppression or final benign/artifact claims.
4. Keep `merge_gap_bp=2_000_000` as explicit overlay only; default remains
   `merge_gap_bp=0`.

## 9. Core File Sync

- `CURRENT_STATE.md`: updated with Branch B V2 gap2m benchmark results.
- `PLANS.md`: updated so the next Branch B work is FP/review burden refinement,
  not first benchmark implementation.
- `docs/CURRENT_CONTEXT_INDEX.md`: updated to make this benchmark the current
  Branch B context.
- `AGENTS.md`: not updated; hard constraints did not change.
- `REPO_MAP.md`: not updated; stable repo structure did not change.

## 10. Do Not Misread

- Do not treat this as Branch B V2 final promotion.
- Do not treat 2026-06-15 as TP/FN/FP evidence.
- Do not treat Branch S as final SCA.
- Do not use legacy/current-code Branch B kept/artifact/final-disposition fields
  as V2 decision evidence.
- Do not promote gap2m to default until downstream A/B/S/report gates are
  benchmarked and accepted.
