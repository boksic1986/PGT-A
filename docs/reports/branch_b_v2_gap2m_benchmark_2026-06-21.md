# Branch B V2 Gap2m Benchmark

Date: 2026-06-21

Status: `materialized_shadow_benchmark`

## Scope

This report records the first materialized Branch B V2-only benchmark using the
explicit Branch A gap2m overlay:

```text
reference_id=h_r0_shadow_ref_20260619
branch_a.merge_gap_bp=2_000_000
branch_a.strong_z=10.0
branch_a.output_dir=wisecondorx/cnv/a_branch_gap2m
postprocess.output_dir=wisecondorx/cnv/postprocess_gap2m
```

This benchmark is not a final report gate and does not promote Branch B V2. It
only verifies that the V2 classifier/benchmark path can consume Branch A gap2m
candidates, ignore legacy Branch B decision fields, and preserve locked truth
overlap candidates without hard suppression.

## Workflow Contract

New workflow target:

```text
branch_b_v2_benchmark
```

The benchmark consumes:

- `CNV_B_V2_CLASSIFIER` rows.
- Locked truth table when configured.
- Current run metadata.

The benchmark does not consume:

- legacy `CNV_B_FINAL_EVENTS`;
- legacy artifact summary;
- matched-negative percentile outputs;
- report outputs.

Ignored legacy decision fields are recorded in summary JSON:

```text
final_disposition
branch_b_keep_event
branch_b_report_class
branch_b_artifact_status
```

`final_report_impact` remains:

```text
none_shadow_only
```

## Important Fix During Materialization

The first materialized benchmark exposed a statistics bug: the benchmark counted
`V2_REVIEW_NO_HARD_SUPPRESSION` as hard suppression because it used a broad
substring match on `SUPPRESS`.

The benchmark now counts hard suppression only for explicit hard-action or
artifact tokens, such as:

```text
V2_SUPPRESS_*
SUPPRESS_*
V2_HARD_SUPPRESS*
HARD_SUPPRESS*
*ARTIFACT*
```

`NO_HARD_SUPPRESSION` is therefore treated as review-preserved, not as a hard
suppression.

## Remote Result Files

Y1-Y8:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419/results_y_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m/v2_benchmark/summary.json
/data/project/CNV/PGT-A/refactor_validation_20260419/results_y_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m/v2_benchmark/truth_metrics.tsv
/data/project/CNV/PGT-A/refactor_validation_20260419/results_y_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m/v2_benchmark/sample_summary.tsv
```

H1-H16:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419/results_h_20260608_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m/v2_benchmark/summary.json
/data/project/CNV/PGT-A/refactor_validation_20260419/results_h_20260608_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m/v2_benchmark/truth_metrics.tsv
/data/project/CNV/PGT-A/refactor_validation_20260419/results_h_20260608_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m/v2_benchmark/sample_summary.tsv
```

2026-06-15 five-sample context:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419/results_20260615_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m/v2_benchmark/summary.json
/data/project/CNV/PGT-A/refactor_validation_20260419/results_20260615_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m/v2_benchmark/sample_summary.tsv
```

## Summary Metrics

| cohort | samples | candidates | truth events | truth preserved | FN | hard suppressed truth | positive-support candidates | status |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| Y1-Y8 | 8 | 97 | 10 | 10 | 0 | 0 | 88 | ready |
| H1-H16 | 16 | 105 | 10 | 10 | 0 | 0 | 67 | ready |
| 2026-06-15 | 5 | 165 | 0 | 0 | 0 | 0 | 138 | skipped_no_truth |

The 2026-06-15 cohort has no locked truth table. It is burden/context only and
cannot support TP/FN/FP claims.

## Truth Notes

Y1-Y8:

- All 10 locked truth events are preserved.
- Y3 chrX loss is carried as `V2_NO_CALL_CONTRACT_RISK` /
  `V2_REVIEW_NO_HARD_SUPPRESSION`, not as a positive-support call.

H1-H6 truth within H1-H16:

- All 10 locked truth events are preserved.
- H6 chr21 gain remains preserved with `top_a_abs_zscore=7.1135`.
- H5/H6 chrX loss truth events are carried as `V2_NO_CALL_CONTRACT_RISK` /
  `V2_REVIEW_NO_HARD_SUPPRESSION`, not as positive-support calls.

This means V2 currently preserves truth-overlap candidates but still needs a
separate SCA/Branch S truth gate before sex-chromosome/SCA calls can be final.

## Interpretation

What this benchmark supports:

- Branch B V2 benchmark path is workflow-materialized.
- Legacy/current-code Branch B decision fields are excluded from V2 benchmark
  decisions.
- Under the Branch A gap2m overlay, locked Y/H truth candidates are not hard
  suppressed by V2.
- H6 chr21 weak positive remains preserved.

What this benchmark does not support:

- It does not prove Branch B V2 reduces FP burden enough for final reporting.
- It does not make Branch S final.
- It does not make the 2026-06-15 samples report-release ready.
- It does not promote `merge_gap_bp=2_000_000` to default.
- It does not promote `h_r0_shadow_ref_20260619` to production reference.

## Remote Commands

Unit tests:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 /biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest tests/unit/test_branch_b_v2_classifier.py tests/unit/test_branch_b_v2_benchmark.py tests/unit/test_branch_ab_phase12_workflow_contract.py -q'
```

Result:

```text
27 passed in 0.77s
```

Materialization was run on `ssh fengxian` with:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile <branch_b_v2_gap2m_config> --cores 24 --rerun-incomplete branch_b_v2_benchmark
```

The corrected aggregate benchmark outputs were then forced with:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile <branch_b_v2_gap2m_config> --cores 4 --forcerun cnv_branch_b_v2_benchmark branch_b_v2_benchmark
```

All three configs completed:

- `config_predict_y_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`
- `config_predict_h_20260608_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`
- `config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`

Post-sync mirror validation:

- After syncing the committed files to the non-git remote mirror, Snakemake
  dry-run initially failed because `Snakefile` had CRLF line endings again.
- The remote workflow/Python/YAML sources touched by this loop were rewritten
  to LF using remote Python.
- After LF normalization, all three `branch_b_v2_benchmark` dry-runs parsed and
  returned `Nothing to be done`.

Dry-run command pattern:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile <branch_b_v2_gap2m_config> --cores 1 -n branch_b_v2_benchmark
```

## Next Gate

The next gate is Branch B V2 evidence/disposition refinement for FP/review
burden, still under shadow/no-final-impact mode. The immediate focus should be:

1. Keep truth-overlap preservation and no hard suppression as non-regression
   checks.
2. Separate autosomal CNV evidence from SCA/sex-chromosome review evidence.
3. Quantify positive-support burden on H7-H16 and 2026-06-15 context samples.
4. Decide which V2 evidence fields can support report review labels without
   hard suppression or final benign/artifact claims.
