# Branch S Materialization Check

Date: 2026-06-19

## Scope

This note records the first formal materialization of Branch S shadow outputs beyond the original Y1-Y8 run.

Branch S remains shadow-only. These outputs do not change:

- WisecondorX predict;
- Branch A candidate generation;
- legacy Branch B final events;
- current sex calling;
- `cnv_report`.

## Remote Environment

- remote: `fengxian`
- code path: `/data/project/CNV/PGT-A/refactor_validation_20260419`
- python: `/biosoftware/miniconda/envs/snakemake_env/bin/python`
- snakemake: `/biosoftware/miniconda/envs/snakemake_env/bin/snakemake`

## Commands And Results

2026-06-15 five-sample run:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && /biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile /data/project/CNV/PGT-A/refactor_validation_20260419/config_predict_20260615_mask_only.yaml --cores 8 cnv --rerun-triggers mtime'
```

Result:

```text
16 of 16 steps completed.
Complete log: .snakemake/log/2026-06-19T103219.027740.snakemake.log
```

H1-H16 run:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && /biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile /data/project/CNV/PGT-A/refactor_validation_20260419/config_predict_h_20260608_mask_only.yaml --cores 8 cnv --rerun-triggers mtime'
```

Result:

```text
49 of 49 steps completed.
Complete log: .snakemake/log/2026-06-19T103240.765088.snakemake.log
```

Final H dry-run:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && /biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile /data/project/CNV/PGT-A/refactor_validation_20260419/config_predict_h_20260608_mask_only.yaml --cores 1 -n cnv --rerun-triggers mtime'
```

Result:

```text
Nothing to be done (all requested files are present and up to date).
```

Final 2026-06-15 dry-run:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && /biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile /data/project/CNV/PGT-A/refactor_validation_20260419/config_predict_20260615_mask_only.yaml --cores 1 -n cnv --rerun-triggers mtime'
```

Result:

```text
Nothing to be done (all requested files are present and up to date).
```

## Materialized Output Summary

| run | Branch S dir exists | evidence tables | state-score tables | summary JSON | all shadow-only | replaces final report |
|---|---:|---:|---:|---:|---:|---:|
| Y1-Y8 | true | 8 | 8 | 8 | true | false |
| H1-H16 | true | 16 | 16 | 16 | true | false |
| 2026-06-15 | true | 5 | 5 | 5 | true | false |

Output roots:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419/results_build_ref_v2_mask_only/wisecondorx/cnv/postprocess/branch_s
/data/project/CNV/PGT-A/refactor_validation_20260419/results_h_20260608_mask_only/wisecondorx/cnv/postprocess/branch_s
/data/project/CNV/PGT-A/refactor_validation_20260419/results_20260615_mask_only/wisecondorx/cnv/postprocess/branch_s
```

## SCA Truth Check

Current chrX-loss truth samples show strong Branch A chrX loss candidates, but the current Branch S state-score sign is not biologically validated:

| run | sample | truth context | X non-PAR mean calibrated z | X non-PAR median calibrated z | X Branch A candidate count | current X scores |
|---|---|---|---:|---:|---:|---|
| Y1-Y8 | Y3 | chrX loss / XO-like | 2.4326 | -0.0036 | 2 | `X_GAIN=2.433`, `X_LOSS=-2.433` |
| H1-H16 | H5 | chrX loss segments, including mosaic | 2.1368 | -0.0044 | 2 | `X_GAIN=2.137`, `X_LOSS=-2.137` |
| H1-H16 | H6 | whole-chrX loss / XO-like | 1.9344 | -0.0050 | 2 | `X_GAIN=1.934`, `X_LOSS=-1.934` |

The formal Branch A candidates for these same samples are loss calls:

| sample | Branch A chrX candidates |
|---|---:|
| Y3 | 2 loss candidates, max abs z 77.06 |
| H5 | 2 loss candidates, max abs z 78.48 |
| H6 | 2 loss candidates, max abs z 73.65 |

Interpretation:

- Branch S is materialized and report-safe.
- Branch S region evidence is useful for SCA review design.
- The current `X_GAIN` / `X_LOSS` state scores must not be promoted to SCA classification because their sign does not currently align with known chrX-loss truth.
- Before any Branch S promotion, score direction must be reworked or redefined against a locked SCA truth set and sex-aware baseline.

## Current Conclusion

This completes materialization of Branch S shadow outputs for Y1-Y8, H1-H16, and the 2026-06-15 five-sample run.

The output is safe to keep in the workflow as shadow evidence because every summary remains `final_report_impact=none_shadow_only` and `replaces_final_report=false`.

It is not ready for clinical/reporting SCA decisions. The next technical task is to fix or explicitly redefine SCA state-score direction, then validate on additional locked SCA positives and clean XX/XY negatives.
