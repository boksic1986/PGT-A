# H R0 Shadow Reference Post-Rebuild Predict Status

Date: 2026-06-19

## Scope

This report records the current H R0 shadow reference rebuild result and the
post-rebuild predict rerun setup. It is not a final validation report yet.

## Reference Rebuild Result

Remote project:

`/data/project/CNV/PGT-A/refactor_validation_20260419`

Reference result root:

`/data/project/CNV/PGT-A/refactor_validation_20260419/results_h_r0_shadow_ref_20260619`

Snakemake completion evidence:

- PID `61980` is no longer running.
- Driver log ended with `18 of 18 steps (100%) done`.
- Primary reference outputs exist under:
  - `reference/XX/result/ref_xx_best.npz`
  - `reference/XY/result/ref_xy_best.npz`
  - `reference/gender/result/ref_gender_best.npz`
  - `reference/gender/common_best_binsize.txt`

Selected reference samples:

- selected/inlier count: `38`
- H-batch R0 samples included: `H9,H10,H11,H12,H15,H16`
- H7, H8, H13, and H14 were not included in this first-pass R0 shadow ref.

Best tuning output:

```yaml
best_binsize: 750000
best_pca_components: 7
best_selected_pca: 7
best_elbow_pca: 7
best_cum_explained_variance: 0.497554
best_cv_reconstruction_mse: 0.8337433596
best_inliers: 38
best_outliers: 0
best_selection_method: maximize_inliers_then_minimize_outlier_fraction_then_minimize_cv_reconstruction_mse_then_smaller_binsize
signal_key: sample_dict_chr1_22_aligned
```

Interpretation:

- `best_binsize=750000` is the selected WisecondorX reference binsize for this
  shadow ref.
- The PCA fields are project-level tuning/QC diagnostics. They must not be
  described as WisecondorX `newref` CLI parameters unless the executable is
  separately verified to expose and use such a parameter.

## BAM Reuse Contract

The post-rebuild predict configs must reuse existing BAMs and must not request
`mapping`.

Implemented workflow change:

- `Snakefile` now allows predict overlays to merge `samples.<sample_id>.bam`
  and `samples.<sample_id>.bai` into the base sample config.
- `build_reference` carryover is still removed from predict overlays.

Post-rebuild predict config files:

- `config_predict_y_h_r0_shadow_20260619.yaml`
- `config_predict_h_20260608_h_r0_shadow_20260619.yaml`
- `config_predict_20260615_h_r0_shadow_20260619.yaml`

Existing BAM sources:

- Y1-Y8: `results_build_ref_v2_mask_only/mapping`
- H1-H16: `results_h_20260608_mask_only/mapping`
- 2026-06-15 five samples: `results_20260615_mask_only/mapping`

Remote BAM/BAI check:

- Y1-Y8: missing BAM `0`, missing BAI `0`
- H1-H16: missing BAM `0`, missing BAI `0`
- 2026-06-15 five samples: missing BAM `0`, missing BAI `0`

## Predict/Reference Asset Boundary Fix

Issue found during the first post-rebuild launch:

- The predict DAG did not include `mapping` or `fastp_bwa`.
- However, it did include `build_reference_bin_annotations` and
  `build_reference_masks`.
- That happened because the `cnv` target added reference asset files to
  `REFERENCE_ASSET_TARGET_FILES`, making predict a producer of reference assets.
- The unsafe driver was stopped: PID `108616`.

Impact:

- The stopped run removed incomplete reference asset bin files before it was
  killed.
- The reference assets were restored with the H R0 reference config by running
  only:
  - `build_reference_bin_annotations`
  - `build_reference_masks`
  - `export_reference_hard_mask_bed`
- No BAM or mapping rule was involved in asset restoration.

Fix:

- `rules/target_assembly.smk` now treats reference annotation/mask files as
  inputs consumed by predict, not as targets that predict is allowed to build.
- If these assets are missing, predict should fail fast rather than regenerate
  or overwrite a reference asset set.

Post-fix dry-run result:

- Y config: success; no `mapping`; no `fastp_bwa`; no reference asset build.
- H config: success; no `mapping`; no `fastp_bwa`; no reference asset build.
- 2026-06-15 config: success; no `mapping`; no `fastp_bwa`; no reference asset
  build.

## Completed Downstream Refresh

Remote executable:

`/biosoftware/miniconda/envs/snakemake_env/bin/snakemake`

Forced downstream-refresh driver:

`/data/project/CNV/PGT-A/refactor_validation_20260419/logs/driver/post_rebuild_branch_b_artifact_rules_refresh_20260619.sh`

Log:

`/data/project/CNV/PGT-A/refactor_validation_20260419/logs/driver/post_rebuild_branch_b_artifact_rules_refresh_20260619.log`

Completion evidence:

- Driver log contains `DONE_ALL 2026-06-19T21:32:33+08:00`.
- The final 2026-06-15 block reached `28 of 28 steps (100%) done`.
- The refresh was intentionally limited to downstream Branch B/report artifacts.

Log token check:

| token | count |
|---|---:|
| `rule map_reads` | 0 |
| `rule sort_bam` | 0 |
| `rule index_bam` | 0 |
| `wisecondorx_convert_for_cnv` | 0 |
| `wisecondorx_predict_cnv` | 0 |
| `cnv_artifact_rules_branch_b` | 32 |

Interpretation:

- No BAM generation was run.
- No WisecondorX convert or predict was rerun in this forced refresh.
- The refreshed outputs test only the Branch B artifact-rule change and
  downstream evaluation/report propagation.

Command shape:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake \
  -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile \
  --configfile <post-rebuild-config> \
  --directory /data/project/CNV/PGT-A/refactor_validation_20260419 \
  --cores 32 \
  --forcerun cnv_artifact_rules_branch_b \
  --rerun-incomplete \
  --rerun-triggers mtime \
  --printshellcmds
```

The refresh driver ran these configs sequentially:

1. `config_predict_y_h_r0_shadow_20260619.yaml`
2. `config_predict_h_20260608_h_r0_shadow_20260619.yaml`
3. `config_predict_20260615_h_r0_shadow_20260619.yaml`

## Branch B Artifact-Rule Adjustment

The post-rebuild Y/H truth scan exposed two FN-prone cases under the previous
hard artifact rules:

- `Y7 chr8 loss`: strong Branch A evidence (`a_z=-40.90`) was suppressed by
  boundary/blacklist/weak calibrated support flags.
- `H6 chr21 gain`: sensitive Branch A evidence (`a_z=7.11`) was suppressed by
  strict high-risk-fraction and broad-chromosome filters.

Implemented policy:

- Strong or sensitive Branch A evidence can preserve an event as `review` when
  Branch B calibrated evidence is weak but not contradicting direction.
- This does not convert the event to `PASS`.
- This does not create B-only final events.
- The affected event remains visible for review/reporting so FN is avoided while
  matched-negative evidence is still unavailable.

Key parameter changes:

| parameter | old | new | purpose |
|---|---:|---:|---|
| `a_branch_sensitive_review_max_high_risk_fraction` | `0.05` | `0.20` | allow moderate-risk-dominant small-chromosome positive candidates to remain review |
| `a_branch_boundary_protect_min_abs_z` | inherited `50` | `30` | separate boundary protection from discordant-protect threshold |
| `narrow_boundary_artifact_protect_min_a_abs_z` | `50` | `30` | preserve strong A-supported narrow boundary events as review |

## Post-Refresh Truth Metrics

| set | samples | events | kept/review | artifacts | truth detected | truth recall | A recall | B recall |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Y1-Y8 | 8 | 131 | 25 | 106 | 10/10 | 1.000 | 1.000 | 1.000 |
| H1-H16 | 16 | 221 | 20 | 201 | 10/10 | 1.000 | 1.000 | 1.000 |

Recovered FN-sensitive events:

| sample | truth event | A support | refreshed Branch B disposition | notes |
|---|---|---:|---|---|
| Y7 | `chr8:3000001-14250000 loss` | `a_z=-40.90` | `review`, `keep_event=1` | low calibrated signal, but strong A evidence preserved for review |
| H6 | `chr21:15000001-42000000 gain` | `a_z=7.11` | `review`, `keep_event=1` | sensitive A evidence preserved for review despite weak empirical support |

Y1-Y8 report burden after refresh:

| sample | total B events | kept review | top event |
|---|---:|---:|---|
| Y1 | 20 | 2 | `chr21:15000001-42000000 loss` |
| Y2 | 37 | 6 | `chr14:22500001-106500000 gain` |
| Y3 | 6 | 3 | `chrX:61500001-147750000 loss` |
| Y4 | 15 | 2 | `chr13:20250001-114000000 gain` |
| Y5 | 10 | 3 | `chr16:2250001-32250000 gain` |
| Y6 | 18 | 5 | `chr7:63000001-141750000 gain` |
| Y7 | 3 | 1 | `chr8:3000001-14250000 loss` |
| Y8 | 22 | 3 | `chr13:20250001-38250000 gain` |

H1-H16 report burden after refresh:

| sample | total B events | kept review | top event / conclusion |
|---|---:|---:|---|
| H1 | 26 | 2 | `chr16:1500001-31500000 gain` |
| H2 | 26 | 4 | `chr2:95250001-238500000 gain` |
| H3 | 3 | 1 | `chr13:20250001-114000000 gain` |
| H4 | 8 | 3 | `chr15:25500001-84000000 gain` |
| H5 | 5 | 4 | `chrX:62250001-147750000 loss` |
| H6 | 6 | 5 | `chr21:15000001-42000000 gain` |
| H7 | 14 | 0 | A-branch strong sensitive signal only |
| H8 | 15 | 0 | A-branch strong sensitive signal only |
| H9 | 12 | 0 | A-branch strong sensitive signal only |
| H10 | 12 | 0 | A-branch strong sensitive signal only |
| H11 | 12 | 0 | A-branch strong sensitive signal only |
| H12 | 13 | 0 | A-branch strong sensitive signal only |
| H13 | 23 | 1 | `chr21:15000001-42000000 gain` |
| H14 | 22 | 0 | A-branch strong sensitive signal only |
| H15 | 12 | 0 | A-branch strong sensitive signal only |
| H16 | 12 | 0 | A-branch strong sensitive signal only |

2026-06-15 five-sample current-workflow report burden:

| sample | total B events | kept review | top event / conclusion |
|---|---:|---:|---|
| JZ26125843-56-56 | 28 | 0 | A-branch strong sensitive signal only: `chr12:37500001-128250000 gain` |
| JZ26125844-59-59 | 38 | 3 | `chr4:52500001-186000000 gain` |
| JZ26125845-60-60 | 57 | 26 | `chr14:60750001-106500000 gain` |
| JZ26125846-61-61 | 44 | 3 | `chr12:37500001-124500000 gain` |
| JZ26125847-62-62 | 34 | 1 | `chr4:121500001-186000000 loss` |

## Remaining Constraints

Negative bank status is unchanged:

```json
{
  "version": "branch_ab_v2_h_r0_shadow_ref_20260619",
  "matched_negative_ready": false,
  "matched_negative_blocking_reason": "no_n0_locked_clean_negative_samples",
  "matched_negative_eligible_count": 0,
  "label_counts": {"N1": 6, "N2": 4}
}
```

Interpretation:

- The artifact-rule adjustment restores Y/H known-positive recall to 1.0 in
  this post-rebuild shadow reference context.
- It does not solve matched-negative calibration.
- H7-H16 should not be promoted to N0 based only on this report.
- The 2026-06-15 outputs can be used as current-workflow exploratory reports,
  not as locked validation proof.
- Branch B V2 final-report promotion remains blocked until post-rebuild N0 or
  another locked negative-background strategy is available and tested without
  recall regression.

## Remote Verification After Sync

Remote code mirror:

`/data/project/CNV/PGT-A/refactor_validation_20260419`

Unit test:

```bash
cd /data/project/CNV/PGT-A/refactor_validation_20260419
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
  /biosoftware/miniconda/envs/snakemake_env/bin/python \
  -m pytest tests/unit/test_branch_b_artifact_rules.py -q
```

Result:

`44 passed in 0.86s`

Post-rebuild dry-run check:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake \
  -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile \
  --configfile <post-rebuild-config> \
  --directory /data/project/CNV/PGT-A/refactor_validation_20260419 \
  --cores 1 \
  -n \
  --rerun-incomplete \
  --forcerun cnv_artifact_rules_branch_b \
  --rerun-triggers mtime
```

Dry-run results:

| config | jobs | mapping / BAM rules | WisecondorX convert/predict | reference asset build | expected downstream rule |
|---|---:|---:|---:|---:|---:|
| `config_predict_y_h_r0_shadow_20260619.yaml` | 45 | 0 | 0 | 0 | `cnv_artifact_rules_branch_b=11` |
| `config_predict_h_20260608_h_r0_shadow_20260619.yaml` | 85 | 0 | 0 | 0 | `cnv_artifact_rules_branch_b=19` |
| `config_predict_20260615_h_r0_shadow_20260619.yaml` | 28 | 0 | 0 | 0 | `cnv_artifact_rules_branch_b=8` |
