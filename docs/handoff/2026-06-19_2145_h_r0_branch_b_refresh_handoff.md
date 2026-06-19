# 2026-06-19 21:45 H R0 Branch B refresh handoff

## 本轮目标

在不重新生成 BAM 的前提下，继续 H R0 shadow reference rebuild 后的
WisecondorX/Branch A/Branch B 验证链路，修复 BAM 输入契约和 Branch B FN
风险，并把研发进度落到文档。

## 已读取上下文

- `docs/handoff/2026-06-19_2055_h_r0_post_rebuild_predict_handoff.md`
- `AGENTS.md`
- `skills/conversation_handoff/SKILL.md`
- `skills/pgta_reference_modeling_analysis/SKILL.md`
- `REPO_MAP.md`
- `PLANS.md`
- `CURRENT_STATE.md`
- `Snakefile`
- `rules/target_assembly.smk`
- `rules/predict_layout.smk`
- `rules/predict_workflow.smk`
- `pgta/predict/branch_b/artifact_rules.py`
- `tests/unit/test_branch_b_artifact_rules.py`
- `docs/reports/h_r0_shadow_post_rebuild_predict_2026-06-19.md`
- `docs/reports/branch_ab_v2_rnd_constraints_2026-06-18.md`

## 当前有效远端路径

Remote project:

`/data/project/CNV/PGT-A/refactor_validation_20260419`

H R0 shadow reference:

`/data/project/CNV/PGT-A/refactor_validation_20260419/results_h_r0_shadow_ref_20260619`

Post-rebuild predict outputs:

- Y1-Y8:
  `/data/project/CNV/PGT-A/refactor_validation_20260419/results_y_h_r0_shadow_ref_20260619`
- H1-H16:
  `/data/project/CNV/PGT-A/refactor_validation_20260419/results_h_20260608_h_r0_shadow_ref_20260619`
- 2026-06-15 five samples:
  `/data/project/CNV/PGT-A/refactor_validation_20260419/results_20260615_h_r0_shadow_ref_20260619`

## 已确认事实

1. H R0 shadow reference rebuild 已完成。
   - selected/inlier samples: `38`
   - H R0 samples included: `H9,H10,H11,H12,H15,H16`
   - best binsize: `750000`
   - project PCA diagnostic: `7`; 不得写成 WisecondorX `newref` CLI 参数。

2. BAM 不需要重做。
   - 当前 mapping contract 没有改变。
   - post-rebuild configs 使用现有 BAM/BAI：
     - Y: `results_build_ref_v2_mask_only/mapping`
     - H: `results_h_20260608_mask_only/mapping`
     - 0615: `results_20260615_mask_only/mapping`

3. 已修复 predict overlay 合并问题。
   - `Snakefile` 不再在 predict overlay 中丢弃 `samples`。
   - `build_reference` carryover 仍会从 predict overlay 中移除。

4. 已修复 predict/reference asset 边界。
   - `rules/target_assembly.smk` 不再让 predict target 生产
     reference annotation/mask assets。
   - predict 只消费已发布的 reference asset contract。

5. 已完成 Branch B artifact-rule 下游强制刷新。
   - Driver:
     `/data/project/CNV/PGT-A/refactor_validation_20260419/logs/driver/post_rebuild_branch_b_artifact_rules_refresh_20260619.sh`
   - Log:
     `/data/project/CNV/PGT-A/refactor_validation_20260419/logs/driver/post_rebuild_branch_b_artifact_rules_refresh_20260619.log`
   - Result:
     `DONE_ALL 2026-06-19T21:32:33+08:00`

## 远端命令记录

### 远端同步

Command:

```powershell
scp <changed files> fengxian:/data/project/CNV/PGT-A/refactor_validation_20260419/<path>
```

Normalization executable:

`/biosoftware/miniconda/envs/snakemake_env/bin/python`

Result:

`normalized 12 files`

### 结果解析

Executable:

`/biosoftware/miniconda/envs/snakemake_env/bin/python`

Command:

```bash
ssh fengxian "/biosoftware/miniconda/envs/snakemake_env/bin/python -" < parse script
```

Result:

- Y1-Y8: `truth_recall=1.0`, `A recall=1.0`, `B recall=1.0`
- H1-H16: `truth_recall=1.0`, `A recall=1.0`, `B recall=1.0`
- 0615 report outputs exist for all five samples.
- `negative_bank.matched_negative_ready=false`
- `matched_negative_blocking_reason=no_n0_locked_clean_negative_samples`

### 下游刷新日志 token check

Executable:

remote log read with `/biosoftware/miniconda/envs/snakemake_env/bin/python`

Result:

| token | count |
|---|---:|
| `rule map_reads` | 0 |
| `rule sort_bam` | 0 |
| `rule index_bam` | 0 |
| `wisecondorx_convert_for_cnv` | 0 |
| `wisecondorx_predict_cnv` | 0 |
| `cnv_artifact_rules_branch_b` | 32 |

Interpretation:

- This forced refresh did not regenerate BAM.
- This forced refresh did not rerun WisecondorX convert/predict.
- It only refreshed Branch B artifact rules and downstream report/evaluation
  products.

### Unit test

Executable:

`/biosoftware/miniconda/envs/snakemake_env/bin/python`

Command:

```bash
cd /data/project/CNV/PGT-A/refactor_validation_20260419
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
  /biosoftware/miniconda/envs/snakemake_env/bin/python \
  -m pytest tests/unit/test_branch_b_artifact_rules.py -q
```

Result:

`44 passed in 0.86s`

### Snakemake dry-run

Executable:

`/biosoftware/miniconda/envs/snakemake_env/bin/snakemake`

Command:

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

Result:

| config | dry-run | jobs | mapping/BAM rules | WisecondorX convert/predict | reference asset build |
|---|---|---:|---:|---:|---:|
| `config_predict_y_h_r0_shadow_20260619.yaml` | success | 45 | 0 | 0 | 0 |
| `config_predict_h_20260608_h_r0_shadow_20260619.yaml` | success | 85 | 0 | 0 | 0 |
| `config_predict_20260615_h_r0_shadow_20260619.yaml` | success | 28 | 0 | 0 | 0 |

## Branch B 修改摘要

Files:

- `pgta/predict/branch_b/artifact_rules.py`
- `rules/predict_layout.smk`
- `rules/predict_workflow.smk`
- `tests/unit/test_branch_b_artifact_rules.py`

Behavior:

- Strong/sensitive same-direction Branch A evidence can preserve an event as
  `review` when Branch B calibrated support is weak but not directionally
  contradictory.
- This is not a `PASS` promotion.
- This does not emit B-only final report events.

Key parameter changes:

- `a_branch_sensitive_review_max_high_risk_fraction`: `0.05 -> 0.20`
- New `a_branch_boundary_protect_min_abs_z`: `30.0`
- `narrow_boundary_artifact_protect_min_a_abs_z`: `50.0 -> 30.0`

Recovered FN-sensitive cases:

| sample | event | A support | refreshed Branch B result |
|---|---|---:|---|
| Y7 | `chr8:3000001-14250000 loss` | `a_z=-40.90` | `review`, `keep_event=1` |
| H6 | `chr21:15000001-42000000 gain` | `a_z=7.11` | `review`, `keep_event=1` |

## Post-refresh metrics

| set | events | kept/review | artifacts | truth detected | truth recall |
|---|---:|---:|---:|---:|---:|
| Y1-Y8 | 131 | 25 | 106 | 10/10 | 1.000 |
| H1-H16 | 221 | 20 | 201 | 10/10 | 1.000 |

H7-H16 post-rebuild report burden:

| sample | total B events | kept review |
|---|---:|---:|
| H7 | 14 | 0 |
| H8 | 15 | 0 |
| H9 | 12 | 0 |
| H10 | 12 | 0 |
| H11 | 12 | 0 |
| H12 | 13 | 0 |
| H13 | 23 | 1 |
| H14 | 22 | 0 |
| H15 | 12 | 0 |
| H16 | 12 | 0 |

2026-06-15 current-workflow report burden:

| sample | total B events | kept review | top event |
|---|---:|---:|---|
| JZ26125843-56-56 | 28 | 0 | A-branch signal only: `chr12 gain` |
| JZ26125844-59-59 | 38 | 3 | `chr4 gain` |
| JZ26125845-60-60 | 57 | 26 | `chr14 gain` |
| JZ26125846-61-61 | 44 | 3 | `chr12 gain` |
| JZ26125847-62-62 | 34 | 1 | `chr4 loss` |

## 当前仍未完成 / 不得误写为完成

1. Branch B V2 matched-negative calibration 仍未 ready。
   - `matched_negative_ready=false`
   - `N0=0`
   - 当前不能把 matched-negative percentile 用作 hard filter。

2. H7-H16 不能直接 promote 为 N0。
   - H R0 shadow ref 已纳入 H9/H10/H11/H12/H15/H16。
   - 这只证明可作为 shadow reference experiment 的一部分，不等于
     production reference 或 N0 locked negative。

3. 0615 五个样本只能作为 current-workflow exploratory report。
   - 不能作为 locked validation proof。

4. 当前 artifact-rule 修改是 recall-preserving review policy。
   - 它解决 FN 风险暴露问题。
   - 它没有完成 Branch B V2 final promotion。

## 已更新文档

- `docs/reports/h_r0_shadow_post_rebuild_predict_2026-06-19.md`
- `docs/reports/branch_ab_v2_rnd_constraints_2026-06-18.md`
- 本 handoff:
  `docs/handoff/2026-06-19_2145_h_r0_branch_b_refresh_handoff.md`

## 下一步建议

优先完成当前代码验证、提交和同步 main/worktree。提交时不要纳入本地历史
untracked helper scripts / sample Excel / 临时报告目录。

下一轮研发再处理：

- H7-H16 的 post-rebuild R0/R1/R2 与 N0/N1/N2 分离复核；
- 0615 五样本报告的人工 review；
- matched-negative calibration 的可用 N0 或替代 locked negative-background
  方案。
