# 2026-06-19 17:10 H R0 shadow BAM reuse handoff

## 本轮目标

修复 H R0 shadow reference rebuild 误触发 `fastp_bwa` 的问题，并在不重跑 BAM 的前提下继续启动 reference rebuild。

## 已读取上下文

- Latest handoff:
  `docs/handoff/2026-06-19_1603_branch_ab_v2_plan_sync_handoff.md`
- `AGENTS.md`
- `skills/conversation_handoff/SKILL.md`
- `skills/pgta_reference_modeling_analysis/SKILL.md`
- `REPO_MAP.md`
- `PLANS.md`
- `CURRENT_STATE.md`
- `Snakefile`
- `rules/common_preprocess.smk`
- `rules/qc_workflow.smk`
- `rules/reference_workflow.smk`
- `rules/predict_workflow.smk`
- `rules/target_assembly.smk`
- `config_reference_h_r0_shadow_20260619.yaml`

## 根因

`SORTED_BAM` 原先固定为：

`core.project_path/mapping/{sample}.sorted.bam`

H R0 shadow reference 使用新的 `core.project_path` 后，该目录没有 BAM，
Snakemake 因而把 `fastp_bwa` 加入 DAG。由于 mapping reference 和 read
preprocessing contract 未改变，这一步不应重做 BAM。

## 已完成修改

新增统一输入契约：

- `samples.<sample_id>.bam`
- `samples.<sample_id>.bai`
- `samples.<sample_id>.fastp_json`

涉及文件：

- `Snakefile`
- `rules/qc_workflow.smk`
- `rules/reference_workflow.smk`
- `rules/predict_workflow.smk`
- `pgta/qc/report.py`
- `config_reference_h_r0_shadow_20260619.yaml`
- `docs/reports/branch_ab_v2_rnd_constraints_2026-06-18.md`
- `docs/reports/h_r0_shadow_rebuild_bam_input_2026-06-19.md`

## 远端同步

已同步到：

`/data/project/CNV/PGT-A/refactor_validation_20260419`

注意：第一次 scp 多文件时曾把部分文件展平成远端根目录副本，已删除：

- `/data/project/CNV/PGT-A/refactor_validation_20260419/qc_workflow.smk`
- `/data/project/CNV/PGT-A/refactor_validation_20260419/reference_workflow.smk`
- `/data/project/CNV/PGT-A/refactor_validation_20260419/predict_workflow.smk`
- `/data/project/CNV/PGT-A/refactor_validation_20260419/report.py`

远端文本文件已用 `perl -pi -e 's/\r$//'` 转为 LF，避免 Snakemake DSL
解析 CRLF 时出现假性 `IndentationError`。

## Remote validation

Dry-run command:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake \
  -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile \
  --configfile /data/project/CNV/PGT-A/refactor_validation_20260419/config_reference_h_r0_shadow_20260619.yaml \
  --directory /data/project/CNV/PGT-A/refactor_validation_20260419 \
  --cores 1 \
  -n \
  reference_qc reference
```

Result:

- success
- 18 jobs
- `fastp_bwa` absent from DAG
- `baseline_bam_uniformity_qc` inputs are existing BAMs from:
  - `results_build_ref_v2_mask_only/mapping`
  - `results_h_20260608_mask_only/mapping`

Artifact contract check:

- `ref_sample_count=42`
- `external_bam_count=42`
- `missing_bam_count=0`
- `missing_bai_count=0`

## Running task

PID:

`61980`

PID file:

`/data/project/CNV/PGT-A/refactor_validation_20260419/results_h_r0_shadow_ref_20260619/logs/driver/h_r0_shadow_reference_20260619.pid`

Log:

`/data/project/CNV/PGT-A/refactor_validation_20260419/results_h_r0_shadow_ref_20260619/logs/driver/h_r0_shadow_reference_20260619.snakemake.log`

Command:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake \
  -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile \
  --configfile /data/project/CNV/PGT-A/refactor_validation_20260419/config_reference_h_r0_shadow_20260619.yaml \
  --directory /data/project/CNV/PGT-A/refactor_validation_20260419 \
  --cores 60 \
  --rerun-incomplete \
  --latency-wait 60 \
  --printshellcmds \
  reference_qc reference
```

Current observed resource state shortly after launch:

- load average: `0.45, 1.59, 4.23`
- memory: about `10G` used of `1.4T`
- active job: `baseline_bam_uniformity_qc`
- no `fastp` or `bwa` job observed

Historical runtime evidence:

- Previous baseline BAM QC benchmark was about 7 minutes on a smaller cohort.
- Previous tuning benchmark was about 4 minutes.
- Previous sex-specific reference build benchmark was about 26 seconds.
- Current run has 42 samples and a larger binsize grid, so estimate is
  approximate; it should be much shorter than the original 9-11 hour mapping
  estimate.

## Next required steps

1. Check PID/log/targets for H R0 shadow reference completion.
2. If successful, record:
   - selected samples
   - inliers/outliers
   - best binsize
   - project QC PCA fields
   - reference output paths
3. Build predict config against the new reference.
4. Rerun:
   - WisecondorX predict
   - Branch A candidates
   - Branch B evidence/classification
   - negative-bank labels tied to the new reference ID
   - evaluation
   - benchmark
   - report
5. Do not reuse old `N0=0` as the post-rebuild Phase 3 result.

## Git state

Tracked/local files modified or added for this task should be committed after
the current verification checkpoint:

- `Snakefile`
- `rules/qc_workflow.smk`
- `rules/reference_workflow.smk`
- `rules/predict_workflow.smk`
- `pgta/qc/report.py`
- `docs/reports/branch_ab_v2_rnd_constraints_2026-06-18.md`
- `docs/reports/h_r0_shadow_rebuild_bam_input_2026-06-19.md`
- `docs/reports/h7_h16_reference_rebuild_eligibility_2026-06-19.tsv`
- `config_reference_h_r0_shadow_20260619.yaml`

`docs/handoff/` is ignored by default; force-add this handoff only if the
current commit policy requires handoff files in Git.
