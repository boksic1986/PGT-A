---
status: active_handoff
decision_use: current_context
created_at: 2026-06-22 01:01 Asia/Shanghai
active_reference_id: h_r0_shadow_ref_20260619
branch_a_overlay: merge_gap_bp_2000000
final_impact: development_review_only
lowres_used: false
---

# Branch S / SCA V2 Sex-Aware Handoff

## 1. 本轮目标

修正 Branch S 的 SCA review 逻辑，避免正常 XY 样本中强 chrX Branch A
candidate 被单独解释成 `X_GAIN` / `SCA_REVIEW_STRONG`。

本轮只处理 Branch S / SCA：

- 不改 Branch A；
- 不改 reference；
- 不改 autosomal CNV report-layer filtering；
- 不读取、不引用、不等待 2Mb/3Mb lowres 结果；
- 不停止后台 lowres Snakemake 任务。

## 2. 已确认的项目状态

- Current context entrypoint:
  `docs/CURRENT_CONTEXT_INDEX.md`.
- Previous handoff:
  `docs/handoff/2026-06-21_2253_branch_b_v2_lowres_ref_evidence_handoff.md`.
- Active reference remains `h_r0_shadow_ref_20260619`, still
  `fixed_shadow_baseline_not_production`.
- Active Branch A input remains the explicit `merge_gap_bp=2_000_000`
  overlay; default Branch A behavior is unchanged.
- Branch B V2 report-layer filtering remains development-only.
- Branch S remains `review_reportable_with_limitations`, not final SCA.

Lowres background task was not touched. The PID file was present at:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419/logs/driver/lowres_ref_2mb_3mb_20260622.pid
```

with recorded PID:

```text
71117
```

## 3. 已完成事项

Code:

- `pgta/predict/branch_s.py`
  - added non-PAR median/robust direction corroboration;
  - Branch A sex-chrom candidates are anchors only;
  - uncorroborated XY chrX gain is demoted to no-call / sex-consistent audit;
  - strong `XX + X_LOSS` Branch A support is preserved as review to protect
    current Y3/H5/H6 truth evidence;
  - Branch S summaries now emit `sca_report_layer_class` and
    `sca_report_layer_reason`.
- `pgta/predict/report.py`
  - loads and reports Branch S report-layer fields;
  - adds `SCA_report_layer=...` to technical conclusion text.

Tests:

- `tests/unit/test_branch_s_shadow.py`
  - normal XY + strong chrX Branch A + neutral X non-PAR evidence must not
    produce strong SCA;
  - X-loss with non-PAR support remains report-review visible;
  - strong `XX + X_LOSS` Branch A support remains visible even when bin median
    is neutral.
- `tests/unit/test_cnv_report.py`
  - report loader and technical conclusion carry Branch S report-layer fields.

Report:

- `docs/reports/branch_s_sca_v2_sex_aware_review_2026-06-22.md`.

## 4. 当前结果

Remote materialization used active gap2m Branch S outputs:

```text
results_y_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m/branch_s
results_h_20260608_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m/branch_s
results_g1_g8_h_r0_shadow_ref_20260621/wisecondorx/cnv/postprocess_gap2m/branch_s
results_20260615_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m/branch_s
```

Cohort summary:

| cohort | samples | strong SCA review | weak SCA review | no-call | report review | internal review | sex-consistent/audit |
|---|---:|---:|---:|---:|---:|---:|---:|
| Y1-Y8 | 8 | 1 | 0 | 7 | 1 | 0 | 5 |
| H1-H16 | 16 | 2 | 1 | 13 | 2 | 1 | 12 |
| G1-G8 | 8 | 0 | 2 | 6 | 0 | 2 | 3 |
| 2026-06-15 | 5 | 0 | 2 | 3 | 0 | 2 | 3 |

Locked SCA truth visibility:

- Y3: `X_LOSS`, `SCA_REVIEW_STRONG`, `sca_report_review_event`.
- H5: `X_LOSS`, `SCA_REVIEW_STRONG`, `sca_report_review_event`.
- H6: `X_LOSS`, `SCA_REVIEW_STRONG`, `sca_report_review_event`.

H7-H16 normal XY context:

- all are `SCA_NO_CALL`;
- all chrX Branch A-only support is routed to
  `sca_filtered_or_sex_consistent_event`;
- no H7-H16 sample remains `X_GAIN` / `SCA_REVIEW_STRONG`.

## 5. 已验证内容

Remote unit tests:

```text
cd /data/project/CNV/PGT-A/refactor_validation_20260419
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_s_shadow.py \
  tests/unit/test_cnv_report.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py \
  tests/unit/test_current_context_index.py -q
```

Result:

```text
38 passed in 1.25s
```

Remote dry-run:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile <active gap2m config> --cores 1 -n branch_s_review cnv_report
```

Configs:

```text
config_predict_y_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml
config_predict_h_20260608_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml
config_predict_g1_g8_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml
config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml
```

Result:

```text
RC=0 for all four configs; all requested files were present and up to date.
```

## 6. 当前结论

本轮修正可以作为当前 Branch S/SCA development report-layer evidence：

- Branch A sex-chrom candidate 不再单独决定 strong SCA。
- 正常 XY 的 chrX Branch A-only 信号已从 strong SCA 降级为
  sex-consistent/audit 层。
- Y3/H5/H6 的 X-loss evidence 仍保持可见。

本轮不能作为：

- final SCA production promotion；
- SCA PASS/FAIL release；
- 2026-06-15 TP/FN/FP 结论；
- 2Mb/3Mb lowres evidence 结论。

## 7. 建议下一步

下一步应分两条线：

1. SCA truth gate：
   - 补 X gain、XXY、XYY、Y loss、mosaic SCA fraction series、跨批次 clean
     XX/XY negative；
   - 在 broader truth panel 上评估 `XX + X_LOSS` 保护逻辑是否需要进一步收紧。
2. Report integration:
   - 将 Branch S section 继续保持为 `review_reportable_with_limitations`；
   - 常染色体 CNV report 与 Branch S/SCA section 分离展示；
   - 2026-06-15 仍只作为 burden/context，直到有 locked truth。

## 8. 关键文件与路径

Changed code:

- `pgta/predict/branch_s.py`
- `pgta/predict/report.py`

Changed tests:

- `tests/unit/test_branch_s_shadow.py`
- `tests/unit/test_cnv_report.py`

New report:

- `docs/reports/branch_s_sca_v2_sex_aware_review_2026-06-22.md`

Remote result paths:

- `results_y_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m/branch_s`
- `results_h_20260608_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m/branch_s`
- `results_g1_g8_h_r0_shadow_ref_20260621/wisecondorx/cnv/postprocess_gap2m/branch_s`
- `results_20260615_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m/branch_s`

## 9. 环境约束

- All real validation must run on `ssh fengxian`.
- Remote mirror:
  `/data/project/CNV/PGT-A/refactor_validation_20260419`.
- Python:
  `/biosoftware/miniconda/envs/snakemake_env/bin/python`.
- Snakemake:
  `/biosoftware/miniconda/envs/snakemake_env/bin/snakemake`.
- Do not use local test results as pass evidence.
- Do not stage local temp shell scripts, Excel files, local `reports/`, images,
  or `validate_npz.py`.

## 10. Core file sync

- `docs/CURRENT_CONTEXT_INDEX.md`: updated to this handoff.
- `CURRENT_STATE.md`: updated with Branch S/SCA V2 status.
- `PLANS.md`: updated with the next Branch S truth/report gates.
- `REPO_MAP.md`: not updated; no stable directory-level ownership change.
- `AGENTS.md`: not updated; no hard-rule change.
