# PLANS.md

## 2026-06-22 0615 High-Confidence Review Next Gate

Current handoff:
`docs/handoff/2026-06-22_0907_0615_high_confidence_report_handoff.md`.

Current report:
`docs/reports/0615_high_confidence_report_candidates_2026-06-22.md`.

The current 0615 high-confidence review is complete as a read-only development
summary. It identifies 10 conservative autosomal review/report candidates, all
in `JZ26125845-60-60`, and excludes the 5/5 shared
`chr4:67.50-101.25Mb gain` region from high-confidence interpretation.

Immediate next steps:

1. Manually inspect the 10 `JZ26125845-60-60` high-confidence rows against the
   corresponding calibrated-z CNV plots.
2. Keep chr4 shared gains as batch/context review, not high-confidence report
   candidates.
3. Keep Branch S/SCA output as a separate development-only context section.
4. Do not use 0615 to derive new Branch B V2 thresholds or production filters.
5. Any future filter/demotion proposal must be ablated against Y1-Y8, H1-H6,
   and G1-G8 locked truth with FN=0, H6 chr21 retained, and G2 truth visible.

## 2026-06-22 Report Main And CNV Plot Next Gate

Current handoff:
`docs/handoff/2026-06-22_0437_report_main_cnv_plot_handoff.md`.

Current report:
`docs/reports/report_main_convergence_cnv_plot_2026-06-22.md`.

The report layer now has a usable development contract:

- final autosomal main-table events are `report_strong_event` and
  `report_weak_event`;
- internal review, filtered audit, and Branch S are separate from the autosomal
  main table;
- zero-report-event samples remain visible in sample summaries;
- every materialized sample has a calibrated-z WisecondorX-style CNV plot and
  plot-bin TSV.

The next gate is candidate-level interpretation of the remaining report burden,
not another broad threshold.

Immediate next steps:

1. Review `report_strong_event` and `report_weak_event` rows for Y/H/G/0615 by
   candidate-level evidence:
   - A z and length tier;
   - B-side signal context;
   - lowres same-direction support;
   - ref-MAD stability;
   - clean support;
   - acrocentric/telomere/centromere/high-repeat risk;
   - matching visual trend in the `.final_cnv.svg`.
2. Use Y1-Y8, H1-H6, and G1-G8 locked truth as the no-FN guardrail.
3. Keep H6 chr21 and G2 locked truth visible.
4. Treat 2026-06-15 as burden/context only until locked truth labels exist.
5. Do not use sample-level report event counts, truth labels, or 2026-06-15
   burden counts to derive candidate filters.
6. If a new demotion/filter candidate is proposed, run an ablation table first:
   affected truth rows, affected 0615 rows, and affected plot-visible events.

Current open issue:

- 2026-06-15 still has 71 autosomal report events across 5 samples. This is too
  high for a final-style report, but broad filtering risks FN. The next
  reduction must be candidate-level and truth-safe.

## 2026-06-22 Lowres Branch B/S Next Gate

Current handoff:
`docs/handoff/2026-06-22_0930_lowres_branch_bs_integration_handoff.md`.

Current report:
`docs/reports/branch_b_s_lowres_integration_2026-06-22.md`.

The 2Mb/3Mb shadow references are built, lowres predict is materialized, and the
Branch B/S workflow now has a lowres auxiliary evidence contract. The next gate
is interpretation of lowres/ref-MAD evidence, not new threshold design.

Immediate next steps:

1. Inspect lowres support labels for each locked truth event in Y/H/G.
2. Inspect lowres/ref-MAD labels for remaining `report_event` and
   `internal_review_event` rows.
3. Identify which events gain explanation from 2Mb/3Mb same-direction support.
4. Separately list events with no lowres support but stable ref-MAD; treat them
   as candidates for review, not automatic filtering.
5. Keep H6 chr21 and short/mosaic-sensitive events protected from lowres
   absence-based demotion.
6. Keep 2026-06-15 as burden/context only until locked truth labels exist.

Design constraints:

- Lowres same-direction support can improve confidence/explanation.
- Lowres absence cannot independently demote or filter.
- High ref-MAD weakens negative lowres interpretation.
- Branch S global non-PAR median is for whole-SCA trends; local segment
  non-PAR evidence protects small X/Y CNV review.
- PAR is secondary context only and cannot independently create SCA calls.

## 2026-06-22 Branch S / SCA Next Gate

Current handoff:
`docs/handoff/2026-06-22_0101_branch_s_sca_v2_handoff.md`.

Branch S/SCA V2 now has a sex-aware report-layer correction:

- normal XY chrX Branch A-only signals are routed to no-call / sex-consistent
  audit instead of strong SCA review;
- Y3/H5/H6 X-loss evidence remains visible;
- Branch S remains `review_reportable_with_limitations`, not final SCA.

Immediate next work should not change Branch A or reference. The next Branch S
gate is:

1. Build a broader locked SCA validation set:
   - X gain;
   - XXY;
   - XYY;
   - Y loss;
   - mosaic SCA fraction series;
   - PAR / XY-homology edge cases;
   - cross-batch clean XX/XY negatives.
2. Re-evaluate the conservative `XX + X_LOSS` Branch A preservation rule once
   the broader truth panel exists.
3. Keep Branch S as a separate report section:
   - `sca_report_review_event`;
   - `sca_internal_review_event`;
   - `sca_filtered_or_sex_consistent_event`;
   - `sca_no_call`.
4. Do not use 2026-06-15 for TP/FN/FP until locked truth labels exist.
5. Do not use 2Mb/3Mb low-resolution evidence for SCA decisions unless a
   separate lowres/SCA validation gate is explicitly opened.

## 2026-06-21 Whole R&D Context And Report Gate Override

Current context entrypoint:
`docs/CURRENT_CONTEXT_INDEX.md`.

This plan applies to the whole current R&D cycle, not only G1. The next work
must preserve the corrected interpretation:

- P2 Branch A no-FN under `h_r0_shadow_ref_20260619` is valid current
  sensitivity evidence.
- Branch A still needs burden optimization under FN=0 and H6 chr21 preservation
  constraints.
- Branch B V2 must be evaluated without legacy/current-code Branch B decision
  fields.
- Branch S is not final SCA, but it must be carried into report development as
  review-reportable-with-limitations rather than omitted.
- P6/report remains the delivery target after A/B/S contracts are strengthened;
  the old 0615 package is historical development-only evidence, not the final
  stopping point.

Immediate execution order:

1. Keep `docs/CURRENT_CONTEXT_INDEX.md` current and commit each completed
   logical loop.
2. Preserve LF line endings for workflow sources when syncing the remote mirror;
   the 2026-06-21 Snakemake parse blocker has been repaired and verified.
3. Branch A burden optimization Phase 1 is complete as config plumbing,
   ablation evidence, and isolated `merge_gap_bp=2_000_000` materialization.
   Defaults remain unchanged. The next benchmark can use the explicit gap2m
   overlay, but promotion still requires the fixed A/B/S chain and report gate.
4. Branch B V2-only benchmark has been materialized on the explicit gap2m
   overlay. It preserves Y1-Y8 and H1-H6 truth-overlap candidates without hard
   suppression, while still remaining `none_shadow_only`.
5. The first Branch B V2 burden refinement has been materialized: `chrX`/`chrY`
   candidates now route to `V2_SEX_CHROMOSOME_REVIEW` /
   `V2_ROUTE_BRANCH_S_REVIEW`, preserving evidence tiers and no-final-impact
   status while separating them from autosomal Branch B positive-support
   burden.
6. Next Branch B work should refine remaining autosomal FP/review burden without
   reintroducing legacy/current-code Branch B decision fields.
7. The remaining autosomal burden audit shows that Branch B-side direction
   support is useful as review evidence but unsafe as a hard filter or universal
   positive-support downgrade, because it would hit multiple locked truth top
   candidates.
8. The review-label-only direction-support contract has now been implemented
   and materialized. `v2_direction_support_label` and
   `v2_direction_support_reason` are review evidence only; they do not change
   candidate class, action, hard-suppression behavior, or final report impact.
9. The unresolved `UNKNOWN_BACKGROUND + NO_NULL_SUPPORT` burden is now
   explicit in Branch B V2 outputs through `v2_background_context_label`,
   `v2_background_context_reason`, and `background_context_label_counts`.
   These fields are review context only and must not be interpreted as benign,
   background-compatible, or final-reportable evidence.
10. Branch B V2 method reset is now materialized. The active V2 interpretation
   is Branch-A-anchored evidence stratification, not legacy/current-code Branch
   B filtering. The classifier now emits signal-strength, length, clean-support,
   GC/RC context, B-side signal context, and conservative V2 disposition fields.
11. Current first-round dispositions are `background_unknown_review` and
   `sca_branch_s_review` in the materialized cohorts. This preserves truth but
   does not yet reduce FP/review burden enough for final report promotion.
12. Branch B V2 now has an explicit truth-safe filter-action contract that has
   passed remote unit tests, dry-runs, and targeted materialization. The
   classifier emits `v2_filter_action`,
   `v2_filter_reason`, `v2_filter_scope`, and
   `v2_filter_hard_suppression_allowed`; the benchmark records filter actions
   in truth metrics and counts hard-suppressed candidates. Current materialized
   Y1-Y8/H1-H16 truth remains 10/10 with FN=0 and hard-suppressed truth=0.
13. The next Branch B V2 refinement should use the new disposition/context and
   filter-action fields to reduce burden cautiously. Do not convert background
   context, B-side signal context, length tier, GC/RC attenuation, or low clean
   support into a hard filter without locked-truth ablation. Current hard
   suppression is limited to workflow/reference contract risk.
14. Branch B V2 burden stratification is now materialized. It adds
   `v2_burden_reduction_*` fields and explicit CNVpro/CNVseq evidence tags, but
   it does not yet reduce total review burden: current outputs remain
   `background_unknown_review` and `branch_s_review`. Treat it as a reporting
   and auditability improvement, not FP-reduction proof.
15. Branch B V2 report-contract integration is now materialized:
   `cnv_report` displays V2 burden/status fields as
   `development_review_only` evidence without changing Branch A, adding hard
   filters, using legacy/current-code Branch B decision fields, or promoting
   Branch B V2 / Branch S.
16. Branch B V2 report-layer filtering is now materialized. It splits Branch A
   gap2m candidates into `report_event`, `internal_review_event`,
   `filtered_event`, and `branch_s_event`; `filtered_event` rows are audit-only.
   Current locked truth remains preserved with FN=0 and H6 chr21 retained.
17. G1-G8 current-scheme validation is now materialized as an additional
   positive cohort check under the same active reference and explicit gap2m
   overlay. It preserves G1-G8 truth 10/10 with FN=0 and no report-layer
   filtered truth, while exposing a high G2 review/report burden that should be
   reviewed before the next demotion pass.
18. The second Branch B V2 report-burden pass has been retracted as a current
   candidate-level decision rule. The condition `sample report_event count >=3`
   is now sample-level burden audit only and must not demote/filter individual
   candidates.
19. Corrected rematerialization is complete for Y/H/G/2026-06-15. Y/H/G truth
   remains preserved with FN=0; H6 chr21 and G2 truth remain visible; sample
   report burden is now audit-only in summary/report outputs.
20. The next burden-reduction design must begin with candidate-level evidence
    audit. Do not use truth labels, sample report counts, or 2026-06-15 burden
    counts to reverse-engineer demotion/filtering rules.
21. Any future demotion/suppression rule must be ablated against Y1-Y8/H1-H6
    and G1-G8, explicitly retain H6 chr21 and G1-G8 locked truth, and keep
    H7-H16/0615 as context unless locked truth labels are added.
22. The next candidate-level evidence audit can use the newly implemented
    low-resolution/ref-MAD interface as auxiliary context only. First build
    named 2Mb/3Mb shadow refs and rerun low-res WisecondorX predict for
    Y1-Y8, H1-H16, G1-G8, and 2026-06-15 as a separate long-task gate. Do not
    use low-resolution absence as a standalone filter; H6 chr21-like short weak
    positives must remain at least internal review.
23. Use ref-MAD stability to interpret whether a candidate region is unstable
    in the reference itself before treating low-res absence as meaningful.
    High-MAD regions should weaken negative low-res evidence, not strengthen
    filtering.
24. Upgrade Branch S toward `review_reportable_with_limitations`: visible
    SCA/sex-chromosome report section, controlled ref/negative-like FP burden,
    explicit uncertainty, and no final SCA promotion without locked truth.
25. Generate the next P6/report package from workflow outputs after fixed Branch
    A/B/S contracts are represented. The report should expose evidence levels and
    limitations rather than silently treating unknown evidence as benign.

## 2026-06-21 Current Gate Override

The immediate next gate is no longer P6 report review or Branch B V2 promotion.
After the P1-P6 credibility audit, the active next step is:

```text
G2/P2 evidence review -> V2-only Branch B truth benchmark design
```

Use:

- `docs/reports/p1_p6_result_credibility_audit_2026-06-21.md`
- `docs/reports/branch_ab_v2_gate_phase_plan_2026-06-20.md`
- `docs/reports/branch_ab_v2_rnd_constraints_2026-06-18.md`

Current planning constraints:

1. Do not use legacy/current-code Branch B kept/artifact/final-disposition
   fields as Branch B V2 performance evidence.
2. Do not use current Branch B V2 shadow classifier output as V2-only
   performance evidence; it still contains legacy-derived shadow classes.
3. Treat P3 as an evidence-ledger engineering artifact until a V2-only truth
   benchmark exists.
4. Treat Branch S as review/development-only until locked SCA truth validation
   passes.
5. Treat the 2026-06-15 report package as development-only because it has no
   truth table and cannot produce TP/FN/FP performance.
6. Before making new workflow-validation claims, fix or re-sync the remote
   mirror issue that currently causes Snakemake dry-run to fail at Snakefile
   parse time.

Minimum next implementation plan:

- Keep the explicit Branch A `merge_gap_bp=2_000_000` overlay as the current
  Branch B V2 benchmark input while keeping default `merge_gap_bp=0` as
  rollback/control.
- Treat the materialized V2-only benchmark as a preservation gate, not a final
  performance gate: Y1-Y8 and H1-H6 truth are preserved 10/10 with FN=0 and no
  hard-suppressed truth, but FP/review burden is not yet solved.
- Treat the sex-route refinement as a routing/burden clarification gate, not a
  final SCA gate: sex-chromosome candidates now route to Branch S review, while
  autosomal Branch B positive-support burden remains the next refinement target.
- Do not implement a Branch B direction-support hard filter. Direction support
  is now present in workflow outputs as a review label only and preserves
  truth-overlap candidates and H6 chr21 under the materialized gap2m benchmark.
- Do not treat `UNKNOWN_BACKGROUND_NO_NULL_SUPPORT` as benign. It now appears
  as explicit V2 review context and marks missing matched-negative plus missing
  calibration-null support.
- Do not treat `B_SIGNAL_DISCORDANT_WITH_A_DIRECTION` as Branch B calling an
  opposite event. It is a B-side auxiliary signal discordance review label only.
- Use the new V2 disposition fields as review/report-support context:
  `background_unknown_review`, `sca_branch_s_review`, `technical_risk_review`,
  `review_candidate`, and `report_candidate`. In current materialized cohorts,
  only `background_unknown_review` and `sca_branch_s_review` are present.
- Continue excluding `final_disposition`, `branch_b_keep_event`, legacy artifact
  status, and legacy kept counts from V2 decision and metric computation.
- Refine Branch B V2 evidence/disposition so it can reduce FP/review burden
  without hard suppression of truth-overlap candidates. Any future reduction
  rule needs a locked Y1-Y8/H1-H6 truth-safe ablation and explicit H6 chr21
  preservation.
- Keep H7-H16 and 0615 as burden/context only unless locked truth labels are
  added.

## 2026-06-20 Current Gate Override

The current execution order is defined by
`docs/reports/branch_ab_v2_gate_phase_plan_2026-06-20.md`.

Older Branch B V2 / N1 shadow-report sections below are retained as historical
implementation notes only. They must not be used as the active plan when they
conflict with the 2026-06-20 gate plan or
`docs/reports/branch_ab_v2_rnd_constraints_2026-06-18.md`.

Current active sequence:

1. G1/P1 reference cohort audit and R0/R1/R2 decision. Completed for
   `h_r0_shadow_ref_20260619`.
2. P2 Branch A/reference contract freeze and no-FN validation. Materialized for
   Y1-Y8, H1-H16, and 2026-06-15 exploratory burden review; current truth
   FN=0 for Y1-Y8 and H1-H6.
3. P3 Branch B candidate-level evidence refinement. Implemented and
   materialized as review-safe evidence ledger under `h_r0_shadow_ref_20260619`;
   current output remains shadow/review-only and is not report promotion.
4. P5 Branch S/SCA report-boundary validation. Implemented and materialized as
   review/development-only summary schema under `h_r0_shadow_ref_20260619`;
   G5-review boundary is represented, but G5-final SCA promotion is not passed.
5. P6 report package development/report-boundary integration is implemented for
   the 2026-06-15 five-sample package. The package now visibly carries P3/P5
   review/development status and does not claim final SCA or Branch B
   promotion. Release remains blocked until an explicit G6 release decision and
   upstream promotion gates pass.

Current next step:

- Review the generated P6 report package as an internal exploratory artifact and
  decide whether to continue toward G6 release validation or return to upstream
  reference/Branch B/Branch S promotion gates. Do not treat the package as a
  final clinical/report release.

## 文档定位
本文件只描述当前周期计划、优先级和执行顺序，不承担长期事实快照职责。

## 当前总目标
在保留 A branch 的前提下，继续完成 Branch B 的工程化收敛，并把 FAIL 样本 re-sequencing 后重新进入 reference 建模的路径正式纳入 workflow 设计，使其在 remote 环境下同时满足：
1. 保持独立于 A branch 的推断链路。
2. 维持当前已获得的 `truth_recall = 1.0`。
3. 控制 review 泛化与 report/top-event 副作用。
4. 清掉与 detector 结论无关的尾部流程故障。
5. 让重测样本回流通过明确 manifest、QC、cohort 决策和 reference build 验证完成，而不是通过手工补丁加入。

## 当前阶段判断
当前周期仍然处于：

`Branch B 收敛 + re-sequencing cohort 回流设计阶段`

补充说明：
- 本轮做的代码结构审查、runtime tracking 修复、report 输入闭环修复和 modular dry-run，都属于中间代码质控与工程收口。
- 这些工作服务于流程稳定性与可验证性，但不改变项目算法主线，也不替代 Branch B 收敛计划。

## 当前周期目标
### 目标 1：稳住当前 Branch B recall
- 保持 remote `truth_recall = 1.0`
- 保持 remote `branch_b_detection_rate = 1.0`
- 不把 A branch 信息重新引入 B branch 主推断链路

### 目标 2：维持 report/top-event 副作用受控
- 保持 `sex-aware ranking / report suppression` 当前效果
- 当前 remote `cnv_summary.tsv` 中继续不出现 `chrX / chrY` review 事件占据 `branch_b_top_event`
- 不回退 sex chromosome truth recall

### 目标 3：维持 runtime tracking 与流程尾部闭环稳定
- `collect_runtime_tracking` 已在 remote targeted rerun 中恢复成功
- 不回退当前 benchmark / runtime tracking 兼容修复
- 后续若继续改 rule 或 report，不能重新破坏 targeted rerun 的 clean success

### 目标 4：评估是否还需要进一步收敛 review
- 当前新增恢复的 `Y6 chr7 loss` 仍是高风险 kept review
- 需要确认是否还要继续收敛这类高风险 kept review
- 任何收敛动作都不能回退 `1.0` recall
- `segmental_duplication_overlap` 类 kept review 已完成首轮收缩，不再是当前最高优先级

### 目标 5：建立 re-sequencing 样本回流路径
- FAIL 样本需要 re-sequencing 后重新评估，不能被静默视为永久排除
- 重测样本必须先进入显式 manifest / config / QC 流程，再进入 reference cohort 候选
- 原始失败样本与重测样本之间的替代或比较策略必须可追踪
- 重新纳入 reference 建模后，需要用 remote 真实运行或 targeted validation 证明没有破坏 Branch B 结论

### 目标 6：收敛补丁式实现
- 将明确重复的 Branch B rule 定义收敛为单一实现
- 将 predict evaluation / benchmark / report 的 truth/event 公共逻辑收敛到共享模块
- 对兼容入口和固定参数 rule 写清保留理由与退出条件
- 不新增 `_fix` / `_tmp` / `_new` / `_final` 类型平行脚本或 workflow 分叉

## 当前执行顺序
1. 先继续围绕 Branch B 结果做收敛判断。
   - 优先看当前 kept review 的剩余来源是否还需要收缩
   - 所有判断都以 `evaluation / benchmark / report` 三套结果回核
2. 补齐 re-sequencing 样本回流设计。
   - 先建立 manifest / cohort 决策 / reference build 触发条件
   - 再决定是否修改 `init_project.py`、`reference/cohort.py` 或配置 schema
3. 并行或穿插完成中间代码质控。
   - 包括 runtime tracking、report 输入闭环、remote dry-run 这类工程稳定性工作
   - 但这些工作不应升级成新的算法主线
4. 收敛已确认重复实现。
   - 优先处理 predict truth/event 公共逻辑和 Branch B rule 重复
   - 暂不动核心 CNV calling、mosaic、sex calling 和结果 schema
5. 然后再决定是否要处理与当前 workflow 主线无关的陈旧测试或兼容层问题。
6. 最后再判断当前 cycle 是否可以冻结。

## 当前里程碑
### 已完成
1. Branch B detector 已完成最小可运行链路。
2. WisecondorX NPZ sex chromosome numeric key 已接入 Branch B。
3. `sex-aware ranking / report suppression` 已在 remote report rerun 中生效。
4. `Y6 chr7 loss` 已通过 `masked_gap_rescue_detector` 恢复进入 Branch B kept review。
5. `benchmark.py` 与 `evaluation.py` 已统一 truth overlap 边界语义。
6. remote `fengxian` 复测下，Branch B 已达到 `10/10` truth detected。
7. `collect_runtime_tracking` 已恢复，targeted workflow 已 clean success。
8. `segmental_duplication_overlap` 类 kept review 已从 `16` 降到 `0`，且 remote recall 仍保持 `1.0`。
9. 本轮代码结构审查发现的两个工程问题已修复并推送：
   - `e784187` `fix: align runtime tracking and report inputs`
10. 上述修复已通过 remote dry-run 复核。
11. 已确认需要把 FAIL 样本 re-sequencing 回流作为当前周期正式任务线，而不是后续手工补丁。
12. 已完成 Branch B 后处理 rule 重复定义收敛，四个 Branch B 后处理 rule 现在各只有单一实现。
13. 已完成 `pgta/predict/evaluation.py` 与 `pgta/predict/benchmark.py` 的 truth/event 公共逻辑收敛，公共实现位于 `pgta/predict/truth.py`。
14. 已完成 re-sequencing manifest 的首轮正式接入：
   - `build_reference.resequencing.manifest`
   - `build_reference.resequencing.allowed_statuses_for_reference`
   - 默认只有 `promoted` 状态可进入 reference cohort selection
   - `replace_original` 会在 cohort selection 中替换原样本
15. 已完成 `build_reference.groups` 缺样本 fail-fast 校验。
16. 已完成 Branch A/B V2 Phase 1 shadow evidence ledger：
   - 每个 Branch A candidate 一行；
   - 不改变 WisecondorX predict、Branch A inclusion、Branch B final events 或 `cnv_report`。
17. 已完成 Branch A/B V2 Phase 2 N0/N1/N2 negative-bank labeling：
   - 当前 H7-H16 seed 中 `N0=0`、`N1=6`、`N2=4`；
   - 只有 N0 可用于 matched-negative empirical null。
18. 已完成 Branch A/B V2 Phase 3 matched-negative percentile shadow path：
   - 输出 `wisecondorx/cnv/postprocess/matched_negative/{sample}.candidate_evidence.tsv`；
   - 当前因为没有 N0，Y1-Y8 全部输出 `UNKNOWN_BACKGROUND` / `REVIEW_NO_CALL`；
   - 不进入 `cnv_report`，不作为 hard artifact filter。
19. 已完成 Branch S / SCA 验证计划文档：
   - `docs/reports/branch_s_sca_validation_plan_2026-06-19.md`
   - 明确 Branch S 只做 shadow review；
   - 明确当前 SCA truth set 不足以 promotion；
   - 明确 2026-06-15 五个样本只能按 current-workflow report 输出，不能作为 locked validation proof。
20. 已完成 Branch S shadow 输出物化：
   - Y1-Y8、H1-H16、2026-06-15 五个样本均已有 Branch S evidence / state-score / summary 输出；
   - 所有 summary 仍为 `final_report_impact=none_shadow_only`、`replaces_final_report=false`；
   - 物化报告：`docs/reports/branch_s_materialization_2026-06-19.md`。
21. 已完成 negative-bank readiness contract hardening：
   - `negative_bank_summary.json` 现在明确记录 `matched_negative_ready` 和阻塞原因；
   - 当前远端结果为 `matched_negative_ready=false`、`matched_negative_blocking_reason=no_n0_locked_clean_negative_samples`；
   - 详细报告：`docs/reports/branch_ab_v2_negative_bank_readiness_2026-06-19.md`。

### 待完成
1. review 收敛是否继续推进的决策与验证。
2. 基于当前 `main` 的 remote 非 dry-run 实跑验证。
3. 是否要进一步清理与当前 repo 布局不一致的陈旧测试断言。
4. re-sequencing 真实样本到位后，执行 `mapping / baseline_qc / reference_qc / reference` 非 dry-run 验证。
5. 记录 selected cohort metadata 到 reference 发布产物旁边。
6. 收敛 predict truth/event/io 公共逻辑中的剩余 IO 重复。
7. 建立 `R0/R1/R2` reference-rebuild eligibility 表，先把 H7-H16 作为新批次阴性候选重新分层；不要用旧 32-ref 下的 Branch A 信号作为唯一排除标准。
8. 构建至少一个 named shadow reference variant，并重跑 WisecondorX predict、Branch A candidate、Branch B、evaluation、benchmark 和 report。
9. 若 H7-H16 或其子集进入 named shadow reference，必须在新 ref 后重新运行 WisecondorX predict、Branch A、Branch B evidence、negative-bank labeling、evaluation、benchmark 和 report，并生成新的 negative-bank version；旧 32-ref 下的 `N0=0` 只能作为历史证据，不能作为新 ref 后的 Phase 3 结论。
10. 定义真正可用的 post-rebuild N0 locked clean negative cohort 后，再重新运行 Phase 3 以产生可解释 percentile；在此之前不得用 N1/H7-H16 构造 matched-negative empirical null。
11. Branch B V2 final-report promotion 需要单独门槛：post-rebuild FN=0、H6 chr21 仍被跟踪、FP/review burden 不增加或下降、unknown background 不做 hard suppression、所有 disposition rule 有 rollback path。
12. Branch S 仍需 locked SCA truth validation，不能因为已物化或 autosomal/reference rebuild 完成就 promotion。
13. 补充 SCA truth set：X gain、XXY、XXX、XYY、Y loss、mosaic fraction series、PAR/XY-homology edge cases 和跨批次 clean XX/XY negative。

## 当前不优先事项
1. 不直接 promote 扩 reference 样本；允许先做 named shadow reference rebuild。
2. 不优先做仓库级重构。
3. 不优先引入新的大规模模型或额外复杂 detector。
4. 不优先替换 A branch。
5. 不把中间代码质控升级成项目算法主线。
6. 不在 re-sequencing manifest 设计前直接把重测样本手工写进 reference cohort。

## 进入下一阶段的条件
只有同时满足以下条件，才考虑把当前 cycle 视为基本完成：
1. Branch B remote recall 稳定维持在 `1.0`。
2. report/top-event 的 sex-chromosome 副作用继续保持受控。
3. `collect_runtime_tracking` 持续不再阻塞 targeted rerun 的 clean success。
4. 当前需要保留的 review 是否继续收敛，已经形成明确结论。
5. 当前 `main` 上的模块化代码至少完成一轮 remote 非 dry-run 真实验证。
6. re-sequencing 样本回流路径已经形成明确 manifest / QC / cohort / reference build 契约。
7. 已确认的重复实现已收敛，或已有明确暂缓理由和后续验证计划。

## 2026-06-18 Branch A/B V2 研发计划更新

当前研发计划以以下文档为准：

- `docs/reports/branch_ab_v2_rnd_constraints_2026-06-18.md`
- `docs/reports/branch_b_refmap_calibration_risk_2026-06-17.md`

本轮方向调整：

1. 冻结 legacy Branch B，停止继续叠加一次性 artifact patch；legacy Branch B 只作为 ablation 对照。
2. 保持 WisecondorX-first：
   - Branch A 负责高灵敏度候选发现；
   - Branch B V2 只做候选级 evidence / classification shadow mode；
   - Branch B V2 默认不得输出 B-only final report event。
3. 第一阶段只实现 shadow evidence ledger：
   - 每个 Branch A candidate 一行；
   - 记录 raw/corrected amplitude、attenuation、flank contrast、region-risk、refmap/calibration-null 状态；
   - 记录最终 disposition；
   - 不改变当前 final report。
4. 在使用阴性样本做 matched-negative 统计前，先建立 N0/N1/N2 分层。
5. Branch S 作为后续 SCA shadow model 单独推进，不并入 Branch B V2 第一阶段。

H7-H16 reference 决策：

- H7-H16 当前均为 QC PASS/XY，是新批次阴性样本和 reference-rebuild 候选。
- 旧 32-ref 下的 Branch A signal 可能代表 reference mismatch、batch shift、高重复区域不稳定或真实个体异常；不能单独作为 reference rebuild 排除理由。
- reference rebuild 使用 `R0/R1/R2` 标签；Branch B matched-negative calibration 使用 `N0/N1/N2` 标签，两者不能混用。
- H8、H13、H14 在 Branch B calibration 中仍不应作为 `N0`，但可作为 `R1` shadow-rebuild 候选进入 ablation。
- H9、H10、H11、H12、H15、H16 可优先作为 `R0/R1` shadow reference evaluation 候选。
- 不得把任何 H7-H16 样本直接 promote 到 production reference 或 N0 high-confidence negative bank；必须先完成 named shadow ref rebuild 和重跑验证。新 ref 后需要重新计算 `N0/N1/N2`，不能沿用旧 32-ref 的 `N0=0` 作为最终 Phase 3 状态。

Reference rebuild / Branch A NPZ 契约：

- 若 binsize、mask、preprocess contract 不变，已有 predict NPZ 可以复用。
- 构建新 reference 后必须重跑 WisecondorX predict、Branch A candidate、Branch B、evaluation、benchmark 和 report。
- 若 binsize、mask 或 preprocessing 改变，reference 和 predict 样本的 raw/mask-only NPZ 都必须重建。
- 不同 reference ID 下生成的 candidate/evidence 不得直接混比，必须记录 reference ID、binsize、mask/preprocess version 和 candidate source。
- 新 reference 后必须单独评估 Branch A 性能：Y1-Y8/H1-H6 recall、H6 chr21 状态、candidate burden、建模阴性样本是否仍有样本特异强阳性信号，以及 residual signal 是否落在 CNVpro 明确特殊处理的区域。当前已确认 CNVpro 对 chr13/14/15/21/22 采用只分段 qter 的 acrocentric segmentation contract；这不是猜测性的 chr13/14/15 短臂 hard filter。

Branch A 迭代约束必须以现有代码为准：

- Branch A 当前输入是 WisecondorX predict 生成的 `{sample}_aberrations.bed`。
- WisecondorX predict 当前默认参数来自 workflow/config：`zscore=5`、`alpha=0.001`、`maskrepeats=5`、`minrefbins=150`、`seed=1`。
- `pgta/predict/branch_a.py` 只做候选表组装：标准化字段、剔除非法行、同样本/同染色体/同方向邻接或重叠合并、按 `a_abs_zscore` 排名。
- `strong_z=10.0` 当前只是 `strong/sensitive` 标签阈值，不是 candidate inclusion threshold。
- Branch A 不应加入 artifact filter、region-risk filter 或 Branch B calibrated evidence filter。
- 若要调整 WisecondorX predict 参数或 Branch A assembly 行为，必须先写 RFC 并在 remote 对 Y1-Y8、H1-H6、H6 chr21 和 H7-H16 A-sensitive burden 做 ablation。
