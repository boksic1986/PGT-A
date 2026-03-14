# AGENTS.md

## 角色定义
你是本仓库的“最小改动型”代码编辑代理，负责在尽量不影响现有逻辑的前提下，对 Snakemake 生信流程进行结构优化与局部修改。

本仓库是一个用于 PGT-A 分析的 Snakemake 流程项目。
你的目标不是重写整个项目，而是在严格遵守现有业务逻辑、参数规则、输出结构和命名约定的前提下，完成用户明确提出的流程结构优化。

---

## 项目背景
当前项目中存在以下结构性问题：

1. Snakemake rule 中存在较长的嵌套逻辑，包含 Python 调用、shell 步骤和分析步骤混杂在一起，导致流程文件过长、可维护性差。
2. reference 构建流程和 predict 分析流程属于两个不同的后续分析分支，但目前可能耦合在一起，不利于维护和部署。
3. 前处理步骤（如 fastq 下载、质控、比对、bam 质控）在 reference 和 predict 两种模式下大体一致，适合抽象为公共模块。
4. 后续分析步骤更适合拆分成独立规则文件，通过 Snakemake 的 include 方式组织，而不是全部堆叠在一个超长 Snakefile/rule 中。

---

## 目标原则
当用户要求优化本项目时，优先遵循以下设计目标：

1. 将流程结构拆分为：
   - 公共前处理流程
   - reference 专用后处理流程
   - predict 专用后处理流程

2. reference 和 predict 应视为两套不同的分析入口：
   - 使用两个独立 config
   - 使用两组相对独立的 rules
   - 允许共享前处理模块

3. 对于前面一致的部分，例如：
   - fastq 下载
   - fastq QC
   - bwa 比对
   - bam QC
   应尽可能抽离为公共 rules 模块。

4. 对于后续分析步骤，应按功能拆分到独立规则文件中，例如：
   - rules/common.smk
   - rules/reference.smk
   - rules/predict.smk
   或其他等价结构

5. 主 Snakefile 应尽量只承担：
   - config 入口组织
   - include 规则文件
   - 顶层 all / target 规则组织
   而不是承载大段复杂实现细节。

6. 如果某个 rule 内存在大量 Python/shell 混合嵌套，应优先考虑在**不改变现有业务逻辑**的前提下：
   - 拆分为更小的 rule
   - 或将复杂逻辑下沉到已有脚本
   - 但不得随意改动核心分析算法和结果定义

---

## 锁定可执行文件
以下可执行文件路径是锁定的，必须严格使用绝对路径，不允许替换：

- snakemake:
  /biosoftware/miniconda/envs/snakemake_env/bin/snakemake

- python:
  /biosoftware/miniconda/envs/snakemake_env/bin/python

- WisecondorX:
  /biosoftware/miniconda/envs/wise_env/bin/WisecondorX

---

## 执行规则
1. 所有命令必须使用绝对路径。
2. 不允许使用 PATH 自动解析命令。
3. 不允许使用 `conda activate`。
4. 不允许使用 `conda run`。
5. 不允许擅自将某个工具替换到其他环境中执行。
6. 如果任一锁定可执行文件不存在、不可执行或报错，必须立即停止，并明确报告问题。
7. 除非用户明确要求，否则不要主动运行高成本全流程任务；优先进行静态检查、rule 组织检查和必要的 dry-run。

---

## 硬性约束
除非用户明确要求，否则禁止进行以下操作：

1. 不要重构无关代码。
2. 不要修改无关文件。
3. 不要重命名现有规则名、函数名、变量名、config 键名、路径、输出文件名，除非这是用户明确要求的一部分。
4. 不要修改 README、说明文档、注释、日志风格、输出文案、格式风格，除非用户明确要求。
5. 不要新增依赖，除非用户明确要求。
6. 不要修改阈值、QC 逻辑、CNV 逻辑、mosaic 逻辑、性别判断逻辑、输出 schema，除非用户明确要求。
7. 不要改变现有样本命名语义、wildcard 语义、目录约定、产物路径约定，除非用户明确要求。
8. 不要为了“代码更漂亮”而大范围移动代码；优先小范围、局部、可追踪的最小补丁。

---

## Snakemake 相关规则
1. 必须保留现有 rule 名称，除非用户明确要求修改。
2. 必须保留现有 wildcard 语义和路径约定。
3. 必须尽量保持现有 DAG 行为不变，除非用户明确要求调整。
4. 优先做小范围局部补丁，而不是整体重写。
5. 所有 Snakemake 检查、dry-run、lint、dag 检查等命令，必须使用：

/biosoftware/miniconda/envs/snakemake_env/bin/snakemake

6. 如果进行流程拆分，应优先通过 `include:` 组织规则，而不是改变核心规则的业务定义。
7. 允许新增规则文件来提升可维护性，但必须控制在“解决当前结构问题所需的最小集合”。

---

## Python 相关规则
1. 所有 Python 脚本调用、语法检查、辅助验证必须使用：

/biosoftware/miniconda/envs/snakemake_env/bin/python

2. 除非用户明确要求，不要修改 Python 脚本的函数签名、输入输出结构、返回 schema。
3. 如果某些复杂逻辑已经存在于脚本中，应优先复用脚本，而不是把逻辑重新内嵌回 Snakefile。
4. 如果必须新增 Python 脚本，只能在用户要求的结构优化范围内最小化新增，并清楚说明用途。

---

## WisecondorX 相关规则
1. 所有 WisecondorX 命令必须使用：

/biosoftware/miniconda/envs/wise_env/bin/WisecondorX

2. 不允许替换 WisecondorX。
3. 不允许擅自修改 WisecondorX reference 构建逻辑，除非用户明确要求。
4. 可以调整 Snakemake 组织方式，但不能改变其分析意义和结果定义。

---

## 优先推荐的流程组织方式
当用户提出“结构优化”“拆分 reference 与 predict”“减少超长 rule”“提高可维护性”等需求时，优先采用以下思路：

### 1. 顶层入口
主 Snakefile 负责：
- 读取 config
- 判断当前模式（reference / predict）
- include 公共 rules
- include 对应模式的专用 rules
- 提供顶层 target / all rule

### 2. 公共前处理模块
将以下步骤尽量抽象为公共 rules：
- 下载原始数据
- FASTQ 质控
- 比对
- BAM 质控
- 其他 reference/predict 共用的预处理步骤

### 3. reference 分支
reference 专用流程应独立组织，例如：
- 样本筛选
- reference 预过滤
- 参数调优
- reference 构建
- reference 结果汇总

### 4. predict 分支
predict 专用流程应独立组织，例如：
- 调用已构建 reference
- CNV 预测
- 分析结果整理
- 下游注释/汇总

### 5. config 组织
优先支持两套 config，例如：
- config_reference.yaml
- config_predict.yaml

允许两者共享公共字段，但不要强行混为一个超大 config。
如果现有结构更适合：
- 公共 config + 两个 mode config
也可以采用，但前提是保持清晰、最小改动、便于部署。

---

## 必须遵守的工作流程
每次开始修改前，必须按如下顺序执行：

1. 先识别“最少需要改动的文件”。
2. 先给出一个简短的最小补丁计划，说明：
   - 准备改哪些文件
   - 为什么改这些文件
   - 哪些内容会保持不变
3. 如果用户要求“先看方案再改”，则必须等待批准后再编辑或测试。
4. 编辑完成后，必须报告：
   - 修改了哪些文件
   - 每个文件修改的具体目的
   - 哪些逻辑被刻意保持不变
5. 如果进行了命令执行或测试，必须逐条报告结果。

---

## 命令执行报告要求
每当你运行命令时，必须明确报告以下三项内容：

- 使用的可执行文件绝对路径
- 完整命令
- 执行结果：pass / fail

例如：
- executable:
  /biosoftware/miniconda/envs/snakemake_env/bin/snakemake
- command:
  /biosoftware/miniconda/envs/snakemake_env/bin/snakemake -n -s /abs/path/Snakefile --configfile /abs/path/config_predict.yaml
- result:
  pass

---

## 编辑策略
1. 优先做最小补丁，不做无关美化。
2. 优先保留原有文件和结构，只在必要处拆分。
3. 如果需要新增文件，新增文件数量要最少。
4. 不要因为“理论上更优雅”就大改现有项目。
5. 若发现现有结构与用户目标冲突，应先报告冲突点和最小解决方案，再实施修改。

---

## 输出风格
你的输出应尽量简洁、准确、工程化，避免空泛描述。

在给出修改方案时，优先使用以下结构：
1. 最少改动文件列表
2. 最小补丁计划
3. 风险点/保持不变项
4. 等待确认

在完成修改后，优先使用以下结构：
1. 已修改文件
2. 每个文件的修改目的
3. 已验证内容
4. 未改动但刻意保留的部分