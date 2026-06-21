# Handoff: Branch B V2 Report-Layer Filtering

Date: 2026-06-21 17:30 Asia/Shanghai

Status: `active_current_handoff`

## 1. Goal

Make Branch B V2 perform actual report-layer filtering while keeping Branch A
as the sensitive WisecondorX/CBS discovery layer.

## 2. Restored Context

Context was restored from:

- `docs/handoff/2026-06-21_1614_branch_b_v2_report_contract_handoff.md`
- `AGENTS.md`
- `skills/conversation_handoff/SKILL.md`
- `skills/pgta_reference_modeling_analysis/SKILL.md`
- `CURRENT_STATE.md`
- `PLANS.md`
- `docs/CURRENT_CONTEXT_INDEX.md`

Active reference:

```text
h_r0_shadow_ref_20260619
```

Active Branch A overlay:

```text
merge_gap_bp=2_000_000
strong_z=10.0
postprocess.output_dir=wisecondorx/cnv/postprocess_gap2m
```

## 3. Completed Changes

Code and workflow changes:

- `pgta/predict/branch_b/v2_classifier.py`
  - Added report-layer classes:
    `report_event`, `internal_review_event`, `filtered_event`,
    `branch_s_event`.
  - Added report-layer visibility and filter-reason fields.
  - Added conservative combination filtering:
    non-strong A plus B-side non-support plus GC/RC attenuation plus
    sensitive/short/low-clean context.
  - Added report-event routing for strong autosomal, B-supported, clean,
    reportable-length candidates.
- `pgta/predict/branch_b/v2_benchmark.py`
  - Counts report/internal/filtered/Branch-S event classes.
  - Writes `filtered_events.tsv` and `filtered_events.json`.
  - Treats report-layer filtered truth-overlap candidates as FN.
- `rules/predict_layout.smk`, `rules/predict_workflow.smk`,
  `rules/target_assembly.smk`, `Snakefile`
  - Wire filtered-event benchmark outputs into the workflow target.
- `pgta/predict/report.py`
  - Loads V2 report-layer counts.
  - Uses V2 report-layer counts as the main biological conclusion when
    present.
  - Labels legacy/current-code Branch B top events as legacy context.
- Unit tests were updated for classifier, benchmark, workflow contract, and
  report behavior.

## 4. Current Results

Remote path:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419
```

Y1-Y8:

```text
candidates=97
truth=10
truth_preserved=10
FN=0
truth_hard_suppressed=0
truth_report_layer_filtered=0
report_event=21
internal_review_event=50
filtered_event=13
branch_s_event=13
```

H1-H16:

```text
candidates=105
truth=10
truth_preserved=10
FN=0
truth_hard_suppressed=0
truth_report_layer_filtered=0
report_event=6
internal_review_event=49
filtered_event=8
branch_s_event=42
H6 chr21 retained as internal_review_event
```

2026-06-15 context cohort:

```text
candidates=165
truth status=skipped_no_truth
report_event=52
internal_review_event=76
filtered_event=23
branch_s_event=14
```

Per-sample 2026-06-15:

```text
JZ26125843-56-56: report=3, internal_review=13, filtered=8, branch_s=2
JZ26125844-59-59: report=9, internal_review=19, filtered=6, branch_s=1
JZ26125845-60-60: report=22, internal_review=20, filtered=0, branch_s=5
JZ26125846-61-61: report=14, internal_review=15, filtered=1, branch_s=3
JZ26125847-62-62: report=4, internal_review=9, filtered=8, branch_s=3
```

## 5. Validation

Remote unit test:

```bash
cd /data/project/CNV/PGT-A/refactor_validation_20260419
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_v2_classifier.py \
  tests/unit/test_branch_b_v2_benchmark.py \
  tests/unit/test_cnv_report.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py -q
```

Result:

```text
66 passed in 1.15s
```

Remote dry-run:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile <active_gap2m_config> --cores 1 -n \
  branch_b_v2_benchmark branch_s_review cnv_report
```

Result:

```text
Passed for all three active gap2m configs.
```

Remote materialization:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile <active_gap2m_config> --cores 24 \
  --forcerun cnv_branch_b_v2_classifier cnv_branch_b_v2_benchmark \
  cnv_report_summary branch_b_v2_benchmark cnv_report
```

Result:

```text
Y1-Y8: 14/14 steps done.
H1-H16: 22/22 steps done.
2026-06-15: 11/11 steps done.
```

## 6. Current Interpretation

This loop adds real report-layer filtering. It is stronger than the previous
burden display-only loop.

It is still not final production promotion:

- 2026-06-15 has no truth labels, so it cannot prove FP precision.
- Branch S is still `review_reportable_with_limitations`, not final SCA.
- The active reference is still a fixed shadow baseline.
- 2026-06-15 still has many `report_event` rows, especially
  JZ26125845-60-60.

## 7. Next Recommended Step

Review the new V2 `report_event` rows, especially in 2026-06-15, and design a
second report-layer demotion pass that can move weakly supported report events
to internal review without filtering locked truth or H6 chr21.

Do not use 2026-06-15 to tune truth thresholds directly; it remains burden and
context only.

## 8. Key Files

- `docs/reports/branch_b_v2_report_layer_filter_2026-06-21.md`
- `pgta/predict/branch_b/v2_classifier.py`
- `pgta/predict/branch_b/v2_benchmark.py`
- `pgta/predict/report.py`
- `rules/predict_layout.smk`
- `rules/predict_workflow.smk`
- `rules/target_assembly.smk`
- `Snakefile`

## 9. Core File Sync

- `CURRENT_STATE.md`: updated.
- `PLANS.md`: updated.
- `docs/CURRENT_CONTEXT_INDEX.md`: updated.
- `AGENTS.md`: not updated; no repository-level hard rule changed.
