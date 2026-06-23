# Repository Archive Inventory

Date: 2026-06-23

Status: conservative_archive_inventory_remote_validated

## Context

This archive inventory follows the validated predict BAM QC and sample
reportability QC loop documented in:

- `docs/handoff/2026-06-23_1949_predict_bam_sample_reportability_qc_handoff.md`
- `docs/reports/predict_bam_and_sample_reportability_qc_2026-06-23.md`

It does not change Branch A, Branch B V2, Branch S, reference build, sex
calling, event classification, or report-event filtering.

## Keep In Main

These files and directories remain part of the active workflow or current
context contract:

- `Snakefile`
- `rules/`
- `pgta/`
- `scripts/_compat_entry.py` and formal script entrypoints
- `tests/unit/`
- `tests/server_validation/`
- active `config*.yaml` and cohort YAML files
- `sample_info/*.tsv` truth/cohort files used by workflow tests
- `docs/CURRENT_CONTEXT_INDEX.md`
- `CURRENT_STATE.md`
- `PLANS.md`
- current reports and handoffs for the active R&D loop

## Archive Only

These should be retained as historical evidence and indexed through
`docs/CURRENT_CONTEXT_INDEX.md`, not treated as active implementation state:

- older `docs/handoff/*.md` files
- older `docs/reports/*.md` R&D reports
- historical ablation reports that are explicitly superseded
- historical plot-development reports from 2026-06-22

No historical docs are deleted in this pass because they are useful for audit
trail and many are referenced by current context tests.

## Delete Or Ignore Candidates

These are local/generated artifacts and should not enter Git:

- `/reports/` local plot/report copies
- local Excel metadata under `/sample_info/*.xlsx`
- Python caches and test caches
- root-level temporary shell/status scripts matching `check_*`, `start_*`,
  `kill_*`, `inspect_*`, or `remote_status_*`
- one-off local validation helpers such as `validate_npz.py`
- copied reference packages under `/reference_packages/`

Current scan did not find root-level temporary shell scripts or
`validate_npz.py`. Local `/reports/` is about 62 MB and is kept for now because
it contains recently reviewed plots; it remains ignored and can be deleted or
re-synchronized later.

## Reference Checks

Before deleting any tracked file or remote result, run reference checks against:

```text
Snakefile
rules/
pgta/
scripts/
cli/
tests/
docs/
*.yaml
```

The first pass used `rg` to inspect likely temporary names and report paths.
No workflow dependency required deletion. The cleanup therefore remains
conservative and non-destructive for tracked workflow files and remote results.

## Current Status Sync

`docs/CURRENT_CONTEXT_INDEX.md` now records:

```text
report_status: sample_reportability_qc_remote_validated_0615_materialized
```

This resolves the stale `pending_remote_validation` wording that conflicted
with the 2026-06-23 QC handoff.

## Remote Validation

Remote validation was executed on `ssh fengxian` in:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419
```

Pytest:

```text
tests/unit/test_current_context_index.py
tests/unit/test_branch_b_plot.py
tests/unit/test_predict_bam_qc.py
tests/unit/test_sample_report_qc.py
tests/unit/test_cnv_report.py
tests/unit/test_branch_ab_phase12_workflow_contract.py

73 passed in 4.99s
```

Snakemake dry-run passed for the Y/H/G/0615 active gap2m lowres configs with
targets:

```text
predict_bam_qc cnv_sample_report_qc branch_b_v2_report_ablation branch_s_review cnv_report
```

Materialization was rerun for the same four active configs. The first attempt
was interrupted by the local SSH timeout during H-batch report generation; it
was rerun with a longer timeout and completed. Y and 0615 were already up to
date; H and G report outputs were refreshed successfully.

## Next Recommended Gate

The next non-archive development gate is to decide whether BAM-only predict
batches should keep `FASTP_METRICS_MISSING` as `BAM_QC_REVIEW` or gain a
BAM-derived GC/insert-size substitute.
