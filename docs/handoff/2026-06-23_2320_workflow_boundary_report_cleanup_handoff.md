# Workflow Boundary, Report Style, And Conservative Cleanup Handoff

Date: 2026-06-23

Status: remote_validated_cleaned

## Context

This handoff supersedes
`docs/handoff/2026-06-23_2052_repo_archive_status_sync_handoff.md` for the
current cleanup/report-boundary loop.

Active reference remains `h_r0_shadow_ref_20260619`.

No Branch A, Branch B V2, Branch S, reference build, sex calling, report-event
classification, or filtering logic was intentionally changed.

## Completed In This Loop

- Added production/test boundary contract coverage.
- Locked active report configs so `branch_b_v2_report_ablation` remains an
  explicit R&D target rather than a default `cnv_report` dependency.
- Locked `tests/server_validation/*.sh` as validation wrappers, not production
  DAG includes.
- Updated z/CN report plot style to the accepted 0615 negative-sample preview
  layout:
  - SVG `3200 x 700`;
  - title `24 bold`;
  - y-axis title `32 bold`;
  - chromosome labels `22 bold`;
  - position ticks `15 bold`, light gray, numbers only;
  - a single leading `Mb` unit label;
  - white background plus light chromosome blocks.
- Added `docs/reports/production_workflow_test_boundary_audit_2026-06-23.md`.

## Validation State

Initial TDD remote check on `fengxian`:

- New plot-style tests failed against the previous implementation because SVG
  height was still `620` and axis labels/ticks used old small fonts.
- After implementation, the targeted plot-style and workflow-boundary tests
  passed: `3 passed in 1.00s`.

Full remote validation after cleanup:

- pytest:
  `tests/unit/test_branch_ab_phase12_workflow_contract.py`,
  `tests/unit/test_branch_b_plot.py`, `tests/unit/test_cnv_report.py`, and
  `tests/unit/test_current_context_index.py`: `67 passed in 4.57s`.
- Active Y/H/G/0615 dry-runs for
  `predict_bam_qc cnv_sample_report_qc branch_s_review cnv_report`: all passed
  and reported requested outputs already present.
- Explicit Y/H/G/0615 dry-runs for `branch_b_v2_report_ablation`: all passed
  and reported requested outputs already present.

## Cleanup Boundary

Allowed local cleanup:

- ignored `reports/`;
- `.pytest_cache`;
- `__pycache__`;
- transient preview outputs.

Allowed remote cleanup:

- `.pytest_cache`;
- `__pycache__`;
- `tests/qc_result`.

Explicitly not deleted remotely:

- `reports/`;
- `.snakemake/`;
- `results_*`;
- workflow result packages.

Actual cleanup:

- local ignored `/reports/` removed: 79 files, 62,238,457 bytes;
- local caches removed;
- remote `.pytest_cache`, recursive `__pycache__`, and `tests/qc_result`
  removed;
- remote `reports/` (`135M`) and `.snakemake/` (`99M`) retained.

## Next Commands

1. Commit the current code, tests, reports, and this handoff.
2. Open PR from `codex/workflow-boundary-report-cleanup` to `main`.
3. After review/merge, sync tracked `main` content to the non-git remote mirror.

## Risks

- `branch_b_v2_report_ablation` is still R&D-only. It must not be interpreted as
  a production report-table filtering step.
- Plot style changes are display-only. They must not be used as evidence that
  CNV/SCA calls changed.
- Cleanup must remain conservative; deleting remote result packages would break
  reproducibility and is outside this pass.
