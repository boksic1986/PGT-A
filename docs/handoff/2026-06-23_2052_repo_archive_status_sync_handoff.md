# Repository Archive And Status Sync Handoff

Date: 2026-06-23 20:52

Status: conservative_archive_inventory_remote_validated

## 1. Goal

Close the current stable-flow loop by recording a conservative repository
archive inventory, syncing stale context status, and keeping report/QC
readiness separate from Branch A/B/S algorithm changes.

## 2. Confirmed Project State

- Active prior handoff:
  `docs/handoff/2026-06-23_1949_predict_bam_sample_reportability_qc_handoff.md`
- Current reference remains `h_r0_shadow_ref_20260619`.
- Reference status remains `fixed_shadow_baseline_not_production`.
- Branch A, Branch B V2, Branch S, sex calling, event classification, and
  report-event filtering are unchanged by this cleanup loop.

## 3. Completed Items

- Added `docs/reports/repo_archive_inventory_2026-06-23.md`.
- Updated `docs/CURRENT_CONTEXT_INDEX.md` so `report_status` matches the
  validated 2026-06-23 QC/reportability handoff.
- Updated context-index unit expectations to the current handoff chain.
- Classified repository artifacts into keep, archive-only, and
  delete-or-ignore candidates without deleting tracked workflow files or remote
  result directories.

## 4. Current Conclusion

The repository should be archived conservatively. Historical handoffs and
reports are retained as audit evidence; local generated reports, Excel files,
cache directories, and one-off scripts remain ignored or cleanup candidates.

## 5. Verification

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_current_context_index.py \
  tests/unit/test_branch_b_plot.py \
  tests/unit/test_predict_bam_qc.py \
  tests/unit/test_sample_report_qc.py \
  tests/unit/test_cnv_report.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py \
  -q
```

Result:

```text
73 passed in 4.99s
```

Active dry-run targets were run for the Y/H/G/0615 gap2m lowres configs:

```text
predict_bam_qc cnv_sample_report_qc branch_b_v2_report_ablation branch_s_review cnv_report
```

Result:

```text
dry-run passed for all four active configs
```

Materialization was rerun for the same targets. The first run was interrupted
by the local SSH timeout during H-batch `cnv_report_summary`, which caused
Snakemake to delete the partial H report outputs. The run was repeated with a
longer timeout and completed; H and G reports were refreshed, and Y/0615 were
already up to date.

## 6. Recommended Next Step

After merge/sync, decide whether BAM-only predict batches should keep
`FASTP_METRICS_MISSING` as `BAM_QC_REVIEW` or implement BAM-derived
GC/insert-size/coverage substitutes.

## 7. Key Files

- `docs/reports/repo_archive_inventory_2026-06-23.md`
- `docs/CURRENT_CONTEXT_INDEX.md`
- `CURRENT_STATE.md`
- `PLANS.md`
- `docs/handoff/2026-06-23_1949_predict_bam_sample_reportability_qc_handoff.md`

## 8. Environment Constraints

- All real workflow validation must run on `ssh fengxian`.
- Use `/biosoftware/miniconda/envs/snakemake_env/bin/python`.
- Use `/biosoftware/miniconda/envs/snakemake_env/bin/snakemake`.
- Do not delete remote result directories in this archive pass.

## 9. Command Record

Local and remote commands:

- `git status --short --branch`
- `git diff -- pgta/predict/branch_b/plot.py tests/unit/test_branch_b_plot.py`
- `rg` scans for temporary script names and report references
- remote pytest listed above
- remote Snakemake dry-run listed above
- remote Snakemake materialization listed above

## 10. Core File Sync

- `AGENTS.md`: not updated; hard rules unchanged.
- `REPO_MAP.md`: not updated; stable structure unchanged.
- `PLANS.md`: updated in this loop.
- `CURRENT_STATE.md`: updated in this loop.
