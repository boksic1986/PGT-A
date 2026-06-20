# 2026-06-21 02:35 Remote Snakemake Parse Repair Handoff

## 1. Goal

Repair the remote mirror / Snakemake parse blocker before continuing Branch
A/B/S workflow validation.

## 2. Context Read

- `docs/CURRENT_CONTEXT_INDEX.md`
- `docs/handoff/2026-06-21_0212_context_index_rnd_contract_handoff.md`
- `AGENTS.md`
- `skills/conversation_handoff/SKILL.md`
- `skills/pgta_reference_modeling_analysis/SKILL.md`
- `CURRENT_STATE.md`
- `PLANS.md`
- `Snakefile`

Active handoff used:

```text
docs/handoff/2026-06-21_0212_context_index_rnd_contract_handoff.md
```

## 3. Root Cause

Remote Snakemake 8.5.4 failed to parse workflow files copied from the Windows
worktree with CRLF line endings.

Evidence:

- local and remote `Snakefile` content SHA256 matched;
- remote byte inspection showed `\r\n`;
- an LF-only temporary `Snakefile` bypassed the original parser error;
- the same parser error then appeared in `rules/pipeline_modes.smk`, proving the
  issue was not a single-file syntax defect.

## 4. Completed

Added:

```text
.gitattributes
tests/unit/test_workflow_line_endings_contract.py
docs/reports/remote_snakemake_parse_repair_2026-06-21.md
docs/handoff/2026-06-21_0235_remote_snakemake_parse_repair_handoff.md
```

Updated:

```text
docs/CURRENT_CONTEXT_INDEX.md
CURRENT_STATE.md
PLANS.md
```

Remote mirror normalized to LF for:

```text
Snakefile
root YAML configs
rules/**/*.smk
pgta/**/*.py
scripts/**/*.py
tests/**/*.py
```

## 5. Verification

Remote CRLF scan:

```text
executable: ssh fengxian
command: scan Snakefile, root YAML, rules, pgta, scripts, and tests for CRLF
result: crlf_file_count=0
```

Remote tests:

```text
executable: ssh fengxian
command: cd /data/project/CNV/PGT-A/refactor_validation_20260419 && PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 /biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest tests/unit/test_workflow_line_endings_contract.py tests/unit/test_current_context_index.py -q
result: 5 passed in 0.16s
```

Remote dry-run that previously failed:

```text
executable: ssh fengxian
command: cd /data/project/CNV/PGT-A/refactor_validation_20260419 && /biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile config_predict_y_h_r0_shadow_branch_b_evidence_20260620.yaml --cores 1 -n branch_b_evidence branch_s_review cnv_report
result: passed; DAG built and dry-run planned cnv_report_summary and cnv_report
```

## 6. Not Changed

- No CNV calling logic.
- No mosaic logic.
- No sex calling logic.
- No result schema.
- No report data.
- No BAM/NPZ/result files.

## 7. Current Conclusion

The remote Snakemake parse blocker is repaired. This restores dry-run
capability for the active remote mirror, but it does not validate any new Branch
A/B/S algorithm behavior.

## 8. Recommended Next Step

Proceed to Branch A burden optimization under `h_r0_shadow_ref_20260619`, with:

- Y1-Y8/H1-H6 FN=0 retained;
- H6 chr21 retained;
- no Branch B artifact filter added into Branch A;
- all workflow validation on `ssh fengxian`.

## 9. Core File Sync

- `AGENTS.md`: not updated; no repository hard-rule change.
- `REPO_MAP.md`: not updated; no stable entrypoint change.
- `PLANS.md`: updated; remote parse repair is no longer the active blocker.
- `CURRENT_STATE.md`: updated; root cause and repair state recorded.
