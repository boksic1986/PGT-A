# Production Workflow Test Boundary Audit

Date: 2026-06-23

Status: remote_validated_cleaned

## Scope

This audit fixes the production/test boundary for the current PGT-A CNV-seq
workflow and locks the formal z/CN report plot style. It does not modify Branch
A, Branch B V2, Branch S, reference build, sex calling, event classification, or
report-event filtering.

## Read Sources

- `docs/handoff/2026-06-23_2052_repo_archive_status_sync_handoff.md`
- `AGENTS.md`
- `skills/conversation_handoff/SKILL.md`
- `skills/pgta_reference_modeling_analysis/SKILL.md`
- `Snakefile`
- `rules/pipeline_modes.smk`
- `rules/target_assembly.smk`
- `rules/predict_workflow.smk`
- `rules/predict_layout.smk`
- `scripts/README.md`
- `cli/pgta.py`
- `tests/server_validation/*.sh`
- active Y/H/G/0615 gap2m lowres predict configs
- `pgta/predict/branch_b/plot.py`

## Production Entry Points

The production workflow remains:

```text
Snakefile
rules/*.smk
scripts/*.py stable dispatchers
pgta/* implementation modules
active config*.yaml files
```

`scripts/*.py` are retained because they are stable dispatcher/wrapper entry
points used by Snakemake and the CLI. They are not one-off patch scripts.

## Validation Entry Points

`tests/server_validation/*.sh` are validation wrappers. They are referenced by
`cli/pgta.py`, `README.md`, and test-plan documentation, but they are not
included by `Snakefile` and are not production DAG includes.

The server validation wrappers are therefore kept in main, but treated as
validation-only utilities.

## Ablation Boundary

`branch_b_v2_report_ablation` is a formal explicit R&D audit target. It remains
available for manual/dedicated validation, but active production report configs
must not request it by default.

`cnv_report` may consume current V2 report events, QC summaries, plot outputs,
and Branch S summaries. It must not automatically consume ablation-only outputs
or truth-only conclusions.

Locked contract:

```text
active report pipeline.targets: no branch_b_v2_report_ablation
cnv_report target assembly: no CNV_B_V2_REPORT_ABLATION_* outputs
tests/server_validation: not a Snakefile include
```

## Formal Report Plot Style

The formal z/CN SVG plot style is pinned to the last accepted negative-sample
preview format:

```text
SVG width: 3200
SVG height: 700
title: 24 bold
y-axis title: 32 bold
chromosome label: 22 bold
position tick: 15 bold, #94a3b8
position unit: one leading "Mb" label, not repeated in every tick
background: white page with light alternating chromosome blocks
chromosome separators: no vertical separator lines
```

This is a display contract only. It does not change plotted data, CNV calls, SCA
classification, filtering, or report-event counts.

## Cleanup Policy

Local cleanup candidates:

- ignored `/reports/`
- `.pytest_cache`
- `__pycache__`
- transient preview outputs

Local files retained:

- `sample_info/*.xlsx` remain ignored but not deleted by this pass.
- tracked TSV truth/sample files remain in Git.

Remote cleanup candidates:

- `.pytest_cache`
- `__pycache__`
- `tests/qc_result`

Remote files explicitly retained:

- `reports/`
- `.snakemake/`
- `results_*`
- workflow result packages

## Verification Plan

Remote validation must run on `fengxian` with fixed executables:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/python
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake
```

Required pytest:

```text
tests/unit/test_branch_ab_phase12_workflow_contract.py
tests/unit/test_branch_b_plot.py
tests/unit/test_cnv_report.py
tests/unit/test_current_context_index.py
```

Required dry-runs:

```text
predict_bam_qc cnv_sample_report_qc branch_s_review cnv_report
branch_b_v2_report_ablation
```

for the active Y/H/G/0615 gap2m lowres configs.

## Validation Results

Remote validation was executed on `fengxian` in:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419
```

Pytest command:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_ab_phase12_workflow_contract.py \
  tests/unit/test_branch_b_plot.py \
  tests/unit/test_cnv_report.py \
  tests/unit/test_current_context_index.py -q
```

Result:

```text
67 passed in 4.57s
```

Active Y/H/G/0615 configs were dry-run with:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile <active_config> --cores 1 -n \
  predict_bam_qc cnv_sample_report_qc branch_s_review cnv_report
```

All four reported:

```text
Building DAG of jobs...
Nothing to be done (all requested files are present and up to date).
```

The explicit R&D audit target was also dry-run for all four configs:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile <active_config> --cores 1 -n branch_b_v2_report_ablation
```

All four reported the same up-to-date dry-run status.

## Cleanup Results

Local cleanup:

- removed ignored `/reports/` copy, 79 files, 62,238,457 bytes;
- removed local `__pycache__` and `.pytest_cache` directories;
- retained ignored `sample_info/*.xlsx`;
- retained tracked sample/truth TSV files.

Remote cleanup:

- removed `.pytest_cache`, recursive `__pycache__`, and `tests/qc_result`;
- retained remote `reports/` (`135M`);
- retained remote `.snakemake/` (`99M`);
- retained remote `results_*` directories.

## Current State

Production/test boundary tests, report plot style tests, remote pytest, active
dry-runs, explicit ablation dry-runs, and conservative cleanup have completed.
