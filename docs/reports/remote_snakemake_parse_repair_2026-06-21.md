# Remote Snakemake Parse Repair

Date: 2026-06-21

Status: `active_current_evidence`

Decision use: records the repair of the remote Snakemake parse blocker before
continuing Branch A/B/S workflow validation.

## Problem

Remote Snakemake dry-runs on `ssh fengxian` failed before DAG construction:

```text
IndentationError in file /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile, line 9:
expected an indented block after function definition on line 10 (Snakefile, line 9)
```

The failure blocked new workflow-validation claims.

## Root Cause Evidence

Local and remote `Snakefile` contents had the same SHA256:

```text
1c6a607867e6bf8d946f500016eee7ac22ed3d8606c5234a6cfd3ffb984be76d
```

Remote byte inspection showed CRLF endings:

```text
configfile: "build_samples.yaml"\r$
...
def merge_config(base, override):\r$
    merged = dict(base)\r$
```

A temporary LF-only copy of `Snakefile` passed the original first-line parser
failure and advanced to later workflow parsing. Running the original dry-run
with an LF `Snakefile` then exposed the same CRLF parser failure in
`rules/pipeline_modes.smk`, confirming the issue affected multiple workflow
source files.

Root cause:

```text
Windows CRLF workflow files were copied into the Linux remote mirror.
Snakemake 8.5.4 on the remote mirror fails to parse these Snakefile/.smk
sources correctly.
```

## Repair

Repository contract added:

```text
.gitattributes
Snakefile text eol=lf
*.smk text eol=lf
*.py text eol=lf
*.yaml text eol=lf
*.yml text eol=lf
```

Remote mirror files normalized to LF under:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419
```

Normalized file classes:

```text
Snakefile
root *.yaml / *.yml
rules/**/*.smk
pgta/**/*.py
scripts/**/*.py
tests/**/*.py
```

The normalization changed 97 remote workflow/config/source files. No result
data, BAM, NPZ, report package, Excel sheet, or image file was touched.

## Verification

Remote CRLF check:

```text
executable: ssh fengxian
command: scan Snakefile, root YAML, rules, pgta, scripts, and tests for CRLF
result: crlf_file_count=0
```

Remote unit tests:

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

## Current Interpretation

The remote parse blocker is repaired. This does not validate new Branch A/B/S
algorithm behavior; it only restores the ability to run Snakemake dry-runs from
the active remote mirror.

The next workflow work can proceed to Branch A burden optimization under
`h_r0_shadow_ref_20260619`, preserving Y1-Y8/H1-H6 no-FN and H6 chr21 detection.
