# P6 Report Package Contract

Date: 2026-06-20

Status: `engineering_artifact_current_contract`

Decision use: historical P6 package contract for the already materialized 0615
development report. Future P6/report remains the delivery target after Branch A
burden optimization and Branch B/S contract strengthening.

## Scope

This report records the P6 report-package integration under:

```text
reference_id: h_r0_shadow_ref_20260619
reference status: shadow-only
```

P6 is a development/report-boundary check. It does not promote the H-augmented
reference, Branch B, Branch S/SCA, or the 2026-06-15 five-sample report to final
release.

## Contract

The report package now carries P3/P5 status visibly while preserving the caller
boundary:

```text
WisecondorX predict / CBS -> Branch A candidates -> current Branch B final events
```

Additional report context is summary-only:

```text
Branch A validation summary JSON
P3 Branch B evidence summary JSON
P5 Branch S summary JSON
```

The report rule does not consume:

```text
P3 raw evidence ledger TSV
Branch S raw evidence TSV
Branch S score TSV
matched-negative percentile TSV
Branch B V2 classifier TSV
```

## Code Changes

- `pgta/predict/report.py`
  - Adds `--reference-id`, `--wisecondorx-predict-command`, and
    `--branch-a-validation-summary`.
  - Adds `--branch-b-evidence-summary` and `--branch-s-summary`.
  - Adds P3 evidence status and SCA report status columns to TSV/JSON/MD/HTML.
  - Adds `report_contract.status=development_only_not_final_release`.
  - Adds `branch_a_no_fn_status` and `same_reference_config_status` to the
    report contract.
- `pgta/predict/branch_b/evidence_ledger.py`
  - Adds `background_source_counts` and `background_status_counts` to P3
    summary JSON so report does not need raw ledger input.
- `rules/predict_workflow.smk`
  - Passes only `CNV_B_EVIDENCE_SUMMARY` and `CNV_BRANCH_S_SUMMARY` into
    `cnv_report_summary`.
- Unit tests updated for P6 summary-only report contract.

## Remote Validation

All commands below ran on `ssh fengxian` in:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419
```

### Unit tests

```bash
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_cnv_report.py \
  tests/unit/test_reference_candidate_audit.py \
  tests/unit/test_branch_a_validation.py \
  tests/unit/test_branch_b_evidence_ledger.py \
  tests/unit/test_branch_s_shadow.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py -q
```

Result:

```text
35 passed in 1.26s
```

### Snakemake dry-runs

Targeted dry-runs passed for:

```text
reference_audit: config_predict_h_reference_audit_20260620.yaml
branch_a_validation: config_predict_y_h_r0_shadow_branch_a_validation_20260620.yaml
branch_a_validation: config_predict_h_20260608_h_r0_shadow_branch_a_validation_20260620.yaml
branch_b_evidence: config_predict_20260615_h_r0_shadow_branch_b_evidence_20260620.yaml
branch_s_review: config_predict_20260615_h_r0_shadow_branch_s_review_20260620.yaml
cnv_report: config_predict_20260615_h_r0_shadow_20260619.yaml
```

The `cnv_report` dry-run includes `cnv_report_summary` and `cnv_report`, with
Branch A validation summaries plus P3/P5 summary JSON inputs.

### 2026-06-15 report package materialization

Command materialized the 2026-06-15 report package:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake \
  -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile \
  --configfile /data/project/CNV/PGT-A/refactor_validation_20260419/config_predict_20260615_h_r0_shadow_20260619.yaml \
  --cores 4 \
  cnv_report
```

Result:

```text
cnv_report: report completed: samples=5
```

## Materialized 2026-06-15 Output

```text
/data/project/CNV/PGT-A/refactor_validation_20260419/results_20260615_h_r0_shadow_ref_20260619/wisecondorx/cnv/report/cnv_summary.tsv
/data/project/CNV/PGT-A/refactor_validation_20260419/results_20260615_h_r0_shadow_ref_20260619/wisecondorx/cnv/report/cnv_summary.json
/data/project/CNV/PGT-A/refactor_validation_20260419/results_20260615_h_r0_shadow_ref_20260619/wisecondorx/cnv/report/cnv_summary.md
/data/project/CNV/PGT-A/refactor_validation_20260419/results_20260615_h_r0_shadow_ref_20260619/wisecondorx/cnv/report/cnv_summary.html
```

Report check:

```text
rows=5
report_contract.status=development_only_not_final_release
reference_id=h_r0_shadow_ref_20260619
wisecondorx_predict_command=WisecondorX predict <npz> <ref> <output-prefix> --zscore 5.0 --alpha 0.001 --maskrepeats 5 --minrefbins 150
branch_a_no_fn_status=passed_no_fn
same_reference_config_status=matched
branch_a_validation_summary_count=2
truth_detected_count=20
truth_event_count=20
FN_count=0
H6_chr21_status=detected
branch_b_evidence_summary_count=5
branch_s_summary_count=5
branch_b_raw_ledger_used=false
branch_s_raw_evidence_used=false
```

Per-sample P3/P5 status:

| sample | P3 candidates | P3 background | SCA mode | SCA status |
|---|---:|---|---|---|
| JZ26125843-56-56 | 28 | UNKNOWN_BACKGROUND=28 | review_development_only | development_only_not_final_reportable |
| JZ26125844-59-59 | 38 | UNKNOWN_BACKGROUND=38 | review_development_only | development_only_not_final_reportable |
| JZ26125845-60-60 | 57 | UNKNOWN_BACKGROUND=57 | review_development_only | development_only_not_final_reportable |
| JZ26125846-61-61 | 44 | UNKNOWN_BACKGROUND=44 | review_development_only | development_only_not_final_reportable |
| JZ26125847-62-62 | 34 | UNKNOWN_BACKGROUND=34 | review_development_only | development_only_not_final_reportable |

## Current Interpretation

P6 now has a workflow-generated report package that visibly carries:

- Branch A / Branch B report state;
- Branch A no-FN evidence under the same `h_r0_shadow_ref_20260619` reference;
- P3 Branch B evidence burden and background status;
- P5 Branch S/SCA review/development status;
- explicit report-contract release status.

This is still not a final report release. The 2026-06-15 report package is an
internal exploratory report package until final reference, Branch B, and Branch
S promotion gates pass.

## Remaining Release Blockers

- `h_r0_shadow_ref_20260619` is still a shadow reference.
- P3 background is `UNKNOWN_BACKGROUND`, not formal N0/cross-fit background.
- Branch S/SCA remains `review_development_only`.
- 2026-06-15 has no locked truth table and must not be used to claim validation
  performance.
- G6 final release requires an explicit release decision after reviewing the
  generated package and upstream gates.
