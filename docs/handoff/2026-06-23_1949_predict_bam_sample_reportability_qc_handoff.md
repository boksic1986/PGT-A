# Predict BAM QC And Sample Reportability QC Handoff

Date: 2026-06-23 19:49

Status: implemented_remote_validated_0615_materialized

## Context

Active prior handoff:

```text
docs/handoff/2026-06-22_2235_report_table_ablation_audit_handoff.md
```

This handoff records the new QC/reportability layer added on top of the current
Branch A/B/S report path.

## Completed Implementation

- Added `pgta.predict.bam_qc`.
- Added `pgta.predict.sample_report_qc`.
- Added dispatcher actions:
  - `predict_bam_qc`
  - `cnv_sample_report_qc`
- Added Snakemake targets:
  - `predict_bam_qc`
  - `cnv_sample_report_qc`
- Added report inputs:
  - `--predict-bam-qc-summary`
  - `--sample-report-qc`
- Added report display fields for sample-level reportability QC.
- Added unit-test contracts for BAM QC, sample reportability QC, report display,
  and workflow target wiring.

## Design Boundary

This change does not modify:

- Branch A discovery;
- Branch B V2 classifier/filter/report-event decisions;
- Branch S classifier;
- reference build;
- sex calling;
- event-level TP/FN/FP metrics.

Sample QC can mark `SAMPLE_QUALITY_REVIEW` or `NO_CALL_RECOMMENDED`, but it does
not reclassify individual CNV/SCA events.

## Validation

Local:

```text
git diff --check: passed
python -m py_compile pgta/predict/bam_qc.py pgta/predict/sample_report_qc.py pgta/predict/report.py scripts/_compat_entry.py: passed
```

Local pytest could not run because the local Python environment lacks `pytest`
and `pandas`.

Remote validation completed on:

```text
ssh fengxian
/data/project/CNV/PGT-A/refactor_validation_20260419
```

Remote pytest:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_predict_bam_qc.py \
  tests/unit/test_sample_report_qc.py \
  tests/unit/test_cnv_report.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py \
  -q

42 passed in 1.25s
```

Dry-run active Y/H/G/0615 configs passed:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile <active-config> --cores 1 -n \
  predict_bam_qc cnv_sample_report_qc cnv_report
```

0615 materialized:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml \
  --cores 8 predict_bam_qc cnv_sample_report_qc cnv_report
```

0615 result:

- `NO_CALL_RECOMMENDED`: 0 samples.
- `SAMPLE_QUALITY_REVIEW`: 5 samples.
- `JZ26125845-60-60` is the expected high-wave outlier and uniquely triggered
  `HIGH_MAD_LOG1P_RELATIVE_OUTLIER` plus
  `AUTOSOMAL_CN_OUTSIDE_1P7_2P3_REVIEW`.
- 56/59/61/62 are not no-called; their review state is driven by missing fastp
  library context and report burden warnings only.

## Next Step

Decide whether BAM-only predict batches should keep `FASTP_METRICS_MISSING` as
`BAM_QC_REVIEW`, or whether a future BAM-derived GC/insert-size substitute is
needed to avoid flagging every BAM-only sample as library review.
