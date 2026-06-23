# Predict BAM QC And Sample Reportability QC

Date: 2026-06-23

Status: implemented_remote_validated_0615_materialized

## Scope

This change adds two QC layers to the predict/report path:

1. `predict_bam_qc`: library/BAM-level reportability evidence.
2. `cnv_sample_report_qc`: post-predict sample-level reportability evidence.

It does not modify Branch A, Branch B V2, Branch S, reference build, sex
calling, event classification, or TP/FN/FP definitions.

## Predict BAM QC

New target:

```text
predict_bam_qc
```

Outputs:

```text
wisecondorx/cnv/qc/predict_bam_qc/{sample}.bam_qc.tsv
wisecondorx/cnv/qc/predict_bam_qc/{sample}.bam_qc.json
wisecondorx/cnv/qc/predict_bam_qc/predict_bam_qc_summary.tsv
wisecondorx/cnv/qc/predict_bam_qc/predict_bam_qc_summary.json
```

Inputs and metrics:

- `samtools flagstat`: mapping rate and proper-pair rate.
- `samtools idxstats`: autosome mapped reads, autosomal chromosome-bias CV,
  chrX/chrY mapped read context.
- `fastp` JSON, if present: clean reads, GC content, duplication rate, insert
  peak.

`fastp` JSON is passed as a path parameter rather than a Snakemake input, so
existing-BAM predict runs do not trigger fastp/BWA remapping just for QC. If the
JSON is missing but BAM metrics are valid, the sample is marked
`BAM_QC_REVIEW`, not failed.

Status values:

```text
BAM_QC_PASS
BAM_QC_REVIEW
BAM_QC_FAIL
```

## Sample Reportability QC

New target:

```text
cnv_sample_report_qc
```

Outputs:

```text
wisecondorx/cnv/postprocess/sample_report_qc/sample_report_qc.tsv
wisecondorx/cnv/postprocess/sample_report_qc/sample_report_qc.json
```

Inputs and metrics:

- predict BAM QC summary;
- existing `cnv_qc` NPZ QC rows;
- Branch B V2 sample summary and report events;
- CN plot bin TSVs for autosomal CN outside `1.7-2.3` context.

Status values:

```text
PASS_REPORTABLE
SAMPLE_QUALITY_REVIEW
NO_CALL_RECOMMENDED
```

First-pass policy is review-first:

- input `cnv_qc=PASS` plus high wave/burden triggers
  `SAMPLE_QUALITY_REVIEW`;
- `BAM_QC_FAIL` or `cnv_qc=FAIL` triggers `NO_CALL_RECOMMENDED`;
- cohort size `<10` relative outliers are review only, not no-call;
- single coherent CNV, SCA, or chrY denominator artifact does not fail a
  sample by itself.

## Report Contract

`cnv_report` now accepts:

```text
--predict-bam-qc-summary
--sample-report-qc
```

The sample table includes reportability QC display fields and the report
contract includes pass/review/no-call sample counts.

This is a sample-level reportability warning layer. It does not suppress,
promote, or reclassify Branch A/B/S events.

## Validation State

Local static check:

```text
git diff --check: passed
python -m py_compile pgta/predict/bam_qc.py pgta/predict/sample_report_qc.py pgta/predict/report.py scripts/_compat_entry.py: passed
```

Local functional pytest was not run because the local Python environment lacks
`pytest` and `pandas`; final validation was run on `ssh fengxian`.

Remote unit validation:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_predict_bam_qc.py \
  tests/unit/test_sample_report_qc.py \
  tests/unit/test_cnv_report.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py \
  -q

42 passed in 1.25s
```

Remote dry-run:

Y1-Y8, H1-H16, G1-G8, and 2026-06-15 active gap2m lowres configs all parsed
and dry-ran successfully for:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile <active-config> --cores 1 -n \
  predict_bam_qc cnv_sample_report_qc cnv_report
```

The dry-runs scheduled only predict BAM QC, sample reportability QC, and report
refresh jobs. They did not schedule reference rebuild.

## 0615 Materialized Check

Materialized command:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml \
  --cores 8 predict_bam_qc cnv_sample_report_qc cnv_report
```

Result summary:

```text
sample_report_qc_status: completed
PASS_REPORTABLE: 0
SAMPLE_QUALITY_REVIEW: 5
NO_CALL_RECOMMENDED: 0
```

All five 0615 samples are BAM-only in the current config, so fastp GC /
duplication / insert-size JSON is unavailable and BAM QC records
`FASTP_METRICS_MISSING` as library context. This does not produce no-call.

60 is the only sample with the expected high-wave/global instability pattern:

```text
sample: JZ26125845-60-60
status: SAMPLE_QUALITY_REVIEW
reasons:
  BAM_QC_REVIEW
  HIGH_MAD_LOG1P_RELATIVE_OUTLIER
  AUTOSOMAL_CN_OUTSIDE_1P7_2P3_REVIEW
  HIGH_MULTI_CHROMOSOME_REPORT_BURDEN_REVIEW
autosomal_cn_outside_1p7_2p3_fraction: 0.728760
report_event_count: 23
```

56/59/61/62 are not no-called. They remain review-level because this batch has
missing fastp library context and multi-chromosome report burden; this is a
sample-level warning only and does not change Branch A/B/S event status.
