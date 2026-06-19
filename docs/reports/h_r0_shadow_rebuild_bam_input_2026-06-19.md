# H R0 Shadow Reference BAM Input Fix 2026-06-19

## Objective

Fix the H7-H16 shadow reference rebuild so it reuses existing BAM/BAI inputs
instead of regenerating BAM files.

## Root Cause

The workflow previously treated `core.project_path/mapping/{sample}.sorted.bam`
as the only BAM input for QC, reference tuning/build, and WisecondorX predict.

When the H R0 shadow reference used a new output directory:

`/data/project/CNV/PGT-A/refactor_validation_20260419/results_h_r0_shadow_ref_20260619`

Snakemake expanded missing outputs under:

`results_h_r0_shadow_ref_20260619/mapping/{sample}.sorted.bam`

and therefore scheduled `fastp_bwa`.

This was not biologically or workflow-contractually required because the mapping
reference and read-level preprocessing contract were unchanged.

## Fix

The workflow now supports per-sample existing artifact overrides:

- `samples.<sample_id>.bam`
- `samples.<sample_id>.bai`
- `samples.<sample_id>.fastp_json`

Rules that consume BAMs now resolve BAM input through this single contract:

- baseline QC
- baseline samtools diagnostics
- reference prefilter
- reference tuning
- reference build
- WisecondorX predict convert
- WisecondorX gender routing

When a sample has external `bam`, missing `fastp_json` no longer pulls
`fastp_bwa` back into the DAG. Existing fastp JSON can still be supplied for
library-level report metrics.

## H R0 Shadow Input Set

Config:

`config_reference_h_r0_shadow_20260619.yaml`

Reference sample count:

`42`

H R0 samples added to the XY reference group:

`H9,H10,H11,H12,H15,H16`

External BAM sources:

- Previous mask-only reference run:
  `/data/project/CNV/PGT-A/refactor_validation_20260419/results_build_ref_v2_mask_only/mapping`
- H batch run:
  `/data/project/CNV/PGT-A/refactor_validation_20260419/results_h_20260608_mask_only/mapping`

H batch fastp JSON source:

`/data/project/CNV/PGT-A/refactor_validation_20260419/results_h_20260608_mask_only/fastp`

## Remote Validation

Remote executable:

`/biosoftware/miniconda/envs/snakemake_env/bin/snakemake`

Dry-run command:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake \
  -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile \
  --configfile /data/project/CNV/PGT-A/refactor_validation_20260419/config_reference_h_r0_shadow_20260619.yaml \
  --directory /data/project/CNV/PGT-A/refactor_validation_20260419 \
  --cores 1 \
  -n \
  reference_qc reference
```

Result:

- Dry-run succeeded.
- DAG contained 18 jobs.
- `fastp_bwa` was absent from the DAG.
- `baseline_bam_uniformity_qc` consumed existing BAM paths from the old
  reference run and the H 2026-06-08 run.

External artifact check:

- `ref_sample_count=42`
- `external_bam_count=42`
- `missing_bam_count=0`
- `missing_bai_count=0`

## Current Running Task

Remote command:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake \
  -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile \
  --configfile /data/project/CNV/PGT-A/refactor_validation_20260419/config_reference_h_r0_shadow_20260619.yaml \
  --directory /data/project/CNV/PGT-A/refactor_validation_20260419 \
  --cores 60 \
  --rerun-incomplete \
  --latency-wait 60 \
  --printshellcmds \
  reference_qc reference
```

PID:

`61980`

Log:

`/data/project/CNV/PGT-A/refactor_validation_20260419/results_h_r0_shadow_ref_20260619/logs/driver/h_r0_shadow_reference_20260619.snakemake.log`

PID file:

`/data/project/CNV/PGT-A/refactor_validation_20260419/results_h_r0_shadow_ref_20260619/logs/driver/h_r0_shadow_reference_20260619.pid`

Expected primary outputs:

- `/data/project/CNV/PGT-A/refactor_validation_20260419/results_h_r0_shadow_ref_20260619/reference/cohorts/pass_warn_h_r0/XX/result/ref_XX_best.npz`
- `/data/project/CNV/PGT-A/refactor_validation_20260419/results_h_r0_shadow_ref_20260619/reference/cohorts/pass_warn_h_r0/XY/result/ref_XY_best.npz`
- `/data/project/CNV/PGT-A/refactor_validation_20260419/results_h_r0_shadow_ref_20260619/reference/cohorts/pass_warn_h_r0/gender/result/ref_gender_best.npz`
- `/data/project/CNV/PGT-A/refactor_validation_20260419/results_h_r0_shadow_ref_20260619/reference/gender/common_best_binsize.txt`

## Next Step

After the shadow reference finishes, rerun against the new reference:

1. WisecondorX predict.
2. Branch A candidates.
3. Branch B evidence/classification.
4. Negative-bank labeling tied to the new reference ID.
5. Evaluation.
6. Benchmark.
7. Report.

The old `N0=0` state must not be reused as the post-rebuild Phase 3 result.
