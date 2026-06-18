# Branch A/B V2 Phase 1/2 Implementation Report

Date: 2026-06-18

## Scope

This report implements the first two phases from `docs/reports/branch_ab_v2_rnd_constraints_2026-06-18.md`.

- Phase 1: candidate evidence ledger.
- Phase 2: N0/N1/N2 negative-bank labeling.

Both phases are shadow-only. They do not change `WisecondorX predict`, Branch A candidate inclusion, Branch B final events, or `cnv_report` output.

## Phase 1 Outputs

New per-sample outputs:

```text
wisecondorx/cnv/postprocess/evidence_ledger/{sample}.candidate_evidence.tsv
wisecondorx/cnv/postprocess/evidence_ledger/{sample}.summary.json
```

Implementation:

```text
pgta/predict/branch_b/evidence_ledger.py
scripts/predict.py cnv_evidence_ledger
rule cnv_branch_b_evidence_ledger
```

The ledger writes one row per Branch A candidate and joins available Branch B evidence by `a_candidate_id`.

Required contract fields include:

```text
sample_id, candidate_id, chrom, start, end, state,
a_zscore, a_ratio, raw_amplitude, corrected_amplitude,
attenuation_ratio, same_direction_fraction, flank_contrast,
hard_region_fraction, sample_noise_mad,
refmap_status, calibration_null_status, evidence_missing_reason,
final_disposition
```

Disposition values are constrained to:

```text
CONFIRMED
MOSAIC_SUSPECT
REVIEW_REQUIRED
LIKELY_ARTIFACT
NO_CALL
```

Missing evidence policy:

- Missing Branch B event -> `branch_b_event_missing`.
- Missing refmap evidence -> `refmap_status=UNKNOWN`.
- Missing calibration-null evidence -> `calibration_null_status=UNKNOWN`.
- Missing evidence is never converted to clean support.

## Phase 2 Outputs

New aggregate outputs when configured:

```text
wisecondorx/cnv/postprocess/negative_bank/negative_bank_labels.tsv
wisecondorx/cnv/postprocess/negative_bank/negative_bank_summary.json
```

Implementation:

```text
pgta/predict/branch_b/negative_bank.py
scripts/predict.py negative_bank_labels
rule cnv_negative_bank_labels
```

Seed manifest:

```text
docs/reports/branch_ab_v2_negative_bank_seed_2026-06-18.tsv
```

Current seed labels:

| label | samples | use |
|---|---|---|
| N0 | none | Matched-negative empirical null only; no H7-H16 sample is promoted to N0. |
| N1 | H9, H10, H11, H12, H15, H16 | Presumed-negative candidates for shadow reference evaluation only. |
| N2 | H7, H8, H13, H14 | Hold-out / pending review. |

Only `N0` is eligible for matched-negative empirical null. `N1` can be used for shadow reference evaluation only, and `N2` is excluded from negative-bank modeling.

## Workflow Contract

The `cnv` target now includes Phase 1 shadow evidence outputs when Branch B is enabled.

The negative-bank rule is config-gated by:

```yaml
core:
  wisecondorx:
    cnv:
      negative_bank:
        samples_tsv: ...
        version: ...
```

The `cnv_report` rule intentionally does not consume the evidence ledger or negative-bank outputs. This preserves the current report logic until ablation is reviewed.

## Remote Validation

RED check:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 /biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest tests/unit/test_branch_b_evidence_ledger.py tests/unit/test_negative_bank.py -q'
```

Initial result:

```text
ModuleNotFoundError for pgta.predict.branch_b.evidence_ledger and pgta.predict.branch_b.negative_bank
```

GREEN check:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 /biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest tests/unit/test_branch_ab_phase12_workflow_contract.py tests/unit/test_branch_b_evidence_ledger.py tests/unit/test_negative_bank.py -q'
```

Result:

```text
9 passed in 0.70s
```

Regression check:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 /biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest tests/unit/test_branch_ab_phase12_workflow_contract.py tests/unit/test_branch_b_evidence_ledger.py tests/unit/test_negative_bank.py tests/unit/test_branch_b_calling.py tests/unit/test_branch_b_artifact_rules.py tests/unit/test_branch_b_correction.py -q'
```

Result:

```text
69 passed in 1.00s
```

Dry-run:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && /biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile /data/project/CNV/PGT-A/refactor_validation_20260419/config_predict_build_ref_v2_mask_only.yaml --cores 1 -n cnv'
```

Initial dry-run result:

```text
Job stats:
cnv_branch_b_evidence_ledger 8
cnv_negative_bank_labels 1
```

After executing the shadow jobs, the final dry-run result was:

```text
Nothing to be done (all requested files are present and up to date).
```

## Generated Remote Outputs

Command:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && /biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile /data/project/CNV/PGT-A/refactor_validation_20260419/config_predict_build_ref_v2_mask_only.yaml --cores 1 cnv'
```

Result:

```text
8 evidence ledger jobs completed.
1 negative-bank label job completed.
Complete log: .snakemake/log/2026-06-18T155912.600677.snakemake.log
```

Evidence ledger row counts:

| sample | A candidate rows | missing evidence rows | disposition summary |
|---|---:|---:|---|
| Y1 | 17 | 0 | LIKELY_ARTIFACT=16; REVIEW_REQUIRED=1 |
| Y2 | 22 | 0 | LIKELY_ARTIFACT=19; REVIEW_REQUIRED=3 |
| Y3 | 15 | 0 | LIKELY_ARTIFACT=13; REVIEW_REQUIRED=2 |
| Y4 | 25 | 0 | LIKELY_ARTIFACT=22; REVIEW_REQUIRED=3 |
| Y5 | 9 | 0 | LIKELY_ARTIFACT=6; REVIEW_REQUIRED=3 |
| Y6 | 15 | 0 | LIKELY_ARTIFACT=12; REVIEW_REQUIRED=3 |
| Y7 | 19 | 0 | LIKELY_ARTIFACT=17; REVIEW_REQUIRED=2 |
| Y8 | 21 | 0 | LIKELY_ARTIFACT=18; REVIEW_REQUIRED=3 |

Negative-bank output:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419/results_build_ref_v2_mask_only/wisecondorx/cnv/postprocess/negative_bank/negative_bank_labels.tsv
/data/project/CNV/PGT-A/refactor_validation_20260419/results_build_ref_v2_mask_only/wisecondorx/cnv/postprocess/negative_bank/negative_bank_summary.json
```

Label counts:

```text
N1: 6
N2: 4
N0: 0
matched_negative_eligible: 0
```

This is expected for the seed manifest. H7-H16 are not promoted to N0.

## Remaining Gates

- Full report promotion remains blocked until Y1-Y8, H1-H6, H7-H16, and the 2026-06-15 shadow-report set are evaluated with no known-positive recall regression.
- Matched-negative percentile remains Phase 3 and must use only N0 samples.
