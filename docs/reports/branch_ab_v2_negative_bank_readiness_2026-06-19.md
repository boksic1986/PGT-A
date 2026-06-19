# Branch A/B V2 Negative-Bank Readiness Audit

Date: 2026-06-19

Remote project path:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419
```

## Purpose

This audit tightens the Phase 2 to Phase 3 contract for Branch B V2.

The negative-bank output now records whether matched-negative percentile evidence is actually ready. This does not change WisecondorX predict, Branch A candidate generation, Branch B final events, Branch S, or the final report.

## Code Contract

`pgta/predict/branch_b/negative_bank.py` now writes these additional summary fields:

```text
matched_negative_ready
matched_negative_blocking_reason
n0_sample_ids
n1_shadow_reference_candidate_count
n2_holdout_count
```

Interpretation:

- `matched_negative_ready=true` only when at least one `N0` / `matched_negative_eligible=1` sample exists.
- `matched_negative_ready=false` with `matched_negative_blocking_reason=no_n0_locked_clean_negative_samples` means Phase 3 must fall back to `UNKNOWN_BACKGROUND`.
- `N1` remains shadow-reference candidate context only.
- `N2` remains hold-out / review context.
- `N1` and `N2` are not eligible empirical-null samples.

## Current Remote Result

Output:

```text
results_build_ref_v2_mask_only/wisecondorx/cnv/postprocess/negative_bank/negative_bank_summary.json
```

Current summary:

```json
{
  "version": "branch_ab_v2_seed_2026-06-18",
  "sample_count": 10,
  "label_counts": {
    "N1": 6,
    "N2": 4
  },
  "matched_negative_eligible_count": 0,
  "matched_negative_ready": false,
  "matched_negative_blocking_reason": "no_n0_locked_clean_negative_samples",
  "n0_sample_ids": [],
  "n1_shadow_reference_candidate_count": 6,
  "n2_holdout_count": 4,
  "n0_only_for_empirical_null": true
}
```

Current H7-H16 interpretation remains unchanged:

- `N1`: H9, H10, H11, H12, H15, H16
- `N2`: H7, H8, H13, H14
- `N0`: none

## Downstream State

After refreshing `config_predict_build_ref_v2_mask_only.yaml`:

```text
matched-negative files: 8
background_status_counts: UNKNOWN_BACKGROUND=143
matched_negative_action: REVIEW_NO_CALL=143

V2 classifier files: 8
class_counts: REVIEW_REQUIRED=143
action_counts: SHADOW_REVIEW_ONLY=143
```

This confirms the current Phase 3 and Phase 4 outputs remain conservative and report-safe.

## Remote Validation

RED test:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 /biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest tests/unit/test_negative_bank.py -q'
```

Result before implementation:

```text
2 failed, 2 passed
KeyError: 'matched_negative_ready'
```

GREEN target test:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 /biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest tests/unit/test_negative_bank.py -q'
```

Result:

```text
4 passed in 0.47s
```

Branch A/B/S regression:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 /biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest tests/unit/test_branch_ab_phase12_workflow_contract.py tests/unit/test_branch_b_evidence_ledger.py tests/unit/test_negative_bank.py tests/unit/test_branch_b_matched_negative.py tests/unit/test_branch_b_v2_classifier.py tests/unit/test_branch_s_shadow.py -q'
```

Result:

```text
25 passed in 0.74s
```

Snakemake dry-run:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && /biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile /data/project/CNV/PGT-A/refactor_validation_20260419/config_predict_build_ref_v2_mask_only.yaml --cores 1 --forcerun cnv_negative_bank_labels -n cnv'
```

Result:

```text
18 jobs in dry-run:
cnv_negative_bank_labels=1
cnv_branch_b_matched_negative=8
cnv_branch_b_v2_classifier=8
cnv=1
```

Snakemake materialization:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && /biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile /data/project/CNV/PGT-A/refactor_validation_20260419/config_predict_build_ref_v2_mask_only.yaml --cores 4 --forcerun cnv_negative_bank_labels cnv'
```

Result:

```text
20 of 20 steps done
```

Final dry-run:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && /biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile /data/project/CNV/PGT-A/refactor_validation_20260419/config_predict_build_ref_v2_mask_only.yaml --cores 1 -n cnv'
```

Result:

```text
Nothing to be done
```

## Decision

Phase 2 is now machine-readable:

- H7-H16 are not a usable N0 empirical-null set.
- Phase 3 matched-negative percentile is not ready for FP reduction.
- Phase 4 V2 classifier remains shadow review only.

Next minimum action remains data/cohort work: define locked N0 samples independent of H7-H16 review candidates, then rerun WisecondorX predict, Branch A candidates, Branch B, evaluation, benchmark, and report.
