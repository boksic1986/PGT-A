# Branch B V2 Burden Stratification

date: 2026-06-21
status: materialized_development_only
active_reference_id: h_r0_shadow_ref_20260619
branch_a_input: explicit_merge_gap_bp_2000000_overlay

## Scope

This loop adds Branch B V2 burden stratification on top of already materialized
Branch A gap2m candidates. It does not change Branch A, does not promote
`merge_gap_bp=2_000_000` to default, and does not promote Branch B V2 or
Branch S to final report logic.

Truth benchmark cohorts remain limited to Y1-Y8 and H1-H6. H7-H16 and
2026-06-15 are burden/context cohorts only.

## Contract

Branch B V2 remains Branch-A-anchored evidence stratification:

- no B-only event creation;
- no legacy/current-code Branch B decision fields in V2 decisions;
- no hard suppression except workflow/reference contract risk;
- unknown background is review context, not benign evidence;
- GC/RC, B-side signal context, length tier, and clean support are evidence
  labels, not standalone truth/falsity calls.

New V2 burden fields:

- `v2_burden_reduction_version`
- `v2_burden_reduction_tier`
- `v2_burden_reduction_action`
- `v2_burden_reduction_reason`
- `v2_burden_evidence_tags`

The benchmark summary now stratifies by:

- `v2_filter_action`
- `v2_disposition`
- `v2_length_tier`
- `v2_clean_support_label`
- `v2_b_signal_context_label`
- `v2_gc_rc_context_label`
- `v2_burden_reduction_tier`
- `v2_burden_reduction_action`

## CNVseq / CNVpro Borrowed Concepts

- `[CNVpro-inspired]` length tiers:
  `large_ge10mb`, `broad_review_ge4mb`, `reportable_candidate_ge2mb`,
  `review_only_ge1mb`, and focal/high-risk below 1 Mb. These tiers affect
  report/review routing only and do not decide truth.
- `[CNVpro-confirmed]` acrocentric qter context for chr13/14/15/21/22. Current
  implementation tags residual/region-risk context only; it does not hard mask
  acrocentric short arms or replace WisecondorX/CBS.
- `[CNVseq-asset]` sex homology/PAR and mask/mappability context. These are
  annotation and Branch S review context only; the hg19 FASTA and WisecondorX
  main chain are not replaced.
- `[CNVpro-like]` GC/RC context. It can indicate attenuation/amplification but
  cannot reverse or hard-suppress a Branch A candidate.
- `[Not used]` CNVpro `CNVcalling.R` / cghFLasso are not used as primary
  segmentation or CNV calling replacements.

## Materialized Results

Remote path:
`/data/project/CNV/PGT-A/refactor_validation_20260419`.

| cohort | candidates | truth events | truth preserved | FN | hard-suppressed truth | filter-suppressed candidates | report burden | review burden | background review | technical risk | Branch S review |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Y1-Y8 | 97 | 10 | 10 | 0 | 0 | 0 | 0 | 84 | 84 | 0 | 13 |
| H1-H16 | 105 | 10 | 10 | 0 | 0 | 0 | 0 | 63 | 63 | 0 | 42 |
| 2026-06-15 | 165 | 0 | 0 | 0 | 0 | 0 | 0 | 151 | 151 | 0 | 14 |

Burden tiers:

| cohort | background_unknown_review | branch_s_review |
|---|---:|---:|
| Y1-Y8 | 84 | 13 |
| H1-H16 | 63 | 42 |
| 2026-06-15 | 151 | 14 |

Length tiers:

| cohort | large >=10Mb | broad >=4Mb | reportable >=2Mb | review-only >=1Mb |
|---|---:|---:|---:|---:|
| Y1-Y8 | 85 | 3 | 5 | 4 |
| H1-H16 | 96 | 4 | 2 | 3 |
| 2026-06-15 | 137 | 11 | 9 | 8 |

B-side signal context:

| cohort | supported A direction | weak B signal | discordant with A direction |
|---|---:|---:|---:|
| Y1-Y8 | 66 | 20 | 11 |
| H1-H16 | 68 | 26 | 11 |
| 2026-06-15 | 97 | 40 | 28 |

GC/RC context:

| cohort | stable | attenuated | attenuated severe | amplified |
|---|---:|---:|---:|---:|
| Y1-Y8 | 24 | 39 | 30 | 4 |
| H1-H16 | 30 | 30 | 16 | 29 |
| 2026-06-15 | 38 | 37 | 35 | 55 |

## Truth-Safety Checks

Y1-Y8 and H1-H6 locked truth remain preserved 10/10 with FN=0 and
hard-suppressed truth=0.

H6 chr21 gain is retained:

- top candidate: `H6.A0003`
- `top_a_abs_zscore=7.113507`
- `top_v2_disposition=background_unknown_review`
- `top_v2_filter_action=keep_background_unknown_review`
- `top_v2_burden_reduction_tier=background_unknown_review`
- `top_v2_burden_reduction_action=stratify_background_unknown_review`

This confirms the burden stratification did not drop the weak-positive /
mosaic-sensitive truth example.

## Interpretation

This loop improves auditability and separates evidence types for later report
burden control. It does not yet reduce total review burden, because current
materialized tiers are still only:

- `background_unknown_review`
- `branch_s_review`

Therefore Branch B V2 remains development-only. The next step is to evaluate
whether any burden tier can move from broad review to a more limited
report/review display without hard-suppressing locked truth candidates.

## Remote Validation

Unit test:

```text
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
tests/unit/test_branch_b_v2_classifier.py \
tests/unit/test_branch_b_v2_benchmark.py -q

30 passed in 0.88s
```

Dry-run:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake \
-s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile \
--configfile <active gap2m config> --cores 1 -n \
branch_b_v2_benchmark branch_s_review cnv_report
```

Passed for:

- `config_predict_y_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`
- `config_predict_h_20260608_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`
- `config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`

Forced materialization refreshed `cnv_branch_b_v2_classifier` and
`cnv_branch_b_v2_benchmark` for all three configs.

## Next Gate

Proceed to burden-display and report-contract integration, not hard filtering.
Any future suppression or demotion rule must pass locked Y1-Y8/H1-H6 truth
preservation, explicitly check H6 chr21, and keep H7-H16/0615 as context only
unless locked truth labels are added.
