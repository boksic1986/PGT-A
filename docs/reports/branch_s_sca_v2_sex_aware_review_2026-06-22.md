---
status: development_review_only
decision_use: branch_s_sca_v2_current_evidence
created_at: 2026-06-22
active_reference_id: h_r0_shadow_ref_20260619
branch_a_overlay: merge_gap_bp_2000000
lowres_used: false
---

# Branch S / SCA V2 Sex-Aware Review

## Purpose

This loop fixes a Branch S failure mode: a strong Branch A chrX candidate in a
normal XY sample could be interpreted as `X_GAIN` / `SCA_REVIEW_STRONG` even
when X non-PAR bin-level evidence was neutral.

The change keeps Branch A as the sex-chromosome candidate anchor, but prevents
Branch A alone from deciding SCA strength. Branch S remains
`review_development_only`; it is not final SCA calling and does not replace sex
calling.

The 2Mb/3Mb low-resolution Snakemake task was not stopped, but this report does
not read or use any low-resolution output.

## Method Change

Branch S now evaluates SCA evidence in this order:

1. Use sex call as context, not as a replacement for SCA evidence.
2. Use Branch A chrX/chrY candidates as anchors only.
3. Require X/Y non-PAR median or robust bin evidence to corroborate the Branch A
   direction before promoting to strong SCA review.
4. Treat PAR / XY-homology as context only, not as a global hard mask.
5. Route output to report-layer classes:
   - `sca_report_review_event`
   - `sca_internal_review_event`
   - `sca_filtered_or_sex_consistent_event`
   - `sca_no_call`

Important exception:

`XX + X_LOSS` with very strong Branch A X-loss support remains visible as SCA
review even if the postprocessed non-PAR median is neutral. This protects the
current locked SCA truth pattern in Y3/H5/H6, where Branch A carries the SCA
signal and the downstream bin median can be near neutral.

## Remote Materialized Results

All materialization and validation below was run in:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419
```

with:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/python
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake
```

Active result paths:

```text
results_y_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m/branch_s
results_h_20260608_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m/branch_s
results_g1_g8_h_r0_shadow_ref_20260621/wisecondorx/cnv/postprocess_gap2m/branch_s
results_20260615_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m/branch_s
```

### Cohort Summary

This table is the Branch S / SCA section after the sex-aware correction. It is
not the autosomal CNV report-event count.

| cohort | samples | SCA review strong | SCA review weak | SCA no-call | report review | internal review | sex-consistent/audit | no-call |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Y1-Y8 | 8 | 1 | 0 | 7 | 1 | 0 | 5 | 2 |
| H1-H16 | 16 | 2 | 1 | 13 | 2 | 1 | 12 | 1 |
| G1-G8 | 8 | 0 | 2 | 6 | 0 | 2 | 3 | 3 |
| 2026-06-15 | 5 | 0 | 2 | 3 | 0 | 2 | 3 | 0 |

### CNV Report-Layer Totals

These counts come from `wisecondorx/cnv/report/cnv_summary.tsv`. They describe
the full Branch B V2 report-layer burden. `branch_s_routed_events` is the number
of sex-chromosome candidates routed out of the autosomal CNV table and into
Branch S context; it is not the number of SCA positives.

| cohort | samples | report events | internal review events | filtered audit-only events | branch-s routed events | sample burden flags |
|---|---:|---:|---:|---:|---:|---:|
| Y1-Y8 | 8 | 21 | 50 | 13 | 13 | 2 |
| H1-H16 | 16 | 6 | 49 | 8 | 42 | 1 |
| G1-G8 | 8 | 15 | 40 | 7 | 13 | 1 |
| 2026-06-15 | 5 | 52 | 76 | 23 | 14 | 5 |

The SCA correction primarily reduces false strong SCA interpretation inside the
Branch S section. It does not by itself reduce autosomal `report_events` or
`internal_review_events`.

### SCA Sample-Level Results

| cohort | sample | sex | Branch A X | Branch A Y | SCA state | SCA tier | report layer | reason |
|---|---|---|---|---|---|---|---|---|
| Y1-Y8 | Y1 | XY | present | absent | none_detected | SCA_NO_CALL | sca_filtered_or_sex_consistent_event | branch_a_only_uncorroborated_by_nonpar_median |
| Y1-Y8 | Y2 | XX | present | absent | none_detected | SCA_NO_CALL | sca_filtered_or_sex_consistent_event | branch_a_only_uncorroborated_by_nonpar_median |
| Y1-Y8 | Y3 | XX | present | absent | X_LOSS | SCA_REVIEW_STRONG | sca_report_review_event | sca_review_strong_with_nonpar_corroboration |
| Y1-Y8 | Y4 | XX | present | absent | none_detected | SCA_NO_CALL | sca_filtered_or_sex_consistent_event | branch_a_only_uncorroborated_by_nonpar_median |
| Y1-Y8 | Y5 | XX | absent | absent | none_detected | SCA_NO_CALL | sca_no_call | insufficient_sca_evidence |
| Y1-Y8 | Y6 | XY | present | absent | none_detected | SCA_NO_CALL | sca_filtered_or_sex_consistent_event | branch_a_only_uncorroborated_by_nonpar_median |
| Y1-Y8 | Y7 | XX | absent | absent | none_detected | SCA_NO_CALL | sca_no_call | insufficient_sca_evidence |
| Y1-Y8 | Y8 | XY | present | absent | none_detected | SCA_NO_CALL | sca_filtered_or_sex_consistent_event | branch_a_only_uncorroborated_by_nonpar_median |
| H1-H16 | H1 | XY | present | absent | none_detected | SCA_NO_CALL | sca_filtered_or_sex_consistent_event | branch_a_only_uncorroborated_by_nonpar_median |
| H1-H16 | H2 | XY | present | absent | none_detected | SCA_NO_CALL | sca_filtered_or_sex_consistent_event | branch_a_only_uncorroborated_by_nonpar_median |
| H1-H16 | H3 | XX | absent | absent | none_detected | SCA_NO_CALL | sca_no_call | insufficient_sca_evidence |
| H1-H16 | H4 | XX | present | absent | X_LOSS | SCA_REVIEW_WEAK | sca_internal_review_event | sca_review_weak_with_nonpar_corroboration |
| H1-H16 | H5 | XX | present | absent | X_LOSS | SCA_REVIEW_STRONG | sca_report_review_event | sca_review_strong_with_nonpar_corroboration |
| H1-H16 | H6 | XX | present | absent | X_LOSS | SCA_REVIEW_STRONG | sca_report_review_event | sca_review_strong_with_nonpar_corroboration |
| H1-H16 | H7 | XY | present | absent | none_detected | SCA_NO_CALL | sca_filtered_or_sex_consistent_event | branch_a_only_uncorroborated_by_nonpar_median |
| H1-H16 | H8 | XY | present | absent | none_detected | SCA_NO_CALL | sca_filtered_or_sex_consistent_event | branch_a_only_uncorroborated_by_nonpar_median |
| H1-H16 | H9 | XY | present | absent | none_detected | SCA_NO_CALL | sca_filtered_or_sex_consistent_event | branch_a_only_uncorroborated_by_nonpar_median |
| H1-H16 | H10 | XY | present | absent | none_detected | SCA_NO_CALL | sca_filtered_or_sex_consistent_event | branch_a_only_uncorroborated_by_nonpar_median |
| H1-H16 | H11 | XY | present | absent | none_detected | SCA_NO_CALL | sca_filtered_or_sex_consistent_event | branch_a_only_uncorroborated_by_nonpar_median |
| H1-H16 | H12 | XY | present | absent | none_detected | SCA_NO_CALL | sca_filtered_or_sex_consistent_event | branch_a_only_uncorroborated_by_nonpar_median |
| H1-H16 | H13 | XY | present | absent | none_detected | SCA_NO_CALL | sca_filtered_or_sex_consistent_event | branch_a_only_uncorroborated_by_nonpar_median |
| H1-H16 | H14 | XY | present | absent | none_detected | SCA_NO_CALL | sca_filtered_or_sex_consistent_event | branch_a_only_uncorroborated_by_nonpar_median |
| H1-H16 | H15 | XY | present | absent | none_detected | SCA_NO_CALL | sca_filtered_or_sex_consistent_event | branch_a_only_uncorroborated_by_nonpar_median |
| H1-H16 | H16 | XY | present | absent | none_detected | SCA_NO_CALL | sca_filtered_or_sex_consistent_event | branch_a_only_uncorroborated_by_nonpar_median |
| G1-G8 | G1 | XX | absent | absent | none_detected | SCA_NO_CALL | sca_no_call | insufficient_sca_evidence |
| G1-G8 | G2 | XY | present | absent | none_detected | SCA_NO_CALL | sca_filtered_or_sex_consistent_event | branch_a_only_uncorroborated_by_nonpar_median |
| G1-G8 | G3 | XX | present | absent | X_LOSS | SCA_REVIEW_WEAK | sca_internal_review_event | sca_review_weak_with_nonpar_corroboration |
| G1-G8 | G4 | XX | absent | absent | none_detected | SCA_NO_CALL | sca_no_call | insufficient_sca_evidence |
| G1-G8 | G5 | XX | present | absent | X_LOSS | SCA_REVIEW_WEAK | sca_internal_review_event | sca_review_weak_with_nonpar_corroboration |
| G1-G8 | G6 | XY | present | absent | none_detected | SCA_NO_CALL | sca_filtered_or_sex_consistent_event | branch_a_only_uncorroborated_by_nonpar_median |
| G1-G8 | G7 | XX | absent | absent | none_detected | SCA_NO_CALL | sca_no_call | insufficient_sca_evidence |
| G1-G8 | G8 | XY | present | absent | none_detected | SCA_NO_CALL | sca_filtered_or_sex_consistent_event | branch_a_only_uncorroborated_by_nonpar_median |
| 2026-06-15 | JZ26125843-56-56 | XX | present | absent | X_LOSS | SCA_REVIEW_WEAK | sca_internal_review_event | sca_review_weak_with_nonpar_corroboration |
| 2026-06-15 | JZ26125844-59-59 | XX | present | absent | X_LOSS | SCA_REVIEW_WEAK | sca_internal_review_event | sca_review_weak_with_nonpar_corroboration |
| 2026-06-15 | JZ26125845-60-60 | XY | present | absent | none_detected | SCA_NO_CALL | sca_filtered_or_sex_consistent_event | branch_a_only_uncorroborated_by_nonpar_median |
| 2026-06-15 | JZ26125846-61-61 | XY | present | absent | none_detected | SCA_NO_CALL | sca_filtered_or_sex_consistent_event | branch_a_only_uncorroborated_by_nonpar_median |
| 2026-06-15 | JZ26125847-62-62 | XY | present | absent | none_detected | SCA_NO_CALL | sca_filtered_or_sex_consistent_event | branch_a_only_uncorroborated_by_nonpar_median |

### Locked SCA Truth Visibility

| sample | expected evidence | sex_call | state after fix | confidence | report layer | status |
|---|---|---|---|---|---|---|
| Y3 | chrX loss / XO-like review truth | XX | X_LOSS | SCA_REVIEW_STRONG | sca_report_review_event | visible |
| H5 | chrX loss truth | XX | X_LOSS | SCA_REVIEW_STRONG | sca_report_review_event | visible |
| H6 | chrX loss truth | XX | X_LOSS | SCA_REVIEW_STRONG | sca_report_review_event | visible |

### H7-H16 XY Context

All H7-H16 normal XY context samples now have:

```text
sca_candidate_state=none_detected
sca_confidence_tier=SCA_NO_CALL
sca_report_layer_class=sca_filtered_or_sex_consistent_event
sca_report_layer_reason=branch_a_only_uncorroborated_by_nonpar_median
```

This means their chrX Branch A support is retained in the audit context, but no
longer becomes `X_GAIN` / `SCA_REVIEW_STRONG`.

## Tests And Dry-Runs

Remote unit tests:

```text
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_s_shadow.py \
  tests/unit/test_cnv_report.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py \
  tests/unit/test_current_context_index.py -q
```

Result:

```text
38 passed in 1.10s
```

Remote dry-run target for all active gap2m configs:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile <config> --cores 1 -n branch_s_review cnv_report
```

Configs checked:

```text
config_predict_y_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml
config_predict_h_20260608_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml
config_predict_g1_g8_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml
config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml
```

Result:

```text
RC=0 for all four configs; all requested files were present and up to date.
```

## Interpretation

This is a Branch S report-layer correction, not production SCA promotion.

Current evidence supports:

- normal XY chrX Branch A-only signals are no longer promoted to strong SCA;
- Y3/H5/H6 locked X-loss evidence remains visible;
- Branch S output is now easier to separate from autosomal CNV report events.

Current evidence does not support:

- final SCA PASS/FAIL calls;
- production SCA confirmation;
- hard suppression of sex-chromosome evidence without a locked SCA truth panel;
- use of 2Mb/3Mb low-resolution evidence in this SCA decision.

## Remaining Risks

- Locked SCA truth is still narrow: current positive evidence is mostly X-loss.
- X gain, XXY, XYY, Y-loss, mosaic SCA fraction series, and cross-batch clean
  XX/XY negatives remain required before final SCA promotion.
- The `XX + X_LOSS` protection is intentionally conservative for current truth
  preservation and must be revisited when a broader SCA truth panel exists.
- 2026-06-15 remains burden/context only because it has no locked truth labels.
