---
status: active_rnd_report
decision_use: branch_b_v2_method_reset
date: 2026-06-21
active_reference_id: h_r0_shadow_ref_20260619
branch_a_overlay: merge_gap_bp_2000000
---

# Branch B V2 Method Reset And Threshold Inventory

## Scope

This report resets the Branch B V2 interpretation after the background-context
and direction-support review-label loops.

Branch A remains the WisecondorX/CBS-derived primary discovery layer. Branch B
V2 does not create B-only report events and does not use GC/RC corrected signal
direction to contradict Branch A. Its current role is candidate-level evidence
stratification and report/review support for Branch A candidates.

This is still a shadow/development contract. It does not promote Branch B V2,
Branch S, the explicit `merge_gap_bp=2_000_000` overlay, or the shadow
reference to production.

## Current Anchor

Active Branch A candidate anchor fields:

```text
sample_id
chrom
start
end
state
a_zscore
a_abs_zscore
a_support_level
```

Active Branch A overlay:

```text
merge_gap_bp=2_000_000
strong_z=10.0
postprocess.output_dir=wisecondorx/cnv/postprocess_gap2m
```

`strong_z=10` is a strong/sensitive label only. It is not an inclusion or
exclusion hard filter. Signals below 10 remain eligible for weak-positive or
mosaic-sensitive review.

## Legacy Threshold Inventory

The following thresholds exist in the current workflow/codebase but must not be
silently carried into Branch B V2 as final hard filters.

| threshold | current default | source | current V2 interpretation |
|---|---:|---|---|
| `CNV_CALLING_MIN_BINS` | 5 | `rules/predict_layout.smk` / Branch B calling | Legacy/current-code Branch B segmentation/calling support, not V2 final decision evidence. |
| `CNV_CALLING_MIN_EVENT_BINS` | 3 | `rules/predict_layout.smk` / Branch B calling | Legacy/current-code event assembly threshold, not V2 hard inclusion. |
| `CNV_CALLING_MIN_EVENT_Z` | 1.5 | `rules/predict_layout.smk` / Branch B calling | Legacy/current-code Branch B signal threshold, not V2 hard inclusion. |
| `CNV_ARTIFACT_MIN_BINS` | 3 | `rules/predict_layout.smk` / artifact rules | Legacy artifact-rule threshold; excluded from V2 performance evidence. |
| `CNV_ARTIFACT_MIN_ABS_Z` | 2.0 | `rules/predict_layout.smk` / artifact rules | Legacy artifact-rule calibrated-z threshold; excluded from V2 final disposition. |
| `CNV_ARTIFACT_HIGH_CONF_Z` | 4.0 | `rules/predict_layout.smk` / artifact rules | Legacy artifact-rule confidence threshold; may be historical context only. |
| `CNV_ARTIFACT_A_BRANCH_REVIEW_MIN_ABS_Z` | 15.0 | `rules/predict_layout.smk` / artifact rules | Legacy A-protection/review threshold; not a V2 pass/fail threshold. |
| `CNV_ARTIFACT_A_BRANCH_SENSITIVE_REVIEW_MIN_ABS_Z` | 7.0 | `rules/predict_layout.smk` / artifact rules | Legacy sensitive-review threshold; overlaps with H6 chr21 range and cannot suppress truth. |
| `CNV_ARTIFACT_A_BRANCH_SENSITIVE_REVIEW_MIN_BINS` | 10 | `rules/predict_layout.smk` / artifact rules | Legacy review-threshold context only. |
| `CNV_ARTIFACT_BRANCH_B_DIRECTION_MIN_ABS_Z` | 0.25 | `rules/predict_layout.smk` / artifact rules | Legacy direction-discordance trigger; replaced in V2 by non-final signal-discordance review label. |
| `CNV_ARTIFACT_CNVSEQ_REPORTABLE_MIN_BP` | 2,000,000 | `rules/predict_layout.smk` / CNVpro-like tier | Review/reportability tier only; not a truth/falsity decision. |
| `CNV_ARTIFACT_CNVSEQ_REVIEW_MIN_BP` | 1,000,000 | `rules/predict_layout.smk` / CNVpro-like tier | Review tier only; not a deletion rule. |
| `CNV_ARTIFACT_CNVSEQ_LARGE_EVENT_MIN_BP` | 10,000,000 | `rules/predict_layout.smk` / CNVpro-like tier | Large-event tier only; not a final pass rule. |
| `CNV_ARTIFACT_CNVSEQ_BOUNDARY_MAX_ABS_Z` | 4.0 | `rules/predict_layout.smk` / artifact rules | Legacy boundary-artifact threshold; may inform risk explanation only after truth-safe ablation. |

Legacy/current-code fields that remain excluded from V2 decisions and metrics:

```text
final_disposition
branch_b_keep_event
branch_b_report_class
branch_b_artifact_status
legacy artifact labels
old N0=0 / N1-only matched-negative promotion
```

## V2 Evidence Layer Contract

Branch B V2 now emits explicit evidence/context fields for each Branch A
candidate:

```text
v2_signal_strength_tier
v2_length_tier
v2_clean_support_label
v2_gc_rc_context_label
v2_background_context_label
v2_b_signal_context_label
v2_b_signal_context_reason
v2_disposition
```

### A Signal Strength

```text
A_STRONG_Z_GE_10
A_SENSITIVE_Z_5_TO_10
A_WEAK_Z_LT_5
A_Z_MISSING
```

The H6 chr21 truth event is expected to remain in
`A_SENSITIVE_Z_5_TO_10` and must not be suppressed by Branch B V2.

### Length Tier

```text
large_ge10mb
broad_review_ge4mb
reportable_candidate_ge2mb
review_only_ge1mb
focal_high_risk_lt1mb
unknown_length
```

These are report/review tiers. They do not decide whether the candidate is true
or false.

### GC/RC Context

GC/RC correction remains auxiliary evidence:

```text
GC_RC_CONTEXT_UNKNOWN
GC_RC_ATTENUATED_SEVERE
GC_RC_ATTENUATED
GC_RC_STABLE
GC_RC_AMPLIFIED
```

It must not define the candidate gain/loss direction and must not hard-filter a
Branch A positive candidate.

### B-Side Signal Context

The old direction-conflict semantics are replaced with a review-context label:

```text
B_SIGNAL_SUPPORTED_A_DIRECTION
A_ANCHORED_WEAK_B_SIGNAL
B_SIGNAL_DISCORDANT_WITH_A_DIRECTION
NO_POSITIVE_A_SIGNAL
```

`B_SIGNAL_DISCORDANT_WITH_A_DIRECTION` means the downstream B-side auxiliary
signal does not support the Branch A candidate direction. It does not mean
Branch B called an opposite event, and it is not a hard artifact call.

### V2 Disposition

First-round V2 dispositions are intentionally conservative:

```text
report_candidate
review_candidate
technical_risk_review
background_unknown_review
sca_branch_s_review
```

`UNKNOWN_BACKGROUND` is never benign. It is routed to review context unless and
until a formal background/null contract is available.

## Validation Contract

Performance validation remains limited to locked truth cohorts:

- Y1-Y8 truth events.
- H1-H6 truth events.

H7-H16 and 2026-06-15 are burden/context cohorts only. They must not be used to
compute TP/FN/FP unless locked truth labels are added.

Required truth-safety checks:

```text
Y1-Y8 truth preserved 10/10
H1-H6 truth preserved 10/10
H6 chr21 retained
hard-suppressed truth=0
```

## Remote Verification

Unit test status before materialization:

```text
tests/unit/test_branch_b_v2_classifier.py
tests/unit/test_branch_b_v2_benchmark.py

22 passed in 0.77s
```

Broader remote unit test status after syncing this loop:

```text
tests/unit/test_branch_b_v2_classifier.py
tests/unit/test_branch_b_v2_benchmark.py
tests/unit/test_branch_ab_phase12_workflow_contract.py
tests/unit/test_cnv_report.py

45 passed in 1.26s
```

Snakemake dry-run status:

```text
branch_b_v2_benchmark branch_s_review cnv_report
```

Dry-runs succeeded for:

- `config_predict_y_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`;
- `config_predict_h_20260608_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`;
- `config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml`.

Forced remote materialization used:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake \
  -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile \
  --configfile <config> \
  --directory /data/project/CNV/PGT-A/refactor_validation_20260419 \
  --cores 8 \
  --forcerun cnv_branch_b_v2_classifier cnv_branch_b_v2_benchmark \
  branch_b_v2_benchmark
```

All three forced materializations completed.

### Benchmark Summary

| cohort | samples | candidates | truth events | truth preserved | FN | hard-suppressed truth | status |
|---|---:|---:|---:|---:|---:|---:|---|
| Y1-Y8 | 8 | 97 | 10 | 10 | 0 | 0 | ready |
| H1-H16 | 16 | 105 | 10 | 10 | 0 | 0 | ready |
| 2026-06-15 | 5 | 165 | 0 | 0 | 0 | 0 | skipped_no_truth |

0615 remains burden/context only. It cannot produce TP/FN/FP performance
claims without locked truth labels.

### Truth Top-Candidate V2 Disposition

Locked Y1-Y8 truth top candidates:

| truth | top candidate | A abs z | V2 disposition | length tier | B signal context | GC/RC context |
|---|---|---:|---|---|---|---|
| Y1 chr21 loss | Y1.A0003 | 112.6040 | background_unknown_review | large_ge10mb | B_SIGNAL_SUPPORTED_A_DIRECTION | GC_RC_ATTENUATED |
| Y2 chr14 gain | Y2.A0001 | 94.7847 | background_unknown_review | large_ge10mb | A_ANCHORED_WEAK_B_SIGNAL | GC_RC_ATTENUATED |
| Y2 chr21 gain | Y2.A0002 | 67.3226 | background_unknown_review | large_ge10mb | B_SIGNAL_SUPPORTED_A_DIRECTION | GC_RC_ATTENUATED_SEVERE |
| Y3 chrX loss | Y3.A0001 | 59.5924 | sca_branch_s_review | large_ge10mb | B_SIGNAL_DISCORDANT_WITH_A_DIRECTION | GC_RC_STABLE |
| Y4 chr13 gain | Y4.A0001 | 123.9036 | background_unknown_review | large_ge10mb | A_ANCHORED_WEAK_B_SIGNAL | GC_RC_AMPLIFIED |
| Y5 chr15 loss | Y5.A0003 | 42.9281 | background_unknown_review | reportable_candidate_ge2mb | B_SIGNAL_SUPPORTED_A_DIRECTION | GC_RC_ATTENUATED_SEVERE |
| Y6 chr7 loss | Y6.A0007 | 8.8477 | background_unknown_review | large_ge10mb | B_SIGNAL_SUPPORTED_A_DIRECTION | GC_RC_ATTENUATED |
| Y6 chr7 gain | Y6.A0003 | 94.0393 | background_unknown_review | large_ge10mb | B_SIGNAL_SUPPORTED_A_DIRECTION | GC_RC_STABLE |
| Y7 chr8 loss | Y7.A0001 | 40.8967 | background_unknown_review | large_ge10mb | B_SIGNAL_SUPPORTED_A_DIRECTION | GC_RC_ATTENUATED_SEVERE |
| Y8 chr13 gain | Y8.A0005 | 26.9159 | background_unknown_review | large_ge10mb | B_SIGNAL_SUPPORTED_A_DIRECTION | GC_RC_ATTENUATED_SEVERE |

Locked H1-H6 truth top candidates:

| truth | top candidate | A abs z | V2 disposition | length tier | B signal context | GC/RC context |
|---|---|---:|---|---|---|---|
| H1 chr16 gain | H1.A0004 | 75.4562 | background_unknown_review | large_ge10mb | B_SIGNAL_SUPPORTED_A_DIRECTION | GC_RC_STABLE |
| H2 chr6 gain | H2.A0005 | 54.7141 | background_unknown_review | large_ge10mb | A_ANCHORED_WEAK_B_SIGNAL | GC_RC_AMPLIFIED |
| H3 chr13 gain | H3.A0002 | 27.9770 | background_unknown_review | large_ge10mb | A_ANCHORED_WEAK_B_SIGNAL | GC_RC_AMPLIFIED |
| H4 chr15 gain | H4.A0001 | 111.5461 | background_unknown_review | large_ge10mb | A_ANCHORED_WEAK_B_SIGNAL | GC_RC_ATTENUATED_SEVERE |
| H4 chr21 gain | H4.A0002 | 75.1369 | background_unknown_review | large_ge10mb | B_SIGNAL_SUPPORTED_A_DIRECTION | GC_RC_ATTENUATED_SEVERE |
| H5 chrX loss part 1 | H5.A0001 | 62.1648 | sca_branch_s_review | large_ge10mb | B_SIGNAL_DISCORDANT_WITH_A_DIRECTION | GC_RC_STABLE |
| H5 chrX loss part 2 | H5.A0002 | 60.5780 | sca_branch_s_review | large_ge10mb | B_SIGNAL_DISCORDANT_WITH_A_DIRECTION | GC_RC_AMPLIFIED |
| H5 chrX loss part 3 | H5.A0002 | 60.5780 | sca_branch_s_review | large_ge10mb | B_SIGNAL_DISCORDANT_WITH_A_DIRECTION | GC_RC_AMPLIFIED |
| H6 chrX loss | H6.A0001 | 62.5303 | sca_branch_s_review | large_ge10mb | B_SIGNAL_DISCORDANT_WITH_A_DIRECTION | GC_RC_STABLE |
| H6 chr21 gain | H6.A0003 | 7.1135 | background_unknown_review | large_ge10mb | B_SIGNAL_SUPPORTED_A_DIRECTION | GC_RC_ATTENUATED |

H6 chr21 is retained. It is below the `A_STRONG_Z_GE_10` tier and must remain
eligible as sensitive/weak-positive evidence.

### Materialized V2 Context Counts

| cohort | background_unknown_review | sca_branch_s_review | B supported A | A-anchored weak B | B signal discordant | large >=10Mb | broad >=4Mb | reportable >=2Mb | review 1-2Mb |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Y1-Y8 | 84 | 13 | 66 | 20 | 11 | 85 | 3 | 5 | 4 |
| H1-H16 | 63 | 42 | 68 | 26 | 11 | 96 | 4 | 2 | 3 |
| 2026-06-15 | 151 | 14 | 97 | 40 | 28 | 137 | 11 | 9 | 8 |

0615 per-sample burden:

| sample | background_unknown_review | sca_branch_s_review |
|---|---:|---:|
| JZ26125843-56-56 | 24 | 2 |
| JZ26125844-59-59 | 34 | 1 |
| JZ26125845-60-60 | 42 | 5 |
| JZ26125846-61-61 | 30 | 3 |
| JZ26125847-62-62 | 21 | 3 |

## Current Decision

Branch B V2 has moved from legacy threshold/filter interpretation to
Branch-A-anchored evidence stratification. The current outputs are still
development-only and shadow-safe:

- no B-only final event;
- no GC/RC direction override;
- no legacy hard filter promotion;
- no final report release;
- no SCA final promotion.
