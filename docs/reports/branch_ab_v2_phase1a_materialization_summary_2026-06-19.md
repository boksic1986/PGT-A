# Branch A/B V2 Phase 1A Materialization Summary

Date: 2026-06-19

## Scope

This report materializes the CNVpro-like Branch B V2 shadow evidence ledger fields for Y1-Y8, H1-H16, and the 2026-06-15 shadow-report set. It is evidence-only and does not change WisecondorX predict, Branch A candidate inclusion, Branch B final report decisions, or report schema.

## Remote Inputs

| dataset                         | truth_scope                                                 |   ledger_files |   classifier_files |   ledger_rows |   classifier_rows |
|:--------------------------------|:------------------------------------------------------------|---------------:|-------------------:|--------------:|------------------:|
| Y1-Y8 locked positives          | locked_positive                                             |              8 |                  8 |           143 |               143 |
| H1-H6 positives + H7-H16 review | H1-H6 locked_positive; H7-H16 presumed_negative_review_only |             16 |                 16 |           307 |               307 |
| 2026-06-15 shadow-report set    | no_locked_truth                                             |              5 |                  5 |           170 |               170 |


## Current Final Truth Metrics

These metrics come from the existing final/report evaluation files, not from promoting Phase 1A fields into decision logic.

| dataset                         |   truth_events |   matched |   detected |   min_top_support_z |   max_top_support_z |
|:--------------------------------|---------------:|----------:|-----------:|--------------------:|--------------------:|
| Y1-Y8 locked positives          |             10 |        10 |         10 |              9.3707 |             223.379 |
| H1-H6 positives + H7-H16 review |             10 |        10 |         10 |              7.3889 |             130.437 |


## Per-sample Candidate Burden

| dataset                         | sample_id        | cohort_role              |   a_candidate_count |   truth_overlap_candidate_count | existing_final_disposition_counts      | v2_candidate_class_counts                     |
|:--------------------------------|:-----------------|:-------------------------|--------------------:|--------------------------------:|:---------------------------------------|:----------------------------------------------|
| 2026-06-15 shadow-report set    | JZ26125843-56-56 | shadow_unknown           |                  29 |                               0 | LIKELY_ARTIFACT=27; REVIEW_REQUIRED=2  | LIKELY_ARTIFACT_SHADOW=27; REVIEW_REQUIRED=2  |
| 2026-06-15 shadow-report set    | JZ26125844-59-59 | shadow_unknown           |                  29 |                               0 | LIKELY_ARTIFACT=27; REVIEW_REQUIRED=2  | LIKELY_ARTIFACT_SHADOW=27; REVIEW_REQUIRED=2  |
| 2026-06-15 shadow-report set    | JZ26125845-60-60 | shadow_unknown           |                  47 |                               0 | LIKELY_ARTIFACT=27; REVIEW_REQUIRED=20 | LIKELY_ARTIFACT_SHADOW=27; REVIEW_REQUIRED=20 |
| 2026-06-15 shadow-report set    | JZ26125846-61-61 | shadow_unknown           |                  38 |                               0 | LIKELY_ARTIFACT=36; REVIEW_REQUIRED=2  | LIKELY_ARTIFACT_SHADOW=36; REVIEW_REQUIRED=2  |
| 2026-06-15 shadow-report set    | JZ26125847-62-62 | shadow_unknown           |                  27 |                               0 | LIKELY_ARTIFACT=25; REVIEW_REQUIRED=2  | LIKELY_ARTIFACT_SHADOW=25; REVIEW_REQUIRED=2  |
| H1-H6 positives + H7-H16 review | H1               | locked_positive          |                   9 |                               2 | LIKELY_ARTIFACT=6; REVIEW_REQUIRED=3   | LIKELY_ARTIFACT_SHADOW=6; REVIEW_REQUIRED=3   |
| H1-H6 positives + H7-H16 review | H10              | presumed_negative_review |                  18 |                               0 | LIKELY_ARTIFACT=18                     | LIKELY_ARTIFACT_SHADOW=18                     |
| H1-H6 positives + H7-H16 review | H11              | presumed_negative_review |                  19 |                               0 | LIKELY_ARTIFACT=19                     | LIKELY_ARTIFACT_SHADOW=19                     |
| H1-H6 positives + H7-H16 review | H12              | presumed_negative_review |                  22 |                               0 | LIKELY_ARTIFACT=22                     | LIKELY_ARTIFACT_SHADOW=22                     |
| H1-H6 positives + H7-H16 review | H13              | presumed_negative_review |                  30 |                               0 | LIKELY_ARTIFACT=28; REVIEW_REQUIRED=2  | LIKELY_ARTIFACT_SHADOW=28; REVIEW_REQUIRED=2  |
| H1-H6 positives + H7-H16 review | H14              | presumed_negative_review |                  22 |                               0 | LIKELY_ARTIFACT=20; REVIEW_REQUIRED=2  | LIKELY_ARTIFACT_SHADOW=20; REVIEW_REQUIRED=2  |
| H1-H6 positives + H7-H16 review | H15              | presumed_negative_review |                  15 |                               0 | LIKELY_ARTIFACT=15                     | LIKELY_ARTIFACT_SHADOW=15                     |
| H1-H6 positives + H7-H16 review | H16              | presumed_negative_review |                  19 |                               0 | LIKELY_ARTIFACT=19                     | LIKELY_ARTIFACT_SHADOW=19                     |
| H1-H6 positives + H7-H16 review | H2               | locked_positive          |                  31 |                               3 | LIKELY_ARTIFACT=28; REVIEW_REQUIRED=3  | LIKELY_ARTIFACT_SHADOW=28; REVIEW_REQUIRED=3  |
| H1-H6 positives + H7-H16 review | H3               | locked_positive          |                  17 |                               1 | LIKELY_ARTIFACT=16; REVIEW_REQUIRED=1  | LIKELY_ARTIFACT_SHADOW=16; REVIEW_REQUIRED=1  |
| H1-H6 positives + H7-H16 review | H4               | locked_positive          |                  19 |                               3 | LIKELY_ARTIFACT=16; REVIEW_REQUIRED=3  | LIKELY_ARTIFACT_SHADOW=16; REVIEW_REQUIRED=3  |
| H1-H6 positives + H7-H16 review | H5               | locked_positive          |                  21 |                               2 | LIKELY_ARTIFACT=19; REVIEW_REQUIRED=2  | LIKELY_ARTIFACT_SHADOW=19; REVIEW_REQUIRED=2  |
| H1-H6 positives + H7-H16 review | H6               | locked_positive          |                  17 |                               3 | LIKELY_ARTIFACT=14; REVIEW_REQUIRED=3  | LIKELY_ARTIFACT_SHADOW=14; REVIEW_REQUIRED=3  |
| H1-H6 positives + H7-H16 review | H7               | presumed_negative_review |                  13 |                               0 | LIKELY_ARTIFACT=13                     | LIKELY_ARTIFACT_SHADOW=13                     |
| H1-H6 positives + H7-H16 review | H8               | presumed_negative_review |                  17 |                               0 | LIKELY_ARTIFACT=16; REVIEW_REQUIRED=1  | LIKELY_ARTIFACT_SHADOW=16; REVIEW_REQUIRED=1  |
| H1-H6 positives + H7-H16 review | H9               | presumed_negative_review |                  18 |                               0 | LIKELY_ARTIFACT=18                     | LIKELY_ARTIFACT_SHADOW=18                     |
| Y1-Y8 locked positives          | Y1               | locked_positive          |                  17 |                               2 | LIKELY_ARTIFACT=16; REVIEW_REQUIRED=1  | REVIEW_REQUIRED=17                            |
| Y1-Y8 locked positives          | Y2               | locked_positive          |                  22 |                               3 | LIKELY_ARTIFACT=19; REVIEW_REQUIRED=3  | REVIEW_REQUIRED=22                            |
| Y1-Y8 locked positives          | Y3               | locked_positive          |                  15 |                               2 | LIKELY_ARTIFACT=13; REVIEW_REQUIRED=2  | REVIEW_REQUIRED=15                            |
| Y1-Y8 locked positives          | Y4               | locked_positive          |                  25 |                               1 | LIKELY_ARTIFACT=22; REVIEW_REQUIRED=3  | REVIEW_REQUIRED=25                            |
| Y1-Y8 locked positives          | Y5               | locked_positive          |                   9 |                               1 | LIKELY_ARTIFACT=6; REVIEW_REQUIRED=3   | REVIEW_REQUIRED=9                             |
| Y1-Y8 locked positives          | Y6               | locked_positive          |                  15 |                               2 | LIKELY_ARTIFACT=12; REVIEW_REQUIRED=3  | REVIEW_REQUIRED=15                            |
| Y1-Y8 locked positives          | Y7               | locked_positive          |                  19 |                               1 | LIKELY_ARTIFACT=17; REVIEW_REQUIRED=2  | REVIEW_REQUIRED=19                            |
| Y1-Y8 locked positives          | Y8               | locked_positive          |                  21 |                               1 | LIKELY_ARTIFACT=18; REVIEW_REQUIRED=3  | REVIEW_REQUIRED=21                            |


## Key Field Missingness

| dataset                         | cohort_role              | field                               |   row_count |   missing_or_unknown_count |   missing_or_unknown_fraction |
|:--------------------------------|:-------------------------|:------------------------------------|------------:|---------------------------:|------------------------------:|
| 2026-06-15 shadow-report set    | shadow_unknown           | cnvpro_like_gc_rc_background_status |         170 |                          0 |                             0 |
| 2026-06-15 shadow-report set    | shadow_unknown           | dynamic_reference_status            |         170 |                          0 |                             0 |
| 2026-06-15 shadow-report set    | shadow_unknown           | matched_negative_source             |         170 |                        170 |                             1 |
| 2026-06-15 shadow-report set    | shadow_unknown           | matched_negative_percentile         |         170 |                        170 |                             1 |
| 2026-06-15 shadow-report set    | shadow_unknown           | copy_number_estimate                |         170 |                          0 |                             0 |
| 2026-06-15 shadow-report set    | shadow_unknown           | mosaic_fraction_proxy               |         170 |                          0 |                             0 |
| 2026-06-15 shadow-report set    | shadow_unknown           | cnvpro_like_evidence_status         |         170 |                          0 |                             0 |
| H1-H6 positives + H7-H16 review | locked_positive          | cnvpro_like_gc_rc_background_status |         114 |                          0 |                             0 |
| H1-H6 positives + H7-H16 review | locked_positive          | dynamic_reference_status            |         114 |                          0 |                             0 |
| H1-H6 positives + H7-H16 review | locked_positive          | matched_negative_source             |         114 |                        114 |                             1 |
| H1-H6 positives + H7-H16 review | locked_positive          | matched_negative_percentile         |         114 |                        114 |                             1 |
| H1-H6 positives + H7-H16 review | locked_positive          | copy_number_estimate                |         114 |                          0 |                             0 |
| H1-H6 positives + H7-H16 review | locked_positive          | mosaic_fraction_proxy               |         114 |                          0 |                             0 |
| H1-H6 positives + H7-H16 review | locked_positive          | cnvpro_like_evidence_status         |         114 |                          0 |                             0 |
| H1-H6 positives + H7-H16 review | presumed_negative_review | cnvpro_like_gc_rc_background_status |         193 |                          0 |                             0 |
| H1-H6 positives + H7-H16 review | presumed_negative_review | dynamic_reference_status            |         193 |                          0 |                             0 |
| H1-H6 positives + H7-H16 review | presumed_negative_review | matched_negative_source             |         193 |                        193 |                             1 |
| H1-H6 positives + H7-H16 review | presumed_negative_review | matched_negative_percentile         |         193 |                        193 |                             1 |
| H1-H6 positives + H7-H16 review | presumed_negative_review | copy_number_estimate                |         193 |                          0 |                             0 |
| H1-H6 positives + H7-H16 review | presumed_negative_review | mosaic_fraction_proxy               |         193 |                          0 |                             0 |
| H1-H6 positives + H7-H16 review | presumed_negative_review | cnvpro_like_evidence_status         |         193 |                          0 |                             0 |
| Y1-Y8 locked positives          | locked_positive          | cnvpro_like_gc_rc_background_status |         143 |                          0 |                             0 |
| Y1-Y8 locked positives          | locked_positive          | dynamic_reference_status            |         143 |                          0 |                             0 |
| Y1-Y8 locked positives          | locked_positive          | matched_negative_source             |         143 |                        143 |                             1 |
| Y1-Y8 locked positives          | locked_positive          | matched_negative_percentile         |         143 |                        143 |                             1 |
| Y1-Y8 locked positives          | locked_positive          | copy_number_estimate                |         143 |                          0 |                             0 |
| Y1-Y8 locked positives          | locked_positive          | mosaic_fraction_proxy               |         143 |                          0 |                             0 |
| Y1-Y8 locked positives          | locked_positive          | cnvpro_like_evidence_status         |         143 |                          0 |                             0 |


## Key Status Field Distributions

| dataset                         | cohort_role              | field                               | value_counts             |
|:--------------------------------|:-------------------------|:------------------------------------|:-------------------------|
| 2026-06-15 shadow-report set    | shadow_unknown           | cnvpro_like_gc_rc_background_status | GC_RC_AVAILABLE=170      |
| 2026-06-15 shadow-report set    | shadow_unknown           | dynamic_reference_status            | OK=170                   |
| 2026-06-15 shadow-report set    | shadow_unknown           | matched_negative_source             | UNKNOWN_BACKGROUND=170   |
| 2026-06-15 shadow-report set    | shadow_unknown           | cnvpro_like_evidence_status         | SHADOW_EVIDENCE_ONLY=170 |
| H1-H6 positives + H7-H16 review | locked_positive          | cnvpro_like_gc_rc_background_status | GC_RC_AVAILABLE=114      |
| H1-H6 positives + H7-H16 review | locked_positive          | dynamic_reference_status            | OK=114                   |
| H1-H6 positives + H7-H16 review | locked_positive          | matched_negative_source             | UNKNOWN_BACKGROUND=114   |
| H1-H6 positives + H7-H16 review | locked_positive          | cnvpro_like_evidence_status         | SHADOW_EVIDENCE_ONLY=114 |
| H1-H6 positives + H7-H16 review | presumed_negative_review | cnvpro_like_gc_rc_background_status | GC_RC_AVAILABLE=193      |
| H1-H6 positives + H7-H16 review | presumed_negative_review | dynamic_reference_status            | OK=193                   |
| H1-H6 positives + H7-H16 review | presumed_negative_review | matched_negative_source             | UNKNOWN_BACKGROUND=193   |
| H1-H6 positives + H7-H16 review | presumed_negative_review | cnvpro_like_evidence_status         | SHADOW_EVIDENCE_ONLY=193 |
| Y1-Y8 locked positives          | locked_positive          | cnvpro_like_gc_rc_background_status | GC_RC_AVAILABLE=143      |
| Y1-Y8 locked positives          | locked_positive          | dynamic_reference_status            | OK=143                   |
| Y1-Y8 locked positives          | locked_positive          | matched_negative_source             | UNKNOWN_BACKGROUND=143   |
| Y1-Y8 locked positives          | locked_positive          | cnvpro_like_evidence_status         | SHADOW_EVIDENCE_ONLY=143 |


## Known-positive Truth-overlap Field Risk

A field is marked risky if any known-positive truth-overlap candidate has a missing/unknown value for that field. This does not mean the event is false; it means that field is unsafe as a hard filter without more evidence.

| dataset                         | field                               |   truth_overlap_candidate_count |   truth_overlap_missing_or_unknown_count |   truth_overlap_missing_or_unknown_fraction | hard_filter_risk_note                  |
|:--------------------------------|:------------------------------------|--------------------------------:|-----------------------------------------:|--------------------------------------------:|:---------------------------------------|
| H1-H6 positives + H7-H16 review | cnvpro_like_gc_rc_background_status |                              14 |                                        0 |                                           0 | no_missing_on_truth_overlap_candidates |
| H1-H6 positives + H7-H16 review | dynamic_reference_status            |                              14 |                                        0 |                                           0 | no_missing_on_truth_overlap_candidates |
| H1-H6 positives + H7-H16 review | matched_negative_source             |                              14 |                                       14 |                                           1 | risky_if_used_as_hard_filter           |
| H1-H6 positives + H7-H16 review | matched_negative_percentile         |                              14 |                                       14 |                                           1 | risky_if_used_as_hard_filter           |
| H1-H6 positives + H7-H16 review | copy_number_estimate                |                              14 |                                        0 |                                           0 | no_missing_on_truth_overlap_candidates |
| H1-H6 positives + H7-H16 review | mosaic_fraction_proxy               |                              14 |                                        0 |                                           0 | no_missing_on_truth_overlap_candidates |
| H1-H6 positives + H7-H16 review | cnvpro_like_evidence_status         |                              14 |                                        0 |                                           0 | no_missing_on_truth_overlap_candidates |
| Y1-Y8 locked positives          | cnvpro_like_gc_rc_background_status |                              13 |                                        0 |                                           0 | no_missing_on_truth_overlap_candidates |
| Y1-Y8 locked positives          | dynamic_reference_status            |                              13 |                                        0 |                                           0 | no_missing_on_truth_overlap_candidates |
| Y1-Y8 locked positives          | matched_negative_source             |                              13 |                                       13 |                                           1 | risky_if_used_as_hard_filter           |
| Y1-Y8 locked positives          | matched_negative_percentile         |                              13 |                                       13 |                                           1 | risky_if_used_as_hard_filter           |
| Y1-Y8 locked positives          | copy_number_estimate                |                              13 |                                        0 |                                           0 | no_missing_on_truth_overlap_candidates |
| Y1-Y8 locked positives          | mosaic_fraction_proxy               |                              13 |                                        0 |                                           0 | no_missing_on_truth_overlap_candidates |
| Y1-Y8 locked positives          | cnvpro_like_evidence_status         |                              13 |                                        0 |                                           0 | no_missing_on_truth_overlap_candidates |


## Output Tables

- /data/project/CNV/PGT-A/refactor_validation_20260419/docs/reports/branch_ab_v2_phase1a_materialization_summary_2026-06-19.burden.tsv
- /data/project/CNV/PGT-A/refactor_validation_20260419/docs/reports/branch_ab_v2_phase1a_materialization_summary_2026-06-19.missingness.tsv
- /data/project/CNV/PGT-A/refactor_validation_20260419/docs/reports/branch_ab_v2_phase1a_materialization_summary_2026-06-19.truth_overlap_field_risk.tsv
- /data/project/CNV/PGT-A/refactor_validation_20260419/docs/reports/branch_ab_v2_phase1a_materialization_summary_2026-06-19.inputs.json

## Interpretation

- Phase 1A materialization is complete for all requested cohorts: 8 Y ledgers, 16 H ledgers, and 5 2026-06-15 ledgers.
- Existing final/report truth metrics remain FN=0 for Y1-Y8 and H1-H6 in the current outputs.
- matched_negative_source and matched_negative_percentile are missing/unknown on all known-positive truth-overlap candidates, so matched-negative evidence is not ready for hard filtering.
- copy_number_estimate, mosaic_fraction_proxy, GC/RC background status, dynamic-reference status, and CNVpro-like evidence status are present on known-positive truth-overlap candidates, but remain shadow-only until ablation proves no FN increase.
