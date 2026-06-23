# Sample Master Manifest Handoff

Date: 2026-06-24

Status: local_static_validated_pending_pr

## Goal

Create a single sample-level fact source for the current PGT-A / CNV-seq R&D
cycle so positive truth samples, presumed negatives, reference candidates, and
context-only samples are not reconstructed from scattered reports.

## Files Added

```text
sample_info/sample_master_manifest_20260624.tsv
docs/reports/sample_master_manifest_2026-06-24.md
```

## Files Updated

```text
CURRENT_STATE.md
PLANS.md
docs/CURRENT_CONTEXT_INDEX.md
```

## Key Decisions

- This is a sample fact-source change only.
- No Branch A, Branch B V2, Branch S, reference build, sex calling, report-event
  classification, filtering, or Snakemake production target was changed.
- `R0/R1/R2` remains reference-rebuild eligibility.
- `N0/N1/N2` remains Branch B background/negative-bank labeling.
- H9/H10/H11/H12/H15/H16 are included in the active shadow reference, but are
  not formal N0 controls.
- 2026-06-15 samples remain burden/context only because no locked truth labels
  exist.

## Current Manifest Summary

| group | rows |
|---|---:|
| Base reference-negative package samples | 32 |
| Y1-Y8 locked truth | 8 |
| H1-H6 locked truth | 6 |
| H7-H16 presumed negatives | 10 |
| G1-G8 locked truth | 8 |
| 2026-06-15 context samples | 5 |
| Total | 69 |

## Validation

Local static validation confirmed:

- no duplicate sample IDs;
- Y1-Y8, H1-H6, and G1-G8 are present as positive truth samples;
- H7-H16 are present with R/N labels;
- 2026-06-15 rows are context only;
- no sample is labeled N0.

No remote Snakemake or pytest run was required because this loop does not touch
workflow code, rules, configs, or report generation logic.

## Known Gap

`G1-G8` truth is currently sourced from:

```text
docs/reports/branch_b_v2_g1_g8_validation_2026-06-21.md
```

The report cites `sample_info/g1_g8_truth_events_20260621.tsv`, but that TSV is
not present in the local working tree. The manifest marks G1-G8 as
`locked_truth_report_only_pending_tsv`.

## Next Step

Create `sample_info/g1_g8_truth_events_20260621.tsv` and update the manifest to
point G1-G8 truth rows at the machine-readable TSV.
