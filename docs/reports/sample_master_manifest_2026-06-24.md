# Sample Master Manifest 2026-06-24

## Scope

This report introduces a single sample-status manifest for the current PGT-A /
CNV-seq R&D cycle:

```text
sample_info/sample_master_manifest_20260624.tsv
```

The manifest is a sample fact source only. It does not change Branch A, Branch B
V2, Branch S, reference build, sex calling, event classification, report-event
filtering, or any Snakemake production target.

## Manifest Contract

Each row is one sample. Key fields:

- `declared_role`: broad biological or operational role.
- `current_use_role`: how the sample may currently be used in this project.
- `reference_label`: reference-rebuild eligibility such as `R0` or `R1`.
- `negative_bank_label`: Branch B background/negative-bank label such as `N1`
  or `N2`.
- `included_in_current_ref`: whether the sample is included in the active shadow
  reference `h_r0_shadow_ref_20260619`.
- `allowed_use` and `forbidden_use`: explicit use boundaries.

Allowed role values used in this first manifest:

```text
declared_role:
  positive_truth
  presumed_negative
  reference_negative
  context_unknown

current_use_role:
  truth_benchmark
  shadow_reference_candidate
  burden_context_only
```

## Current Counts

| group | rows | notes |
|---|---:|---|
| Base reference-negative package samples | 32 | Included in the current shadow reference; not formal N0 controls. |
| Y1-Y8 locked truth | 8 | Machine-readable truth source is `sample_info/positive_truth_events.tsv`. |
| H1-H6 locked truth | 6 | Machine-readable truth source is `sample_info/positive_truth_events.tsv`. |
| H7-H16 presumed negatives | 10 | R/N labels are tracked separately; H9-H12/H15/H16 are included in the current shadow reference. |
| G1-G8 locked truth | 8 | Current truth source is report-only pending a local TSV. |
| 2026-06-15 context samples | 5 | No locked truth; burden/context only. |
| Total | 69 | One row per sample. |

## Positive Truth Sources

`sample_info/positive_truth_events.tsv` is the current machine-readable truth
source for:

- `Y1-Y8`
- `H1-H6`

`sample_info/mosaic_fraction_validation_y7_y8.tsv` remains an auxiliary mosaic
context file for Y7/Y8.

`G1-G8` is recorded from:

```text
docs/reports/branch_b_v2_g1_g8_validation_2026-06-21.md
```

That report cites:

```text
sample_info/g1_g8_truth_events_20260621.tsv
```

but the TSV was not present in the local working tree during this manifest
creation. Therefore G1-G8 rows are marked:

```text
status=locked_truth_report_only_pending_tsv
truth_source=docs/reports/branch_b_v2_g1_g8_validation_2026-06-21.md
```

This avoids pretending the missing machine-readable truth TSV exists.

## Negative And Reference Sources

The base reference-negative package is documented in:

```text
docs/reports/reference_negative_package_sync_2026-06-16.md
```

The committed report records the package-level facts:

- sample count: 32
- FASTQ file count: 64
- XX count: 20
- XY count: 12
- MD5 verified: `fa7036739ecec4685364734762d55d08`

The active shadow reference selected sample list on remote contains those 32
base samples plus H9, H10, H11, H12, H15, and H16, for 38 selected reference
samples total. The manifest records these rows as:

```text
included_in_current_ref=yes
status=active_shadow_reference or shadow_ref_included_not_N0
```

This does not promote them to production reference status and does not make them
formal N0 negatives.

## R Labels Versus N Labels

`R0/R1/R2` and `N0/N1/N2` are intentionally separate:

- `R0/R1/R2` means reference-rebuild eligibility.
- `N0/N1/N2` means Branch B background or negative calibration label.

For H7-H16, the current reference decision source is:

```text
docs/reports/h7_h16_reference_cohort_decision_2026-06-20.tsv
```

The current negative-bank seed source is:

```text
docs/reports/branch_ab_v2_negative_bank_seed_2026-06-18.tsv
```

The first-pass R0 shadow reference additions are:

```text
H9,H10,H11,H12,H15,H16
```

They are included in `h_r0_shadow_ref_20260619`, but their `N1` labels do not
make them matched-clean N0 controls.

## 2026-06-15 Samples

The 2026-06-15 samples are recorded as context only. They have no locked truth
labels in this repository, so they must not be used for:

- TP/FN/FP performance claims
- threshold derivation
- reference-negative promotion
- formal N0 background construction

Current manifest roles:

- `JZ26125843-56-56`, `JZ26125844-59-59`, `JZ26125846-61-61`,
  `JZ26125847-62-62`: `declared_role=presumed_negative`,
  `current_use_role=burden_context_only`
- `JZ26125845-60-60`: `declared_role=context_unknown`,
  `current_use_role=burden_context_only`

The 60 sample is explicitly not locked truth; it remains a quality/review
context sample.

## Validation

Static validation was performed on `sample_info/sample_master_manifest_20260624.tsv`:

- no duplicate `sample_id`;
- all `Y1-Y8`, `H1-H6`, and `G1-G8` rows are present and marked
  `declared_role=positive_truth`;
- all H7-H16 rows are present;
- H7-H16 R labels come from the 2026-06-20 reference decision table;
- H7-H16 N labels come from the 2026-06-18 negative-bank seed table;
- 2026-06-15 rows have `truth_event_count=0` and `current_use_role=burden_context_only`;
- no row has `negative_bank_label=N0`.

## Next Step

Create the missing machine-readable G1-G8 truth file:

```text
sample_info/g1_g8_truth_events_20260621.tsv
```

Then update the G1-G8 rows in the master manifest from
`status=locked_truth_report_only_pending_tsv` to a machine-readable truth source
status.
