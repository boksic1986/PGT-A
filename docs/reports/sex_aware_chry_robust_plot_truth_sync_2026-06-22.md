# Sex-Aware chrY Robust Plot And Truth Sync

Date: 2026-06-22

Status: visualization/support-ledger update only

Reference: `h_r0_shadow_ref_20260619`

## Scope

This loop fixes CNV plot display and plot support ledgers. It does not change
Branch A, Branch B V2 filtering, Branch S classification, reference build,
mapping, report-event classification, or TP/FN/FP definitions.

## Changes

- `plot_bins.tsv` now separates raw z evidence from plot display:
  - `branch_a_ref_z`: raw ref-normalized z for audit.
  - `display_ref_z`: review z before plot clipping.
  - `plot_z`: SVG display z, clipped to `[-8, 8]`.
  - `z_plot_clipped`: marks bins clipped for display only.
- chrY is sex-aware:
  - XX samples keep chrY rows in TSV but hide ordinary z/CN scatter:
    `chrY_display_mode=xx_absent_expected_hidden`.
  - XY samples do not use chrY ref-z or ratio-CN as ordinary CNV evidence when
    the chrY reference denominator is not interpretable. They show a neutral
    Y-presence guide using BAM sex evidence:
    `chrY_display_mode=xy_y_presence_guide_not_ref_z` in z bins and
    `xy_y_presence_guide_not_ratio_cn` in CN bins.
- CN SVG no longer draws grey structural-gap background blocks. Structure/centromere
  bins are represented by absent scatter points. Chromosome blocks use light
  alternating backgrounds and darker separators.
- z SVG legend no longer contains `event ref-z trend`; only final dup/del/neutral
  point classes are shown.

## Remote Validation

Remote path:
`/data/project/CNV/PGT-A/refactor_validation_20260419`

Executables:

- Python: `/biosoftware/miniconda/envs/snakemake_env/bin/python`
- Snakemake: `/biosoftware/miniconda/envs/snakemake_env/bin/snakemake`

Tests:

```bash
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_plot.py \
  tests/unit/test_cnv_report.py \
  tests/unit/test_branch_s_shadow.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py -q
```

Result:

```text
63 passed in 2.63s
```

Dry-runs:

- `config_predict_y_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`
- `config_predict_h_20260608_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`
- `config_predict_g1_g8_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`
- `config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`

Targets:

```text
cnv_branch_s_shadow cnv_branch_ab_plot cnv_report
```

Result: all parsed successfully and did not request mapping/reference rebuild.

Materialization logs:

- `logs/manual/sex_aware_chry_plot_Y_20260622_175637.log`
- `logs/manual/sex_aware_chry_plot_H_20260622_175650.log`
- `logs/manual/sex_aware_chry_plot_G_20260622_175712.log`
- `logs/manual/sex_aware_chry_plot_0615_20260622_180102.log`

Results:

- Y1-Y8: `20 of 20 steps (100%) done`.
- H1-H16: materialized; H1-H6 used as truth gate.
- G1-G8: `20 of 20 steps (100%) done`.
- 0615 context: `14 of 14 steps (100%) done`.

## Acceptance Evidence

Output counts:

| cohort | z SVG | CN SVG | z TSV | CN TSV | support TSV |
|---|---:|---:|---:|---:|---:|
| Y1-Y8 | 8/8 | 8/8 | 8/8 | 8/8 | 8/8 |
| H1-H16 | 16/16 | 16/16 | 16/16 | 16/16 | 16/16 |
| G1-G8 | 8/8 | 8/8 | 8/8 | 8/8 | 8/8 |
| 0615 context | 5/5 | 5/5 | 5/5 | 5/5 | 5/5 |

Truth preservation summaries:

- Y1-Y8: `truth_preserved_count=10`.
- H1-H6: `truth_preserved_count=10`.
- G1-G8: `truth_preserved_count=10`.

Selected truth checks:

- H6 chr21 gain remains visible:
  `chr21:15000001-42000000 gain`, `a_zscore=7.113507`,
  `v2_report_visibility=report_weak_event`,
  `v2_filter_action=keep_report_layer_event`.
- G2 report rows remain visible for its broad truth-like events; only small
  combined-technical-risk rows are filtered.
- G3 Branch S X-loss support:
  `same_direction_ref_z_bin_count=170`,
  `same_direction_ref_z_fraction=0.977011`,
  `support_interpretation_status=Z_DIRECTION_SUPPORTED`.
- G5 Branch S X-loss support:
  `same_direction_ref_z_bin_count=169`,
  `same_direction_ref_z_fraction=0.971264`,
  `support_interpretation_status=Z_DIRECTION_SUPPORTED`.

chrY display mode checks:

- XX samples: chrY TSV rows retained, ordinary z/CN scatter hidden, mode
  `xx_absent_expected_hidden`.
- XY samples: chrY z guide visible at neutral guide position, CN guide visible
  at expected CN=1, modes `xy_y_presence_guide_not_ref_z` and
  `xy_y_presence_guide_not_ratio_cn`.

Local synchronized outputs:

- `D:\Pipeline\PGT-A\reports\truth_y1_y8_cnv_plots\` (`44` files)
- `D:\Pipeline\PGT-A\reports\truth_h1_h6_cnv_plots\` (`34` files)
- `D:\Pipeline\PGT-A\reports\truth_g1_g8_cnv_plots\` (`44` files)
- `D:\Pipeline\PGT-A\reports\0615_cnv_plots\` (`29` files, context only)

## Boundaries

- chrY guide is not a chrY CNV caller.
- `plot_z` clipping is display-only and does not enter filtering/calling.
- 0615 has no locked truth and remains burden/context only.
- This loop does not promote Branch S to final SCA reporting.
