# Report Main Convergence And CNV Plot Contract

Date: 2026-06-22

## Scope

This loop updates the report layer and CNV plot output only.

Unchanged:

- active 1Mb reference: `h_r0_shadow_ref_20260619`
- Branch A caller: WisecondorX/CBS with explicit `merge_gap_bp=2_000_000` overlay in the active validation configs
- Branch B V2: Branch-A-anchored report-layer evidence, no B-only events
- Branch S: review/development-only, no final SCA promotion
- reference build and BWA/mapping

## Implemented Contract

### Branch B V2 report visibility

Branch B V2 now exposes report visibility labels:

- `report_strong_event`
- `report_weak_event`
- `internal_review_event`
- `filtered_event`
- `branch_s_event`

Final autosomal CNV report event input is the V2 benchmark `report_events.tsv`,
which contains only `report_strong_event` and `report_weak_event` rows.
Internal review rows remain in V2/sample summaries and audit ledgers.
Filtered rows remain audit-only. Branch S rows remain in the SCA/sex-chromosome
section and do not enter the autosomal CNV main table.

The weak-report rule is conservative: a sensitive A signal is report-weak only
when it has same-direction B-side support, clean support, and same-direction
lowres support. This protects H6 chr21 while avoiding promotion of every
background-unknown review row to the main table.

### Report sample universe

`cnv_report_summary` now consumes Branch B V2 benchmark summary and sample
summary whenever `branch_b_v2_benchmark` is available in the workflow. It no
longer depends on `branch_b_v2_benchmark` being explicitly requested by the
user.

The report sample table is expanded from the V2 sample summary, so samples with
zero autosomal final report events remain visible with event counts set to zero.
This fixed the prior report-layer gap where Y3 and several H samples were absent
from `cnv_summary.tsv` when refreshing the report alone.

### WisecondorX-style CNV plot

`pgta/predict/branch_b/plot.py` now uses only `calibrated_z` as the bin-level
plot signal.

Rules:

- missing `calibrated_z` column is an error
- non-numeric or non-finite `calibrated_z` bins are skipped
- no fallback to `robust_z`, `raw_robust_z`, or normalized signal
- per-sample plot table is written to
  `wisecondorx/cnv/plots/{sample}.plot_bins.tsv`
- plot TSV columns include `chrom`, `start`, `end`, `genome_pos`, `z`,
  `report_state`
- `report_state` is limited to `dup`, `del`, `neutral`
- only final autosomal report events are highlighted
- SVG legend is limited to `dup`, `del`, `neutral bin`, `smooth z trend`

## Remote Validation

Remote path:
`/data/project/CNV/PGT-A/refactor_validation_20260419`

Executables:

- python: `/biosoftware/miniconda/envs/snakemake_env/bin/python`
- snakemake: `/biosoftware/miniconda/envs/snakemake_env/bin/snakemake`

Unit tests:

```text
tests/unit/test_cnv_report.py
tests/unit/test_branch_ab_phase12_workflow_contract.py
tests/unit/test_branch_b_plot.py
tests/unit/test_branch_b_v2_classifier.py
tests/unit/test_branch_b_v2_benchmark.py
```

Result:

```text
79 passed
```

Dry-run:

```text
branch_b_v2_benchmark branch_s_review cnv_report -n
```

All four active lowres-enabled gap2m configs parsed successfully:

- `config_predict_y_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`
- `config_predict_h_20260608_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`
- `config_predict_g1_g8_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`
- `config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`

Report refresh:

```text
snakemake --forcerun cnv_report_summary cnv_report
```

This was required because `cnv_report` is an aggregate target; forcing only the
aggregate target does not rerun `report.py`.

## Materialized Result Summary

| cohort | samples | candidates | truth | preserved | FN | truth filtered | report | strong | weak | internal | filtered | Branch S | SVG | plot TSV | report rows |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Y1-Y8 | 8 | 97 | 10 | 10 | 0 | 0 | 40 | 21 | 19 | 31 | 13 | 13 | 8 | 8 | 8 |
| H1-H16 | 16 | 105 | 10 | 10 | 0 | 0 | 23 | 6 | 17 | 32 | 8 | 42 | 16 | 16 | 16 |
| G1-G8 | 8 | 75 | 10 | 10 | 0 | 0 | 26 | 15 | 11 | 29 | 7 | 13 | 8 | 8 | 8 |
| 2026-06-15 | 5 | 165 | 0 | 0 | 0 | 0 | 71 | 52 | 19 | 57 | 23 | 14 | 5 | 5 | 5 |

Key locked truth checks:

- Y1-Y8: truth preserved 10/10, FN=0, truth filtered=0.
- H1-H6: truth preserved 10/10, FN=0, truth filtered=0.
- H6 chr21 gain remains `report_weak_event`,
  `keep_report_layer_event`, A z = 7.1135.
- H6 chrX remains routed to Branch S review.
- G1-G8: truth preserved 10/10, FN=0, truth filtered=0.
- G2 chr8 gain remains visible as `internal_review_event`; it is not filtered.
- G2 chr11 loss remains `report_strong_event`.
- 2026-06-15 has no locked truth and remains burden/context only.

Plot contract checks:

- all cohorts have matching `.final_cnv.svg` and `.plot_bins.tsv` per sample
- all plot TSV files contain required columns
- `report_state` values are restricted to `dup`, `del`, `neutral`
- SVG legend contains only `dup`, `del`, `neutral bin`, `smooth z trend`
- SVG does not contain Branch A/B, internal, filtered, mask, or rejected legend
  terms

## Interpretation

This loop makes the report package more usable:

- autosomal final report main table now contains only V2 report-layer events
- strong and weak report events are both visible
- internal review and filtered events are separated from the main report table
- Branch S remains separate from autosomal CNV
- zero-event samples remain present in the sample report summary
- each sample has a WisecondorX-style calibrated-z CNV plot and matching bin
  table

This is still a development report contract. It does not promote the shadow
reference, Branch B V2, Branch S, or lowres evidence to production-final.

## Remaining Work

The 2026-06-15 burden remains high, especially in samples 60 and 61. Further
reduction should not start with another broad threshold. The next gate should
inspect candidate-level evidence among remaining `report_strong_event` and
`report_weak_event` rows:

- same-direction lowres support
- ref-MAD stability
- clean support and region-risk
- B-side signal context
- SCA/sex-chromosome separation
- truth-protecting ablation across Y/H/G

Any future rule must keep Y/H/G FN=0 and must not filter H6 chr21 or G2 locked
truth.
