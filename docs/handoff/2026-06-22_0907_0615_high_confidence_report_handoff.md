# Handoff: 0615 High-Confidence Report Candidate Review

Date: 2026-06-22 09:07

## Context Source

Active previous handoff:
`docs/handoff/2026-06-22_0437_report_main_cnv_plot_handoff.md`

New report:
`docs/reports/0615_high_confidence_report_candidates_2026-06-22.md`

## Scope Completed

This loop recorded a read-only review of the current 2026-06-15 report outputs.
No workflow code, thresholds, Branch A, Branch B V2, Branch S, reference, or
remote result files were modified.

## Remote Inputs

Remote result directory:

`/data/project/CNV/PGT-A/refactor_validation_20260419/results_20260615_h_r0_shadow_ref_20260619`

Inputs reviewed:

- `wisecondorx/cnv/postprocess_gap2m/v2_benchmark/report_events.tsv`
- `wisecondorx/cnv/plots/{sample}.plot_bins.tsv`
- `wisecondorx/cnv/postprocess_gap2m/v2_classifier/*.candidate_classification.tsv`
- `wisecondorx/cnv/report/cnv_summary.tsv`
- `wisecondorx/cnv/postprocess_gap2m/branch_s/*.summary.json`

## Remote Read-Only Check

Executable:

`/biosoftware/miniconda/envs/snakemake_env/bin/python`

Result:

- `report_events.tsv`: 71 autosomal report rows.
- Plot bin TSV files: 5/5 samples present.
- Current plot table columns: `chrom`, `start`, `end`, `genome_pos`, `z`,
  `report_state`.

## Key Result

2026-06-15 has no locked truth labels. The output is therefore
high-confidence review/report candidate context only; no TP/FP/FN conclusion is
allowed.

Conservative high-confidence autosomal candidates:

- 10 rows, all in `JZ26125845-60-60`.
- Other samples have no current conservative high-confidence autosomal row.

Per-sample summary:

| sample | original report events | mechanical pass | conservative high-confidence | batch-shared review | Branch S context |
|---|---:|---:|---:|---:|---:|
| JZ26125843-56-56 | 7 | 1 | 0 | 0 | 2 |
| JZ26125844-59-59 | 11 | 0 | 0 | 0 | 1 |
| JZ26125845-60-60 | 23 | 10 | 10 | 4 | 5 |
| JZ26125846-61-61 | 18 | 3 | 0 | 2 | 3 |
| JZ26125847-62-62 | 12 | 1 | 0 | 2 | 3 |

## Batch-Shared Context

Current classifier-derived shared regions:

- `chr4:67.50-101.25Mb gain`, 5/5 samples: excluded from high-confidence.
- `chr4:52.50-67.50Mb gain`, 4/5 samples: batch-shared review.
- `chr4:101.25-121.50Mb gain`, 4/5 samples: batch-shared review.
- `chr14:60.75-97.50Mb gain`, 4/5 samples: batch-shared review.

## High-Confidence Rows

`JZ26125845-60-60` current conservative high-confidence autosomal candidates:

- `chr1:30.75-118.50Mb gain`
- `chr2:6.00-87.00Mb gain`
- `chr7:143.25-159.00Mb gain`
- `chr8:3.00-42.75Mb loss`
- `chr9:5.25-39.00Mb loss`
- `chr10:18.75-36.00Mb gain`
- `chr12:1.50-33.75Mb loss`
- `chr17:1.50-21.00Mb gain`
- `chr18:23.25-72.00Mb loss`
- `chr21:15.00-42.00Mb loss`

## Important Boundary

This review is not a new filtering rule set. It should not be used to
reverse-engineer Branch B V2 thresholds from 0615 burden. Any future workflow
change still needs locked-truth ablation against Y1-Y8, H1-H6, and G1-G8 while
preserving H6 chr21 and G2 truth visibility.

## Next Recommended Step

Use the high-confidence candidate table to guide manual review of the 0615
development report and CNV plots. Keep chr4 shared gains as batch-effect review
context. Do not promote 0615 report, Branch B V2, Branch S, or
`h_r0_shadow_ref_20260619` to production-final from this evidence alone.
