# Branch B V2 G1-G8 Validation

Date: 2026-06-21

Status: `development_validation`

## Scope

This report records the current-scheme validation of the earlier G1-G8 positive
samples. It uses the active development contract:

- reference: `h_r0_shadow_ref_20260619`
- Branch A overlay: `merge_gap_bp=2_000_000`, `strong_z=10.0`
- postprocess path: `wisecondorx/cnv/postprocess_gap2m`
- Branch B V2 report-layer classes:
  `report_event`, `internal_review_event`, `filtered_event`, `branch_s_event`

This is a validation/readout report. It does not promote the shadow reference,
Branch B V2, or Branch S to production-final status.

## Inputs

Existing G-batch BAM files were reused from:

```text
/data/project/CNV/PGT-A/g_reseq_qc_20260504/mapping/G*.sorted.bam
```

No BWA remapping was performed for this validation. WisecondorX convert,
WisecondorX predict, Branch A, Branch B V2, Branch S, benchmark, and report were
rerun under the current configuration.

Remote working directory:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419
```

Runtime configs:

```text
config_predict_g1_g8_h_r0_shadow_20260621.yaml
config_predict_g1_g8_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml
```

Truth table:

```text
sample_info/g1_g8_truth_events_20260621.tsv
```

Result root:

```text
results_g1_g8_h_r0_shadow_ref_20260621
```

## Overall Result

| Metric | Count |
|---|---:|
| Samples | 8 |
| Branch A candidates | 75 |
| Locked truth events | 10 |
| Truth preserved | 10 |
| FN | 0 |
| Truth hard suppressed | 0 |
| Truth report-layer filtered | 0 |
| Final report events | 15 |
| Internal review events | 40 |
| Filtered audit-only events | 7 |
| Branch S events | 13 |

Interpretation:

- The current sensitivity/preservation gate passes for G1-G8:
  `truth_preserved=10/10`, `FN=0`, and `truth_report_layer_filtered=0`.
- The result supports the current no-FN direction for Branch A plus Branch B V2
  report-layer filtering.
- It does not prove final specificity or production readiness. G2 still has a
  high candidate/review burden and should remain a focus for the next
  truth-safe demotion pass.

## Per-Sample Summary

| Sample | Candidates | Truth | Preserved | FN | Report events | Internal review | Filtered audit-only | Branch S |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| G1 | 4 | 1 | 1 | 0 | 1 | 3 | 0 | 0 |
| G2 | 28 | 2 | 2 | 0 | 6 | 16 | 3 | 3 |
| G3 | 13 | 2 | 2 | 0 | 1 | 10 | 0 | 2 |
| G4 | 4 | 1 | 1 | 0 | 2 | 2 | 0 | 0 |
| G5 | 8 | 1 | 1 | 0 | 1 | 3 | 2 | 2 |
| G6 | 8 | 1 | 1 | 0 | 2 | 3 | 0 | 3 |
| G7 | 4 | 1 | 1 | 0 | 1 | 1 | 2 | 0 |
| G8 | 6 | 1 | 1 | 0 | 1 | 2 | 0 | 3 |

## Truth Event Preservation

| Sample | Truth event | Top candidate | A abs z | V2 report-layer class | Visibility | Notes |
|---|---|---|---:|---|---|---|
| G1 | chr18 gain | G1.A0001 | 137.89 | `report_event` | `final_report` | Preserved |
| G2 | chr8 gain, 3.69-5.92 Mb | G2.A0004 | 29.42 | `internal_review_event` | `internal_review` | Preserved, B-side weak support |
| G2 | chr11 loss, 129.91-134.94 Mb | G2.A0005 | 25.69 | `report_event` | `final_report` | Preserved |
| G3 | chr16 gain | G3.A0001 | 85.58 | `report_event` | `final_report` | Preserved |
| G3 | chrX loss | G3.A0003 | 24.16 | `branch_s_event` | `branch_s_report_section` | Preserved, routed to Branch S |
| G4 | chr5 loss, 0.04-32.21 Mb | G4.A0001 | 69.50 | `report_event` | `final_report` | Preserved |
| G5 | chrX mosaic loss | G5.A0003 | 14.41 | `branch_s_event` | `branch_s_report_section` | Preserved, routed to Branch S |
| G6 | chr9 gain | G6.A0003 | 116.13 | `report_event` | `final_report` | Preserved |
| G7 | chr22 gain | G7.A0001 | 65.97 | `report_event` | `final_report` | Preserved |
| G8 | chr16 mosaic gain | G8.A0004 | 62.05 | `report_event` | `final_report` | Preserved |

## Filtered Events

Seven candidates were filtered from the report/internal-review main flow and
kept only in the audit ledger. None overlaps the locked G1-G8 truth events.

| Sample | Filtered count | Main reason |
|---|---:|---|
| G2 | 3 | combined sensitive/short, B-side non-support, GC/RC attenuation |
| G5 | 2 | B-side non-support plus GC/RC attenuation |
| G7 | 2 | B-side non-support plus GC/RC attenuation |

## Key Output Files

```text
results_g1_g8_h_r0_shadow_ref_20260621/wisecondorx/cnv/postprocess_gap2m/v2_benchmark/summary.json
results_g1_g8_h_r0_shadow_ref_20260621/wisecondorx/cnv/postprocess_gap2m/v2_benchmark/sample_summary.tsv
results_g1_g8_h_r0_shadow_ref_20260621/wisecondorx/cnv/postprocess_gap2m/v2_benchmark/truth_metrics.tsv
results_g1_g8_h_r0_shadow_ref_20260621/wisecondorx/cnv/postprocess_gap2m/v2_benchmark/filtered_events.tsv
results_g1_g8_h_r0_shadow_ref_20260621/wisecondorx/cnv/report/cnv_summary.tsv
results_g1_g8_h_r0_shadow_ref_20260621/wisecondorx/cnv/report/cnv_summary.json
```

## Commands

Dry-run:

```bash
cd /data/project/CNV/PGT-A/refactor_validation_20260419
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile config_predict_g1_g8_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml \
  --cores 1 -n branch_b_v2_benchmark branch_s_review cnv_report
```

Result: passed; 127 jobs planned for G1-G8 using existing BAM inputs.

Materialization:

```bash
cd /data/project/CNV/PGT-A/refactor_validation_20260419
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile config_predict_g1_g8_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml \
  --cores 32 branch_b_v2_benchmark branch_s_review cnv_report
```

Result: completed; `127 of 127 steps (100%) done`.

Log:

```text
results_g1_g8_h_r0_shadow_ref_20260621/logs/g1_g8_current_scheme_20260621_175441.log
```

## Current Conclusion

G1-G8 supports the current sensitivity gate: the active Branch A plus Branch B
V2 report-layer filter preserves all locked truth events and produces no FN in
this cohort. The next refinement should focus on reducing report/internal-review
burden, especially G2, without converting single auxiliary labels into hard
filters and without losing weak or mosaic truth events.
