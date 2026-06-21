# Handoff: G1-G8 Current-Scheme Validation

Date: 2026-06-21 18:10 Asia/Shanghai

Status: `active_current_handoff`

## 1. Goal

Record the G1-G8 current-scheme validation result, update the active R&D
context files, and preserve the current interpretation before the next Branch B
V2 report-burden refinement loop.

## 2. Restored Context

Context was restored from:

- `docs/CURRENT_CONTEXT_INDEX.md`
- `docs/handoff/2026-06-21_1730_branch_b_v2_report_layer_filter_handoff.md`
- `AGENTS.md`
- `skills/conversation_handoff/SKILL.md`
- `skills/pgta_reference_modeling_analysis/SKILL.md`
- `CURRENT_STATE.md`
- `PLANS.md`

Active contract:

```text
reference_id=h_r0_shadow_ref_20260619
reference_status=fixed_shadow_baseline_not_production
branch_a_overlay=merge_gap_bp=2_000_000
strong_z=10.0
postprocess_path=wisecondorx/cnv/postprocess_gap2m
branch_b_v2_final_impact=development_review_only
branch_s_status=review_reportable_with_limitations
```

## 3. Completed Work

Added the G1-G8 validation report:

```text
docs/reports/branch_b_v2_g1_g8_validation_2026-06-21.md
```

Updated context files:

```text
docs/CURRENT_CONTEXT_INDEX.md
CURRENT_STATE.md
PLANS.md
```

Created this handoff:

```text
docs/handoff/2026-06-21_1810_g1_g8_current_scheme_validation_handoff.md
```

No CNV calling, Branch A, Branch B V2, Branch S, mosaic, sex-calling, or result
schema code was changed in this handoff loop.

## 4. Current Result

Remote result root:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419/results_g1_g8_h_r0_shadow_ref_20260621
```

G1-G8 result under the active current scheme:

```text
samples=8
candidates=75
locked_truth_events=10
truth_preserved=10
FN=0
truth_hard_suppressed=0
truth_report_layer_filtered=0
report_events=15
internal_review_events=40
filtered_audit_only_events=7
branch_s_events=13
```

Per-sample counts:

| Sample | Candidates | Truth | Preserved | FN | Report | Internal review | Filtered | Branch S |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| G1 | 4 | 1 | 1 | 0 | 1 | 3 | 0 | 0 |
| G2 | 28 | 2 | 2 | 0 | 6 | 16 | 3 | 3 |
| G3 | 13 | 2 | 2 | 0 | 1 | 10 | 0 | 2 |
| G4 | 4 | 1 | 1 | 0 | 2 | 2 | 0 | 0 |
| G5 | 8 | 1 | 1 | 0 | 1 | 3 | 2 | 2 |
| G6 | 8 | 1 | 1 | 0 | 2 | 3 | 0 | 3 |
| G7 | 4 | 1 | 1 | 0 | 1 | 1 | 2 | 0 |
| G8 | 6 | 1 | 1 | 0 | 1 | 2 | 0 | 3 |

Truth preservation summary:

- All 10 locked G1-G8 truth events are preserved.
- No truth event is hard suppressed.
- No truth event is report-layer filtered.
- chrX truth events are routed to Branch S, not autosomal CNV report events.
- G2 chr8 gain is preserved as `internal_review_event`, not filtered.

## 5. Interpretation

G1-G8 supports the current sensitivity/preservation gate for Branch A plus
Branch B V2 report-layer filtering.

This does not prove production readiness:

- The active reference remains a shadow baseline.
- Branch B V2 remains `development_review_only`.
- Branch S remains review/reportable-with-limitations, not final SCA.
- G2 still has high review/report burden and should be inspected before the
  next truth-safe demotion pass.

## 6. Key Files

Remote benchmark and report outputs:

```text
results_g1_g8_h_r0_shadow_ref_20260621/wisecondorx/cnv/postprocess_gap2m/v2_benchmark/summary.json
results_g1_g8_h_r0_shadow_ref_20260621/wisecondorx/cnv/postprocess_gap2m/v2_benchmark/sample_summary.tsv
results_g1_g8_h_r0_shadow_ref_20260621/wisecondorx/cnv/postprocess_gap2m/v2_benchmark/truth_metrics.tsv
results_g1_g8_h_r0_shadow_ref_20260621/wisecondorx/cnv/postprocess_gap2m/v2_benchmark/filtered_events.tsv
results_g1_g8_h_r0_shadow_ref_20260621/wisecondorx/cnv/report/cnv_summary.tsv
results_g1_g8_h_r0_shadow_ref_20260621/wisecondorx/cnv/report/cnv_summary.json
```

Runtime configs on the remote mirror:

```text
config_predict_g1_g8_h_r0_shadow_20260621.yaml
config_predict_g1_g8_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml
```

Truth table on the remote mirror:

```text
sample_info/g1_g8_truth_events_20260621.tsv
```

## 7. Commands

Dry-run:

```bash
cd /data/project/CNV/PGT-A/refactor_validation_20260419
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile config_predict_g1_g8_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml \
  --cores 1 -n branch_b_v2_benchmark branch_s_review cnv_report
```

Result:

```text
Passed. Planned 127 jobs for G1-G8, using existing BAM inputs.
```

Materialization:

```bash
cd /data/project/CNV/PGT-A/refactor_validation_20260419
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile config_predict_g1_g8_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml \
  --cores 32 branch_b_v2_benchmark branch_s_review cnv_report
```

Result:

```text
Completed. 127 of 127 steps (100%) done.
```

Log:

```text
results_g1_g8_h_r0_shadow_ref_20260621/logs/g1_g8_current_scheme_20260621_175441.log
```

## 8. Next Recommended Step

Inspect G2 and the current `report_event` rows across Y/H/G before designing a
second truth-safe demotion pass. The next demotion pass must preserve:

- Y1-Y8 truth 10/10;
- H1-H6 truth 10/10;
- G1-G8 truth 10/10;
- H6 chr21;
- no report-layer filtered truth.

H7-H16 and 2026-06-15 remain burden/context only unless locked truth labels are
added.

## 9. Core File Sync

- `docs/CURRENT_CONTEXT_INDEX.md`: updated with G1-G8 validation evidence.
- `CURRENT_STATE.md`: updated with G1-G8 materialized result and interpretation.
- `PLANS.md`: updated so future demotion rules must ablate against G1-G8 too.
- `AGENTS.md`: not updated; no repository-level hard rule changed.
