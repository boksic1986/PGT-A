# Branch B V2 Report-Table Ablation Audit

Date: 2026-06-22

## Scope

This audit evaluates whether a plot/support-derived low-confidence class can be
considered for report-table demotion. It does not change Branch A, Branch B V2
classification, Branch S, reference build, sex calling, or production report
filtering.

Current fixed inputs:

- reference ID: `h_r0_shadow_ref_20260619`
- Branch A overlay: explicit `merge_gap_bp=2_000_000`
- truth cohorts: Y1-Y8, H1-H6 within H1-H16, G1-G8
- context-only cohort: 2026-06-15

New workflow target:

```text
branch_b_v2_report_ablation
```

New outputs under each cohort result directory:

```text
wisecondorx/cnv/postprocess_gap2m/v2_benchmark/report_table_ablation_audit.tsv
wisecondorx/cnv/postprocess_gap2m/v2_benchmark/report_table_ablation_summary.json
wisecondorx/cnv/postprocess_gap2m/v2_benchmark/report_table_ablation_audit.md
```

## Candidate Rule Under Audit

Only this candidate rule is audited:

```text
autosomal report_event
AND plot_support_class = Z_SUPPORTED_CN_NOT_SUPPORTED
=> proposed_report_table_action = downgrade_to_internal_review_candidate
```

Guardrails:

- locked-truth overlap candidates are retained as report events;
- no event is hard filtered by this audit;
- `CN_DIRECTION_WEAK_OR_MIXED` is not demoted;
- 0615 remains context-only and is not used for TP/FN/FP;
- legacy Branch B fields such as `final_disposition` and
  `branch_b_keep_event` are not used.

Important limitation:

`Z_AND_CN_NOT_SUPPORTED` is not automatically demoted in this first audit. It
needs separate review because the current support table can reflect sparse or
unavailable display evidence, and truth labels may still be incomplete.

## Remote Validation

Remote path:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419
```

Executed tests:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_plot_manifest_ablation_audit.py \
  tests/unit/test_branch_b_plot.py \
  tests/unit/test_cnv_report.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py \
  tests/unit/test_branch_s_shadow.py -q
```

Result:

```text
73 passed in 4.25s
```

Remote materialization:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile <active gap2m lowres config> \
  --cores 1 branch_b_v2_report_ablation
```

Materialized configs:

- `config_predict_y_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`
- `config_predict_h_20260608_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`
- `config_predict_g1_g8_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`
- `config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`

## Cohort Summary

| cohort | status | original report events | final plot events | proposed report events after audit | proposed internal-review demotions | truth events | truth preserved | FN | hard-suppressed truth |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Y1-Y8 | ready | 40 | 40 | 40 | 0 | 10 | 10 | 0 | 0 |
| H1-H16 | ready | 23 | 22 | 22 | 1 | 10 | 10 | 0 | 0 |
| G1-G8 | ready | 26 | 26 | 26 | 0 | 10 | 10 | 0 | 0 |
| 0615 | context only | 71 | 68 | 68 | 3 | 0 | 0 | 0 | 0 |

0615 has no locked truth, so TP/FN/FP is not computed.

## Proposed Demotion Candidates

| cohort | sample | candidate | region | state | support class | same-direction ref-z median | CN same-direction fraction | truth guard |
|---|---|---|---|---|---|---:|---:|---|
| H1-H16 | H8 | H8.A0004 | chr4:10.50-48.75Mb | gain | Z_SUPPORTED_CN_NOT_SUPPORTED | 1.192013 | 0.0 | no_truth_overlap |
| 0615 | JZ26125846-61-61 | JZ26125846-61-61.A0022 | chr4:52.50-101.25Mb | gain | Z_SUPPORTED_CN_NOT_SUPPORTED | 0.831128 | 0.0 | context_only_no_truth |
| 0615 | JZ26125847-62-62 | JZ26125847-62-62.A0019 | chr10:49.50-135.00Mb | gain | Z_SUPPORTED_CN_NOT_SUPPORTED | 0.769887 | 0.0 | context_only_no_truth |
| 0615 | JZ26125847-62-62 | JZ26125847-62-62.A0012 | chr8:46.50-141.75Mb | gain | Z_SUPPORTED_CN_NOT_SUPPORTED | 0.554839 | 0.0 | context_only_no_truth |

## Interpretation

The audit supports this narrow conclusion:

```text
Z_SUPPORTED_CN_NOT_SUPPORTED can be treated as a report-table demotion
candidate for non-truth autosomal report events.
```

It does not yet support production filtering. If promoted later, the change
should be report-table demotion only:

```text
report_event -> internal_review_event_candidate
```

The evidence should remain in plot/support/audit outputs. It must not become a
hard filter without a separate locked-truth validation gate.

## Next Gate

1. Manually review the four proposed demotion candidates against z/CN plots.
2. Decide whether to implement report-table demotion as a separate code change.
3. If implemented, rerun Y/H/G truth gates and require:
   - Y1-Y8 truth 10/10, FN=0;
   - H1-H6 truth 10/10, FN=0, H6 chr21 visible;
   - G1-G8 truth 10/10, FN=0, G2 truth visible;
   - truth hard-filtered count = 0.
4. Keep 0615 burden/context only until locked truth labels are available.
