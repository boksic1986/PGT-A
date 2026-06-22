# Handoff: Report-Table Ablation Audit

Date: 2026-06-22 22:35

## Context Entry

Use this handoff after:

```text
docs/handoff/2026-06-22_2145_plot_event_manifest_low_confidence_handoff.md
```

Current branch:

```text
codex/cnv-plot-wisecondor-style
```

Current report:

```text
docs/reports/branch_b_v2_report_table_ablation_audit_2026-06-22.md
```

## What Changed

This loop adds a formal workflow target for Branch B V2 report-table ablation
auditing:

```text
branch_b_v2_report_ablation
```

It produces:

```text
report_table_ablation_audit.tsv
report_table_ablation_summary.json
report_table_ablation_audit.md
```

The audit is development-only. It does not change Branch A, Branch B V2
classifier decisions, Branch S, reference build, sex calling, report-event
classification, or production filtering.

## Rule Under Audit

Only this class is proposed for report-table demotion:

```text
autosomal report_event
AND plot_support_class = Z_SUPPORTED_CN_NOT_SUPPORTED
=> downgrade_to_internal_review_candidate
```

Locked-truth overlap candidates remain protected. No hard filtering is
performed.

## Files Modified

Code/workflow:

- `Snakefile`
- `rules/pipeline_modes.smk`
- `rules/predict_layout.smk`
- `rules/predict_workflow.smk`
- `rules/script_entrypoints.smk`
- `rules/target_assembly.smk`
- `scripts/_compat_entry.py`
- `pgta/predict/branch_b/plot_manifest_audit.py`

Tests:

- `tests/unit/test_plot_manifest_ablation_audit.py`
- `tests/unit/test_branch_ab_phase12_workflow_contract.py`

Docs:

- `docs/reports/branch_b_v2_report_table_ablation_audit_2026-06-22.md`
- `docs/handoff/2026-06-22_2235_report_table_ablation_audit_handoff.md`
- `CURRENT_STATE.md`
- `PLANS.md`
- `docs/CURRENT_CONTEXT_INDEX.md`

## Remote Validation

Remote path:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419
```

Unit tests:

```text
73 passed in 4.25s
```

Command:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_plot_manifest_ablation_audit.py \
  tests/unit/test_branch_b_plot.py \
  tests/unit/test_cnv_report.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py \
  tests/unit/test_branch_s_shadow.py -q
```

Materialized target:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile <active gap2m lowres config> \
  --cores 1 branch_b_v2_report_ablation
```

Configs materialized:

- `config_predict_y_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`
- `config_predict_h_20260608_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`
- `config_predict_g1_g8_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`
- `config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`

## Materialized Summary

| cohort | original report events | proposed report events | proposed demotions | truth events | truth preserved | FN |
|---|---:|---:|---:|---:|---:|---:|
| Y1-Y8 | 40 | 40 | 0 | 10 | 10 | 0 |
| H1-H16 | 23 | 22 | 1 | 10 | 10 | 0 |
| G1-G8 | 26 | 26 | 0 | 10 | 10 | 0 |
| 0615 | 71 | 68 | 3 | 0 | 0 | 0 |

0615 is context-only; TP/FN/FP is not computed.

## Proposed Demotion Candidates

- H8 `chr4:10.50-48.75Mb gain`, `H8.A0004`.
- JZ26125846-61-61 `chr4:52.50-101.25Mb gain`,
  `JZ26125846-61-61.A0022`.
- JZ26125847-62-62 `chr10:49.50-135.00Mb gain`,
  `JZ26125847-62-62.A0019`.
- JZ26125847-62-62 `chr8:46.50-141.75Mb gain`,
  `JZ26125847-62-62.A0012`.

## Boundaries

- This is an ablation audit, not a production filter.
- `Z_AND_CN_NOT_SUPPORTED` is not automatically demoted in this loop.
- `CN_DIRECTION_WEAK_OR_MIXED` remains protected from automatic demotion.
- Branch S remains separate from autosomal CNV report-table demotion.
- Legacy Branch B fields are not used.

## Next Recommended Step

Review the four proposed demotion candidates manually against z/CN plots. If
the rule is promoted, implement it as a separate report-table demotion change
and rerun Y/H/G truth gates plus 0615 context burden.
