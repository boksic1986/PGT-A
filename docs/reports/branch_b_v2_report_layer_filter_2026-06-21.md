# Branch B V2 Report-Layer Filtering

Date: 2026-06-21

## Scope

This loop makes Branch B V2 perform report-layer filtering on top of the
existing Branch A gap2m candidate set. It does not modify WisecondorX predict,
Branch A candidate generation, Branch S logic, or the active shadow reference.

Active contract:

- reference: `h_r0_shadow_ref_20260619`
- Branch A overlay: `merge_gap_bp=2_000_000`
- strong/sensitive split: `strong_z=10.0`
- Branch B V2 impact: `development_review_only`

## Report-Layer Classes

Branch B V2 now emits four report-layer classes:

- `report_event`: autosomal CNV event eligible for the development report main
  event count.
- `internal_review_event`: retained for internal review and truth-safety, but
  not counted as a final report event.
- `filtered_event`: excluded from report and internal-review main flow; retained
  in `filtered_events.tsv/json` audit ledgers.
- `branch_s_event`: chrX/chrY route for Branch S review sections.

## Initial Rules

`report_event` requires all of:

- autosomal chromosome;
- `a_abs_zscore >= 10`;
- `v2_b_signal_context_label=B_SIGNAL_SUPPORTED_A_DIRECTION`;
- length tier in `large_ge10mb`, `broad_review_ge4mb`, or
  `reportable_candidate_ge2mb`;
- `v2_clean_support_label=CLEAN_SUPPORT_AVAILABLE`.

`filtered_event` can be assigned only by combination evidence. The current
non-contract technical filter requires all of:

- not strong A signal (`a_abs_zscore < 10`);
- B-side signal is not supportive of Branch A direction;
- GC/RC context is attenuated or severely attenuated;
- at least one of sensitive A signal, short/focal length, or low-clean/high-risk
  support.

Single factors are not sufficient filters:

- background unknown;
- B-side weak/discordant signal;
- GC/RC attenuation;
- short length;
- `a_abs_zscore < 10`;
- blacklist/hard-region overlap.

Strong A signals with B/GC disagreement are retained as
`internal_review_event`, not filtered. H6 chr21 remains protected and retained
as internal review.

## Remote Results

Remote path:
`/data/project/CNV/PGT-A/refactor_validation_20260419`

### Y1-Y8

| Metric | Value |
|---|---:|
| Branch A gap2m candidates | 97 |
| Truth events | 10 |
| Truth preserved | 10 |
| FN | 0 |
| Hard-suppressed truth | 0 |
| Report-layer filtered truth | 0 |
| V2 report events | 21 |
| V2 internal review events | 50 |
| V2 filtered events | 13 |
| V2 Branch S events | 13 |

### H1-H16

| Metric | Value |
|---|---:|
| Branch A gap2m candidates | 105 |
| Truth events | 10 |
| Truth preserved | 10 |
| FN | 0 |
| Hard-suppressed truth | 0 |
| Report-layer filtered truth | 0 |
| V2 report events | 6 |
| V2 internal review events | 49 |
| V2 filtered events | 8 |
| V2 Branch S events | 42 |

H6 chr21 gain remains visible as `internal_review_event` with
`a_abs_zscore=7.113507302991461`.

### 2026-06-15 Context Cohort

This cohort has no locked truth and remains burden/context only.

| Sample | Candidates | Report events | Internal review | Filtered | Branch S |
|---|---:|---:|---:|---:|---:|
| JZ26125843-56-56 | 26 | 3 | 13 | 8 | 2 |
| JZ26125844-59-59 | 35 | 9 | 19 | 6 | 1 |
| JZ26125845-60-60 | 47 | 22 | 20 | 0 | 5 |
| JZ26125846-61-61 | 33 | 14 | 15 | 1 | 3 |
| JZ26125847-62-62 | 24 | 4 | 9 | 8 | 3 |
| Total | 165 | 52 | 76 | 23 | 14 |

## Report Contract

`cnv_report` now prefers Branch B V2 report-layer counts in
`biological_candidate_conclusion` when V2 summaries are present. Legacy/current
code Branch B top events remain available as a labeled legacy context column,
but they are not the main biological conclusion.

Report output fields include:

- `branch_b_v2_report_event_count`
- `branch_b_v2_internal_review_event_count`
- `branch_b_v2_filtered_event_count`
- `branch_b_v2_branch_s_event_count`
- `branch_b_v2_legacy_fields_used=false`
- `branch_b_v2_final_impact=development_review_only`

## Validation

Remote unit tests:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_v2_classifier.py \
  tests/unit/test_branch_b_v2_benchmark.py \
  tests/unit/test_cnv_report.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py -q

66 passed in 1.15s
```

Remote Snakemake dry-run:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile <active_gap2m_config> --cores 1 -n \
  branch_b_v2_benchmark branch_s_review cnv_report

Passed for all three active gap2m configs.
```

Remote materialization:

```text
Forced cnv_branch_b_v2_classifier, cnv_branch_b_v2_benchmark,
cnv_report_summary, branch_b_v2_benchmark, and cnv_report.

Y1-Y8: 14/14 steps done.
H1-H16: 22/22 steps done.
2026-06-15: 11/11 steps done.
```

## Interpretation

This is the first Branch B V2 loop that performs actual report-layer filtering:
filtered events are removed from report/internal-review main flow and retained
only in audit ledgers.

It is still not production promotion:

- 2026-06-15 has no truth, so its burden reduction cannot prove FP precision.
- Branch S remains review/reportable-with-limitations, not final SCA.
- `h_r0_shadow_ref_20260619` remains a fixed shadow baseline, not production.
- Internal review events are still numerous, especially in 2026-06-15.

The next step should inspect V2 `report_event` rows, especially high-burden
samples such as JZ26125845-60-60, and decide whether additional non-overfit
combination rules can move more events from report to internal review or audit
without harming Y/H locked truth and H6 chr21.
