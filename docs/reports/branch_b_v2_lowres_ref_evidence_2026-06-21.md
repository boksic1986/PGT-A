---
status: interface_validated_not_materialized
decision_use: branch_b_v2_auxiliary_evidence_contract
date: 2026-06-21
active_reference_id: h_r0_shadow_ref_20260619
branch_a_overlay: merge_gap_bp_2000000
final_impact: development_review_only
---

# Branch B V2 Low-Resolution Ref Evidence

## Scope

This loop adds an optional Branch B V2 auxiliary evidence layer for
low-resolution reference checks and reference-bin stability. It does not build
or promote a new reference by itself.

Current primary discovery remains the active Branch A/WisecondorX/CBS candidate
chain under `h_r0_shadow_ref_20260619` and the explicit
`merge_gap_bp=2_000_000` overlay. Existing documentation records the active
H-R0 shadow reference best binsize as `750000`; this low-resolution layer does
not assume the primary candidate binsize is exactly 1 Mb.

## Method Contract

The low-resolution evidence layer answers one narrow question:

```text
Does the current Branch A candidate show same-direction support in a more
stable lower-resolution result set?
```

It must not:

- replace WisecondorX predict/CBS as the primary caller;
- create B-only events;
- change Branch A inclusion;
- treat low-resolution absence as a single hard filter;
- suppress H6 chr21-like short or weak positives;
- mix chrX/chrY/SCA candidates into autosomal low-resolution filtering.

## Implemented Files

Code:

- `pgta/predict/branch_b/lowres_evidence.py`
- `pgta/predict/branch_b/ref_stability.py`
- `pgta/predict/branch_b/v2_classifier.py`
- `scripts/_compat_entry.py`

Workflow:

- `rules/script_entrypoints.smk`
- `rules/predict_layout.smk`
- `rules/predict_workflow.smk`
- `rules/target_assembly.smk`

Tests:

- `tests/unit/test_lowres_evidence.py`
- `tests/unit/test_ref_stability.py`
- `tests/unit/test_branch_b_v2_classifier.py`
- `tests/unit/test_branch_ab_phase12_workflow_contract.py`

## Low-Resolution Evidence Output

When `core.wisecondorx.cnv.lowres_evidence.enable=true`, the workflow can
generate:

```text
wisecondorx/cnv/postprocess_gap2m/lowres_evidence/{sample}.candidate_evidence.tsv
wisecondorx/cnv/postprocess_gap2m/lowres_evidence/{sample}.summary.json
```

Per-candidate fields:

```text
lowres_2mb_support_label
lowres_2mb_same_direction
lowres_2mb_overlap_fraction
lowres_2mb_z_or_score
lowres_3mb_support_label
lowres_3mb_same_direction
lowres_3mb_overlap_fraction
lowres_3mb_z_or_score
lowres_consensus_label
```

Current label semantics:

| label | meaning | filter impact |
|---|---|---|
| `LOWRES_2MB_3MB_SAME_DIRECTION_SUPPORT` | 2Mb and 3Mb both support the Branch A direction | confidence context only |
| `LOWRES_2MB_SUPPORT_FOR_3_4MB_EVENT` | 2Mb supports a 3-4Mb candidate where 3Mb may be boundary-diluted | confidence context only |
| `LOWRES_2MB_SAME_DIRECTION_SUPPORT` | 2Mb supports a candidate | confidence context only |
| `LOWRES_3MB_SAME_DIRECTION_SUPPORT` | 3Mb supports a candidate | confidence context only |
| `LOWRES_NO_SUPPORT_INFORMATIVE_BUT_NOT_FILTER` | low-res was informative but lacked same-direction support | review context only; not a single filter |
| `LOWRES_CONTEXT_ONLY_SHORT_OR_BOUNDARY_EVENT` | candidate too short or boundary-sensitive for low-res denial | context only |
| `LOWRES_NOT_APPLICABLE_SEX_CHROM` | chrX/chrY candidate routed outside autosomal low-res logic | Branch S context |
| `LOWRES_NOT_CONFIGURED` | low-res inputs were not configured | no evidence |

## Ref-MAD Stability Output

When low-res evidence is enabled, reference-bin stability is computed from the
configured reference NPZ set:

```text
wisecondorx/cnv/postprocess_gap2m/ref_stability/ref_bin_stability.tsv
wisecondorx/cnv/postprocess_gap2m/ref_stability/summary.json
wisecondorx/cnv/postprocess_gap2m/ref_stability/{sample}.candidate_ref_stability.tsv
wisecondorx/cnv/postprocess_gap2m/ref_stability/{sample}.summary.json
```

Per-bin fields:

```text
ref_sample_count
ref_median
ref_std
ref_cv
ref_mad
ref_mad_z
ref_stability_label
```

Per-candidate fields:

```text
event_ref_mad_median
event_ref_mad_p90
high_ref_mad_bin_fraction
ref_stability_context
```

Ref-MAD interpretation:

- `REF_STABILITY_STABLE`: current ref bins look stable in this event region.
- `REF_STABILITY_MODERATE_REVIEW`: some bin instability; low-res absence should
  be interpreted cautiously.
- `REF_STABILITY_HIGH_MAD_REVIEW`: region is unstable in the reference itself;
  low-res absence cannot be strong negative evidence.
- `REF_STABILITY_NO_BIN_CONTEXT`: no overlapping ref bins were available.

## Workflow Contract

The new path is config-gated:

```yaml
core:
  wisecondorx:
    cnv:
      lowres_evidence:
        enable: true
        reference_npz:
          - /path/to/ref_sample_1.npz
          - /path/to/ref_sample_2.npz
        reference_sample_ids:
          - R1
          - R2
        events_2mb_dir: /path/to/2mb/events
        events_3mb_dir: /path/to/3mb/events
```

`reference_npz` is required when `enable=true`, because ref-MAD stability is
part of the low-res evidence contract. If `reference_sample_ids` is provided,
its length must match `reference_npz`.

Default behavior is unchanged. With `enable=false`, `cnv_branch_b_v2_classifier`
continues to consume the previous Branch B evidence or matched-negative ledger.

## Branch B V2 Classifier Integration

The classifier now emits:

```text
v2_lowres_context_label
v2_lowres_context_reason
v2_ref_stability_context
```

and includes optional audit tags:

```text
[lowres-shadow] consensus=...
[ref-mad-shadow] context=...
```

These fields do not change:

```text
v2_report_layer_class
v2_report_visibility
v2_filter_action
v2_burden_reduction_action
```

This explicitly preserves the rule that low-resolution absence is context, not
a standalone filtering rule.

## Validation Status

Remote unit tests and default-path dry-runs were performed on `ssh fengxian` in:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419
```

Confirmed unit tests:

```text
tests/unit/test_lowres_evidence.py
tests/unit/test_ref_stability.py
tests/unit/test_branch_b_v2_classifier.py
tests/unit/test_branch_b_v2_benchmark.py
tests/unit/test_branch_ab_phase12_workflow_contract.py
tests/unit/test_cnv_report.py
tests/unit/test_current_context_index.py
```

Result:

```text
81 passed in 1.69s
```

Default-path dry-runs passed for `branch_b_v2_benchmark branch_s_review
cnv_report` under all four active gap2m configs:

```text
config_predict_y_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml
config_predict_h_20260608_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml
config_predict_g1_g8_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml
config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_benchmark_20260621.yaml
```

All four returned `Nothing to be done`, confirming the default path remains
parseable and does not request low-res outputs unless low-res evidence is
enabled.

A temporary remote dry-run config also enabled low-res evidence using two
existing reference NPZ inputs. That dry-run planned the expected DAG:

```text
cnv_branch_b_ref_stability_bins: 1
cnv_branch_b_ref_stability: 8
cnv_branch_b_lowres_evidence: 8
cnv_branch_b_v2_classifier: 8
cnv_branch_b_v2_benchmark: 1
branch_b_v2_benchmark: 1
total: 27
```

This validated rule wiring only. It did not materialize low-res outputs.

## Not Yet Done

The following plan items are not completed in this loop:

- build `h_r0_shadow_ref_20260619_2mb`;
- build `h_r0_shadow_ref_20260619_3mb`;
- rerun low-res WisecondorX predict for Y1-Y8, H1-H16, G1-G8, and 2026-06-15;
- materialize low-res event directories for actual Branch B V2 evidence;
- compute truth-level 2Mb/3Mb same-direction support rates;
- use low-res/ref-MAD evidence for a candidate-level demotion ablation.

These require a separate long-task run and should be started in the remote
mirror with explicit PID/log/target reporting.

## Next Gate

Create named low-res shadow configs and launch the 2Mb/3Mb reference/predict
jobs as long tasks. After materialization, run:

- Y1-Y8 truth preservation;
- H1-H6 truth preservation with H6 chr21 visible;
- G1-G8 truth preservation;
- `>=4Mb` truth same-direction support rate at 2Mb and 3Mb;
- 2026-06-15 burden/context only;
- Branch S burden table separately for chrX/chrY/SCA.

Until that gate passes, low-res evidence is an implemented interface and
workflow contract, not a validated biological filtering layer.
