# Handoff: Branch B V2 Autosomal Burden Audit

Date: 2026-06-21 10:21 Asia/Shanghai

Status: `active_current_handoff`

## 1. Goal

Continue Branch B V2 evidence/disposition refinement after the sex-route
refinement by auditing the remaining autosomal `V2_POSITIVE_SUPPORT_REVIEW`
burden.

This loop did not change classifier logic. It fixed the next safe direction:
Branch B-side direction support may be useful as review evidence, but it is not
safe as a hard filter or universal positive-support downgrade.

## 2. Context Restored

Read before execution:

- `C:\Users\11217\.codex\attachments\b29e8565-0516-429a-91ff-991e2ad43c59\goal-objective.md`
- `docs/CURRENT_CONTEXT_INDEX.md`
- `docs/handoff/2026-06-21_1006_branch_b_v2_sex_route_handoff.md`
- `AGENTS.md`
- `skills/conversation_handoff/SKILL.md`
- `skills/pgta_reference_modeling_analysis/SKILL.md`
- `CURRENT_STATE.md`
- `PLANS.md`

Active reference:

```text
h_r0_shadow_ref_20260619
```

Active Branch A overlay:

```text
merge_gap_bp=2_000_000
strong_z=10.0
postprocess.output_dir=wisecondorx/cnv/postprocess_gap2m
```

## 3. Completed Work

Remote result audit:

- Parsed all gap2m `v2_classifier/*.candidate_classification.tsv` files for
  Y1-Y8, H1-H16, and 2026-06-15.
- Audited only:

```text
event_arm_class=autosome
v2_candidate_class=V2_POSITIVE_SUPPORT_REVIEW
```

- Parsed Y/H `v2_benchmark/truth_metrics.tsv` to test candidate-level truth
  safety.

New report:

- `docs/reports/branch_b_v2_autosomal_burden_audit_2026-06-21.md`

Core file updates:

- `CURRENT_STATE.md`
- `PLANS.md`
- `docs/CURRENT_CONTEXT_INDEX.md`

## 4. Current Evidence

Remaining autosomal positive-support rows:

| cohort | autosomal positive rows |
|---|---:|
| Y1-Y8 truth cohort | 78 |
| H1-H6 truth subset | 33 |
| H7-H16 context subset | 24 |
| 2026-06-15 context | 127 |

Dominant unresolved condition:

```text
UNKNOWN_BACKGROUND + NO_NULL_SUPPORT
```

All audited autosomal positive rows share:

```text
matched_negative_background_status=UNKNOWN_BACKGROUND
calibration_null_status=NO_NULL_SUPPORT
refmap_status=OK
sample_noise_status=OK
cnvpro_like_evidence_status=SHADOW_EVIDENCE_ONLY
```

## 5. Direction-Support Safety Finding

A candidate split was audited:

```text
B_DIRECTION_SUPPORTED:
  same_direction_fraction >= 0.5
  OR corrected/raw amplitude supports the same event direction with abs >= 2

A_ONLY_NO_B_DIRECTION_SUPPORT:
  Branch A z/support exists, but Branch B-side direction evidence is weak
```

Impact:

| cohort | autosomal positive rows | B-direction supported | A-only / weak B-direction |
|---|---:|---:|---:|
| Y1-Y8 truth cohort | 78 | 58 | 20 |
| H1-H16 all | 57 | 35 | 22 |
| H7-H16 context subset | 24 | 15 | 9 |
| 2026-06-15 context | 127 | 89 | 38 |

Truth safety blocker:

Several locked autosomal truth top candidates would be labeled
`A_ONLY_NO_B_DIRECTION_SUPPORT`:

- Y2 chr14 gain;
- Y4 chr13 gain;
- H2 chr6 gain;
- H3 chr13 gain;
- H4 chr15 gain.

H6 chr21 remains direction-supported, but this does not make the rule generally
safe.

## 6. Current Conclusion

Do not implement Branch B-side direction support as:

- hard suppression;
- final benign/artifact evidence;
- universal demotion out of positive-support review.

It can be considered in a later loop as a review label or report-visible
evidence tag, but only if truth preservation remains explicit and no final
decision is made from that tag alone.

## 7. Suggested Next Step

The next safe implementation target is a review-label-only output, not a hard
filter:

1. Preserve `V2_POSITIVE_SUPPORT_REVIEW` for truth-safe preservation.
2. Add a non-final review label or summary field distinguishing
   `B_DIRECTION_SUPPORTED` from `A_ONLY_WEAK_B_DIRECTION`.
3. Keep `final_report_impact=none_shadow_only`.
4. Verify Y1-Y8/H1-H6 truth preservation, H6 chr21, and no hard suppression.

If adding that label would require a schema change, write the change as an
explicit workflow contract update first.

## 8. Key Files

- `docs/reports/branch_b_v2_autosomal_burden_audit_2026-06-21.md`
- `docs/reports/branch_b_v2_sex_route_refinement_2026-06-21.md`
- `CURRENT_STATE.md`
- `PLANS.md`
- `docs/CURRENT_CONTEXT_INDEX.md`

Remote result roots:

- `/data/project/CNV/PGT-A/refactor_validation_20260419/results_y_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m`
- `/data/project/CNV/PGT-A/refactor_validation_20260419/results_h_20260608_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m`
- `/data/project/CNV/PGT-A/refactor_validation_20260419/results_20260615_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess_gap2m`

## 9. Command Record

Remote executable:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/python
```

Remote command purpose:

```text
Parse materialized v2_classifier TSVs and truth_metrics TSVs to aggregate
autosomal V2 positive-support burden and direction-support safety.
```

Result:

```text
Completed; outputs summarized in
docs/reports/branch_b_v2_autosomal_burden_audit_2026-06-21.md.
```

## 10. Core File Sync

- `CURRENT_STATE.md`: updated with autosomal burden audit and direction-support
  safety finding.
- `PLANS.md`: updated to block direction-support hard filtering and keep it as
  review-label-only unless validated.
- `docs/CURRENT_CONTEXT_INDEX.md`: updated to include this audit report and the
  next-gate constraint.
- `AGENTS.md`: not updated; hard constraints did not change.
- `REPO_MAP.md`: not updated; stable repository structure did not change.

## 11. Do Not Misread

- Do not treat this audit as a classifier improvement.
- Do not hard-filter candidates with weak Branch B direction support.
- Do not treat 2026-06-15 as TP/FN/FP evidence.
- Do not treat H7-H16 as locked truth.
- Do not promote Branch B V2, Branch S, the gap2m overlay, or the shadow
  reference from this audit.
