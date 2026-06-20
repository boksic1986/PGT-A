# P1-P6 Result Credibility Audit

Date: 2026-06-21

Status: `active_current_evidence`

Decision use: current credibility boundary. Pair with
`docs/CURRENT_CONTEXT_INDEX.md` before using any P1-P6 result for planning.

## Purpose

This audit re-checks the current G1/P1-P6 result boundary after the user raised
concerns that legacy/current-code Branch B thresholds are not reliable evidence.

This audit does not change thresholds, CNV calling, mosaic logic, sex calling,
or report schema. It classifies which existing phase outputs can be used as
current R&D evidence and which outputs must be treated as historical, shadow, or
engineering-only artifacts.

## Context Read

- `docs/handoff/2026-06-20_2025_p6_report_package_handoff.md`
- `AGENTS.md`
- `skills/conversation_handoff/SKILL.md`
- `skills/pgta_reference_modeling_analysis/SKILL.md`
- `CURRENT_STATE.md`
- `PLANS.md`
- `docs/reports/branch_ab_v2_gate_phase_plan_2026-06-20.md`
- `docs/reports/branch_ab_v2_rnd_constraints_2026-06-18.md`
- `docs/reports/h7_h16_reference_candidate_audit_2026-06-20.md`
- `docs/reports/branch_a_validation_h_r0_shadow_2026-06-20.md`
- `docs/reports/branch_b_p3_evidence_contract_2026-06-20.md`
- `docs/reports/branch_s_p5_report_boundary_2026-06-20.md`
- `docs/reports/p6_report_package_contract_2026-06-20.md`
- `pgta/reference/audit.py`
- `pgta/predict/branch_b/evidence_ledger.py`
- `pgta/predict/branch_b/v2_classifier.py`
- `pgta/predict/branch_s.py`
- `pgta/predict/report.py`

Active handoff:

```text
docs/handoff/2026-06-20_2025_p6_report_package_handoff.md
```

## Credibility Classes

| class | meaning |
|---|---|
| `VALID_CURRENT_EVIDENCE` | Can be used as current R&D evidence for the stated scope. |
| `VALID_ENGINEERING_ARTIFACT_ONLY` | Workflow/output contract is useful, but it cannot support biological or performance claims. |
| `HISTORICAL_OR_INVALID_FOR_DECISION` | Historical or shadow output only; must not drive current decisions. |

## Phase Classification

| phase | result path | current status | can support what | cannot support what |
|---|---|---|---|---|
| G1/P1 reference audit | `results_h_20260608_h_r0_shadow_ref_20260619/wisecondorx/cnv/reference_audit/reference_candidate_audit.summary.json` | `VALID_CURRENT_EVIDENCE` for R-label audit; reference remains shadow | H7-H16 R0/R1/R2 reference-rebuild eligibility under `h_r0_shadow_ref_20260619`; old Branch B kept counts were not used | Production reference promotion; N0/N1/N2 background assignment; Branch B V2 performance |
| P2 Branch A validation | `results_y_h_r0_shadow_ref_20260619/wisecondorx/cnv/branch_a_validation/summary.json`; `results_h_20260608_h_r0_shadow_ref_20260619/wisecondorx/cnv/branch_a_validation/summary.json` | `VALID_CURRENT_EVIDENCE` for current shadow-reference Branch A sensitivity | Y1-Y8 and H1-H6 truth are detected under the same reference ID; current known-positive FN evidence is zero | Final reference promotion; Branch B V2 FP reduction; SCA final status |
| P2 0615 exploratory burden | `results_20260615_h_r0_shadow_ref_20260619/wisecondorx/cnv/branch_a_validation/summary.json` | `VALID_ENGINEERING_ARTIFACT_ONLY` | Candidate burden under current shadow reference | TP/FN/FP, because truth_event_count is zero |
| P3 Branch B evidence ledger | `results_*_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess/evidence_ledger/*.summary.json` | `VALID_ENGINEERING_ARTIFACT_ONLY` | One-row-per-Branch-A-candidate evidence ledger and review-safe report impact | Branch B V2 performance; hard PASS/artifact calls; final FP reduction |
| current-code Branch B final/artifact fields | `final_disposition`, `branch_b_keep_event`, `branch_b_artifact_status` in P3 ledger | `HISTORICAL_OR_INVALID_FOR_DECISION` for performance | Historical/shadow context only if explicitly labeled | V2 evidence, V2 benchmark, N0/N1, report release, truth comparison |
| Branch B V2 classifier shadow output | `results_*_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess/v2_classifier/*.summary.json` | `HISTORICAL_OR_INVALID_FOR_DECISION` for current V2 performance | Shows that a shadow classifier path exists and is `none_shadow_only` | V2-only performance; FP/FN claims; final report disposition, because classifier still uses legacy `final_disposition` branches |
| P5 Branch S | `results_*_h_r0_shadow_ref_20260619/wisecondorx/cnv/postprocess/branch_s/*.summary.json` | `VALID_ENGINEERING_ARTIFACT_ONLY` | SCA review/development output boundary and summary schema | Final SCA calls, sex-calling replacement, report release |
| P6 report package | `results_20260615_h_r0_shadow_ref_20260619/wisecondorx/cnv/report/cnv_summary.json` | `VALID_ENGINEERING_ARTIFACT_ONLY` | Development-only report package carrying P3/P5 summary status | Final 0615 report release; locked performance; Branch B/S promotion |

## Static Data-Flow Findings

Local static command:

```powershell
rg -n "final_disposition|branch_b_keep_event|legacy|current-code|v2_candidate_class|report_impact|review_development_only" pgta rules tests docs
```

Key findings:

- `pgta/reference/audit.py` writes
  `legacy_branch_b_kept_counts_used_for_r_label=false`; G1 does not use old
  Branch B kept counts for R-label assignment.
- `pgta/predict/branch_b/evidence_ledger.py` retains
  `final_disposition`, `branch_b_keep_event`, and related current-code Branch B
  fields in the ledger, but P3 disposition is forced to review-safe output and
  `report_impact=none_shadow_only`.
- `pgta/predict/branch_b/v2_classifier.py` still reads `final_disposition` and
  produces `CONFIRMED_SHADOW`, `MOSAIC_SUSPECT_SHADOW`, or
  `LIKELY_ARTIFACT_SHADOW`. Therefore current v2 classifier output is not
  V2-only performance evidence.
- `pgta/predict/branch_s.py` keeps Branch S output as
  `sca_output_mode=review_development_only` and `final_report_impact=none_shadow_only`.
- `pgta/predict/report.py` consumes P3/P5 summary JSONs only. It does not
  consume raw P3 ledger TSV, raw Branch S evidence TSV, Branch S score TSV,
  matched-negative percentile TSV, or Branch B V2 classifier TSV.

## Remote Result Audit

All successful remote checks used:

```text
executable: ssh fengxian
remote cwd: /data/project/CNV/PGT-A/refactor_validation_20260419
python: /biosoftware/miniconda/envs/snakemake_env/bin/python
```

### Summary JSON key fields

Remote result:

```text
G1 reference audit:
reference_id=h_r0_shadow_ref_20260619
sample_count=10
label_counts={'R0': 6, 'R1': 4}
legacy_branch_b_kept_counts_used_for_r_label=False
formal_n0_used=False

P2 Y1-Y8 Branch A:
reference_id=h_r0_shadow_ref_20260619
sample_count=8
truth_event_count=10
truth_detected_count=10

P2 H1-H16 Branch A:
reference_id=h_r0_shadow_ref_20260619
sample_count=16
truth_event_count=10
truth_detected_count=10

P2 2026-06-15 exploratory:
reference_id=h_r0_shadow_ref_20260619
sample_count=5
truth_event_count=0
truth_detected_count=0

P6 2026-06-15 report:
report_contract.status=development_only_not_final_release
reference_id=h_r0_shadow_ref_20260619
branch_a_no_fn_status=passed_no_fn
same_reference_config_status=matched
branch_b_evidence_summary_count=5
branch_s_summary_count=5
branch_b_raw_ledger_used=False
branch_s_raw_evidence_used=False
```

### P3/P5/V2 shadow summary counts

Remote result:

```text
Y_P3 files=8 candidate_count=131 impact={'none_shadow_only': 8}
H_P3 files=16 candidate_count=221 impact={'none_shadow_only': 16}
0615_P3 files=5 candidate_count=201 impact={'none_shadow_only': 5} backgrounds={'UNKNOWN_BACKGROUND': 201}

Y_S files=8 impact={'none_shadow_only': 8} sca_modes={'review_development_only': 8}
H_S files=16 impact={'none_shadow_only': 16} sca_modes={'review_development_only': 16}
0615_S files=5 impact={'none_shadow_only': 5} sca_modes={'review_development_only': 5}

Y_V2 files=8 candidate_count=131 impact={'none_shadow_only': 8} classes={'LIKELY_ARTIFACT_SHADOW': 38, 'REVIEW_REQUIRED': 93}
H_V2 files=16 candidate_count=221 impact={'none_shadow_only': 16} classes={'LIKELY_ARTIFACT_SHADOW': 145, 'REVIEW_REQUIRED': 76}
0615_V2 files=5 candidate_count=201 impact={'none_shadow_only': 5} classes={'REVIEW_REQUIRED': 169, 'LIKELY_ARTIFACT_SHADOW': 32}
```

Interpretation:

- P3 is materialized and review-safe, but it is not a V2 performance result.
- P5 is materialized and review/development-only.
- V2 classifier summary exists, but it includes legacy-derived shadow classes
  and must be excluded from V2-only performance claims.

## Remote Dry-Run Status

Command attempted:

```bash
cd /data/project/CNV/PGT-A/refactor_validation_20260419 && \
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake \
  -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile \
  --configfile config_predict_y_h_r0_shadow_branch_b_evidence_20260620.yaml \
  --cores 1 -n branch_b_evidence branch_s_review cnv_report
```

Result:

```text
FAILED
IndentationError in file /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile, line 9:
expected an indented block after function definition on line 10 (Snakefile, line 9)
```

Follow-up read-only inspection:

```bash
nl -ba Snakefile | sed -n '1,40p'
```

showed the first 40 lines are visually consistent with the local tracked
`Snakefile`, and local `git diff -- Snakefile` is empty.

Interpretation:

- Current remote result files are inspectable and support the audit above.
- The requested remote Snakemake dry-run is not currently passable from the
  active remote mirror and must be fixed or re-synced before any new workflow
  validation claim is made.
- This dry-run failure does not invalidate already materialized result files,
  but it blocks claiming the current remote mirror is workflow-validated today.

## Current Valid Conclusions

1. G1/P1 correctly excludes legacy Branch B kept counts from R0/R1/R2
   reference audit labels.
2. P2 Branch A no-FN evidence is the only current performance-relevant evidence
   for known positives under `h_r0_shadow_ref_20260619`.
3. P3 Branch B evidence ledger is a workflow/evidence contract, not Branch B V2
   performance validation.
4. Current Branch B V2 classifier shadow output must not be used as V2
   performance evidence because it still reads legacy `final_disposition`.
5. P5 Branch S is a report-boundary artifact only; SCA final promotion has not
   passed.
6. P6 0615 report package is development-only and not final-reportable.
7. 2026-06-15 has no truth table and cannot be used for TP/FN/FP performance.

## Current Invalid Or Excluded Conclusions

The project must not claim:

- H-augmented reference is production-final.
- Branch B V2 has improved final-level reporting over legacy/current-code B.
- Current V2 shadow classifier proves FP reduction.
- Legacy/current-code B artifact/kept labels are valid evidence.
- N1/reference-cohort background is formal N0.
- Branch S is final SCA.
- 0615 report package is releasable as a forward clinical/report result.

## Next Gate

The next real gate should be:

```text
G2/P2 evidence review -> V2-only Branch B truth benchmark design
```

Required properties for the next gate:

- Use Y1-Y8 and H1-H6 truth only for initial performance.
- Compare against Branch A truth-overlap candidates, not legacy Branch B kept
  events.
- Exclude `final_disposition`, `branch_b_keep_event`, and legacy artifact status
  from V2 decision and metric computation.
- Keep H7-H16 and 0615 as burden/context only unless a locked truth label exists.
- Keep Branch S separate and review/development-only until locked SCA truth
  coverage passes.
