# Branch A/B/S R&D Constraints

## Purpose

This file is the constraints-only record for the current Branch A/B/S R&D
cycle.

It intentionally removes the older incorrect route that treated the provisional
Branch B V2 / N1 shadow outputs as the next report line. The execution order is
defined by:

```text
docs/reports/branch_ab_v2_gate_phase_plan_2026-06-20.md
```

Historical phase notes, legacy Branch B kept counts, old N0 seed status, and
shadow-report outputs must not be used as current execution plans.

## Context Source

Current controlling documents:

- `docs/reports/branch_ab_v2_gate_phase_plan_2026-06-20.md`
- `docs/reports/branch_b_v2_reference_background_and_sca_design_2026-06-20.md`
- `docs/reports/branch_s_sca_validation_plan_2026-06-19.md`
- `docs/handoff/2026-06-20_1523_branch_ab_plan_correction_handoff.md`

If these conflict with older notes, use the current result files and the
gate/phase plan first.

## Hard Constraints

1. WisecondorX predict/CBS remains the primary CNV caller.
2. Branch A is the sensitivity gate and must be derived from WisecondorX
   predict/CBS evidence.
3. Branch B refines Branch A candidates only.
4. Branch B must not create B-only final report events.
5. Branch B may add LOH, UPD, mosaic, CN-like amplitude, region-risk,
   sample-noise, background, and CNVpro-inspired consistency evidence.
6. Hard filtering must be narrow and must not increase known-positive FN.
7. CNVpro-inspired analysis is a reference/comparison direction, not a
   replacement caller.
8. Branch S/SCA must be developed into a reportable workflow, but current
   unvalidated Branch S output remains review/development-only.
9. All real validation, Snakemake dry-runs, and workflow tests must run on
   `ssh fengxian`.

## Branch A Constraints

Branch A is allowed to:

- expose WisecondorX-derived CNV candidate signals;
- merge adjacent same-direction WisecondorX/CBS signals when this preserves
  known-positive recall;
- report candidate burden and top A signals for reference audit, positives,
  negatives, and shadow-report samples.

Branch A must not:

- include Branch B artifact filters;
- use Branch B thresholds to suppress candidates;
- use legacy Branch B kept/review status as truth;
- be frozen until the named reference/config contract is fixed;
- be promoted if Y1-Y8 or H1-H6 known-positive recall regresses, especially the
  weak H6 chr21 event.

## Reference-Rebuild Constraints

H7-H16 are new-batch negative samples and reference-rebuild candidates, but they
are not automatically clean reference or N0 samples.

Reference expansion must use explicit `R0/R1/R2` labels:

- `R0`: candidate reference sample with acceptable QC and no post-rebuild
  evidence requiring exclusion.
- `R1`: shadow/reference-ablation candidate; can be tested but not silently
  promoted.
- `R2`: holdout or exclusion candidate pending review.

Reference-rebuild limitations:

- do not reject H7-H16 solely because the old 32-sample reference produced
  Branch A calls;
- do not promote the current H-augmented reference as final without whole-cohort
  abnormal-sample audit;
- do not use old-reference Branch B kept counts alone to define `R0/R1/R2`;
- every named reference must record cohort, binsize, mask/preprocess contract,
  WisecondorX command, and reference ID.

New reference IDs have two rerun levels:

- audit level: rerun WisecondorX predict and Branch A candidates; current
  Branch B/S/report outputs may exist only as diagnostic snapshots;
- report-facing level: rerun fixed Branch B, fixed Branch S, evaluation,
  benchmark, and report under the same reference/config contract.

A new reference ID does not by itself require rerunning `fastp_bwa`. If binsize,
mapping, read-filtering, mask, or preprocessing changes, regenerate the affected
NPZ/reference/predict products from valid existing BAM or from upstream data as
appropriate.

## N0/N1/N2 Constraints

`R0/R1/R2` and `N0/N1/N2` are separate labels.

- `R0/R1/R2` controls reference-rebuild eligibility.
- `N0/N1/N2` controls Branch B matched-negative or background use.

Formal `N0` / cross-fit background remains the clean promotion path, but it is
not an immediate blocker for candidate evidence development.

Until a reliable N0/cross-fit set exists:

- background-derived evidence must be labeled limited, reference-context, N1
  shadow, or unknown;
- missing background must not be treated as benign support;
- limited background must not justify hard PASS/artifact calls;
- legacy Branch B kept/review status must not define N0 or benchmark Branch B
  V2 performance;
- reference-included samples can support residual/batch-stability review but
  are not ideal empirical-null controls unless the design is cross-fit or
  leave-one-out.

## Branch B Constraints

Branch B is a candidate-level evidence and disposition layer for Branch A
events.

Allowed Branch B evidence:

- amplitude/CN-like evidence;
- mosaic-like evidence;
- LOH/UPD annotation where supported;
- region-risk and mask/refmap context;
- GC/RC, sample noise, waviness, and batch-risk context;
- CNVpro-inspired consistency evidence;
- matched-negative or percentile evidence only when the background source is
  valid and labeled.

Branch B must not:

- create standalone final candidates;
- act as a replacement for WisecondorX predict/CBS;
- use hard thresholds derived from the current exploratory five-sample set;
- suppress weak positives before no-FN ablation on Y1-Y8 and H1-H6;
- claim final FP reduction from legacy Branch B or N1-only background;
- report PASS/CONFIRMED labels when evidence is limited or unknown.

Branch B report-contract iteration can pass only if:

- every Branch A candidate has a Branch B evidence/disposition row;
- reportable events are derived from Branch A candidates plus Branch B evidence,
  not from the legacy kept-event path alone;
- the report states reference ID, binsize, mask/preprocess contract,
  background source, and evidence version;
- missing refmap, calibration-null, matched-background, or SCA evidence remains
  visible.

## Branch S / SCA Constraints

Branch S is the sex-chromosome/SCA evidence workflow. It must be separated from
autosomal/segmental Branch B promotion.

Current Branch S output remains review/development-only until locked SCA
validation passes.

Branch S must:

- use sex-aware expected ploidy;
- keep X non-PAR and Y non-PAR evidence separate;
- treat PAR / XY-homology evidence as annotation, consistency, and artifact-risk
  context rather than a primary whole-SCA seed;
- state whether each SCA result is final-reportable, review-only, or
  development-only;
- preserve current sex-calling behavior unless a validated replacement is
  explicitly promoted.

Branch S cannot be promoted without locked truth coverage for X gain, XXY, XYY,
Y loss, mosaic SCA, and clean XX/XY negatives across batches.

## Current Cannot-Promote / Cannot-Claim List

The current R&D state cannot claim:

- the H-augmented reference is final;
- Branch A/reference contract is frozen;
- Branch B V2 improves final-level reporting over legacy output;
- N1/reference-cohort background is formal N0;
- matched-negative percentile evidence is ready for hard filtering;
- Branch S/SCA is final-reportable;
- the 2026-06-15 five-sample set is a clean locked validation set;
- old Branch B kept-event counts are a valid benchmark truth table;
- zero Branch B reportable events prove a sample is a clean negative.

The current R&D state can continue only as:

- reference cohort audit;
- Branch A no-FN validation under a named reference;
- Branch B candidate evidence refinement without B-only final calls;
- CNVpro-inspired consistency comparison;
- Branch S/SCA report-boundary development.

## Do-Not-Do List

- Do not use the removed legacy Phase 1-5 route as the current plan.
- Do not promote reports from current Branch B V2/N1 shadow output.
- Do not use old-reference `N0=0` as a post-rebuild blocker.
- Do not treat H7-H16 as automatic N0 or automatic production reference
  samples.
- Do not introduce new hard Branch B thresholds before known-positive no-FN
  ablation.
- Do not report CNVpro-inspired output as a replacement CNV call.
- Do not leave SCA undefined in report schema; if not validated, label it as
  review/development-only with explicit limitations.
