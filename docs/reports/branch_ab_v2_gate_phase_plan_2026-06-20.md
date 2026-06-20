# Branch A/B/S Gate And Phase Plan

Date: 2026-06-20

Status: `active_gate_contract`

Decision use: current gate/phase order. Pair with
`docs/CURRENT_CONTEXT_INDEX.md` for the latest Branch S/P6 report-direction
override.

## Purpose

This file is the current gate/phase control document for the PGT-A CNV-seq R&D
cycle after the 2026-06-20 plan correction.

It supersedes the old interpretation that the next step is simply "promote
Branch B V2 report contract". The active direction is:

```text
Branch A / reference reset
-> Branch B candidate evidence refinement
-> CNVpro-inspired consistency comparison
-> Branch S / SCA workflow finalization
-> report package after no-FN validation
```

Older phase implementation reports remain useful execution history, but this
file defines the current gate order and promotion rules.

## Global Invariants

These apply to every phase:

1. WisecondorX predict remains the primary CNV caller.
2. Branch A is the sensitivity gate and is derived from WisecondorX predict/CBS
   evidence.
3. Branch B refines Branch A candidates only; it must not create B-only final
   report events.
4. Branch B may add LOH / UPD / mosaic / CN-like amplitude, region-risk,
   sample-noise, background, and CNVpro-inspired consistency evidence.
5. Hard filtering must be narrow and must not increase known-positive FN.
6. Current N0/cross-fit background is not complete. It remains the clean
   promotion path, but it is not an immediate blocker for Branch B evidence
   development.
7. Reference-cohort and N1 evidence are limited context only and must be labeled
   as such.
8. Branch S/SCA must become a reportable workflow eventually, but current
   Branch S output is review/development-only until validation gates pass.
9. Every real validation, dry-run, and workflow test must run on `ssh fengxian`.
10. Reference changes have two rerun levels:
    - P1/P2 audit rerun: WisecondorX predict and Branch A are mandatory; current
      Branch B/S/report outputs may be materialized only as diagnostic snapshots
      and must be labeled `legacy/current-code context`.
    - P6 report-facing rerun: WisecondorX predict, Branch A, fixed Branch B,
      fixed Branch S, evaluation, benchmark, and report must all run under the
      same reference/config contract.

## Definitions

| term | meaning |
|---|---|
| Branch A | WisecondorX predict/CBS-derived candidate discovery layer. |
| Branch B | Candidate refinement/evidence/disposition layer for Branch A events. |
| Branch S | Sex chromosome/SCA evidence and reporting workflow. |
| R0/R1/R2 | Reference-rebuild eligibility labels. |
| N0/N1/N2 | Negative/background labels for Branch B calibration or context. |
| CNVpro-inspired | Method reference or consistency comparison; not a replacement caller. |
| final-reportable | Allowed to appear as final report interpretation. |
| review-only | Visible to report/reviewer but not final confirmation or suppression. |
| development-only | Internal R&D evidence, not report interpretation. |

## Phase Overview

| phase | module | goal | entry gate | exit gate |
|---|---|---|---|---|
| P0 | context | Lock the current constraints and inputs. | Latest handoff and R&D docs read. | G0 context ready. |
| P1 | reference / Branch A | Audit candidate reference samples and decide the reference cohort. | G0 | G1 reference cohort decision. |
| P2 | Branch A | Freeze Branch A reference/config and validate no-FN sensitivity. | G1 | G2 Branch A fixed. |
| P3 | Branch B | Add candidate-level evidence refinement without hard promotion. | G2 | G3 Branch B evidence complete. |
| P4 | CNVpro comparison | Run CNVpro-inspired consistency checks against fixed Branch A. | G3 | G4 comparison interpretable. |
| P5 | Branch S / SCA | Define and validate reportable SCA workflow boundaries. | G2 or G3 | G5 SCA report boundary. |
| P6 | report | Generate 0615 five-sample report package under fixed contracts. | G2 + G3 + G5 | G6 report release decision. |

## P0: Context Lock

### Required Inputs

- latest valid handoff;
- `AGENTS.md`;
- `skills/conversation_handoff/SKILL.md`;
- `skills/pgta_reference_modeling_analysis/SKILL.md`;
- current R&D documents;
- current config files for the reference and predict runs;
- current remote output paths.

### G0 Pass Criteria

P0 passes only if the worker can state:

- which handoff is active;
- which reference IDs and result directories are current;
- which files are only historical implementation reports;
- which outputs are current result files versus planned work;
- which work must run remotely.

### If G0 Fails

Do not edit code or run analysis. Read the missing context first and document
the gap.

## P1: Reference Cohort Audit

### Goal

Evaluate the current H-augmented shadow reference and all candidate reference
samples before freezing Branch A.

### Inputs

- original retained normal/reference samples;
- H7-H16 candidate negative samples;
- current H-augmented shadow reference samples;
- QC summaries;
- sex calls and BAM X/Y depth consistency;
- WisecondorX predict output under the evaluated reference;
- Branch A candidate tables;
- region-risk, mask, refmap, acrocentric/high-repeat evidence.

### Required Output

A reference audit table with one row per candidate sample:

```text
sample
batch
qc_status
sex_call
reference_membership
branch_a_candidate_count
branch_a_strong_count
branch_a_sensitive_count
broad_event_count
altered_genome_fraction
acrocentric_or_high_repeat_burden
centromere_telomere_proximal_burden
shared_batch_signal_flag
sample_specific_signal_flag
R_label
R_label_reason
recommended_action
```

### G1 Pass Criteria

G1 passes only if:

- R0/R1/R2 labels are assigned from post-rebuild/ref-specific evidence, not old
  Branch B kept counts alone;
- potentially abnormal samples are identified for exclusion or R1-only ablation;
- H7-H16 are not rejected solely because the old 32-sample reference produced
  Branch A calls;
- the reference cohort decision is written to a stable cohort/config file;
- the chosen reference variant has a unique reference ID.

### If G1 Fails

Do not promote any H-expanded reference. Keep the current H reference as shadow
only and collect the missing audit evidence.

## P2: Branch A Freeze And No-FN Validation

### Goal

Fix the Branch A reference/config contract and prove that known positives remain
detected.

### Inputs

- G1-approved reference cohort;
- WisecondorX reference build output;
- WisecondorX predict config;
- Y1-Y8 truth table;
- H1-H6 truth table;
- H7-H16 and 0615 exploratory samples for burden review.

### Required Output

A Branch A validation summary:

```text
reference_id
binsize
preprocess_mask_version
wisecondorx_predict_command
blacklist_status
minrefbins
zscore
alpha
maskrepeats
sample
truth_event_count
truth_detected_count
branch_a_candidate_count
branch_a_strong_count
branch_a_sensitive_count
top_branch_a_signal
H6_chr21_status
FN_count
```

### G2 Pass Criteria

G2 passes only if:

- Y1-Y8 and H1-H6 known-positive Branch A recall does not regress;
- H6 chr21 weak/focal gain is explicitly tracked;
- WisecondorX predict and Branch A use the same reference ID and config;
- any optional current-code Branch B/S/report diagnostic snapshots are labeled
  diagnostic and are not used as G2 pass/fail criteria;
- Branch A candidate generation remains WisecondorX-first;
- Branch A does not include Branch B artifact filters;
- per-sample candidate burden is reported for H7-H16 and 0615 samples.

### If G2 Fails

Do not continue to FP suppression. Revisit reference cohort, binsize, mask, or
WisecondorX predict settings.

## P3: Branch B Evidence Refinement

### Goal

Turn Branch B into a candidate-level evidence/disposition layer for Branch A
events, not a hard-threshold patch collection.

### Inputs

- fixed Branch A candidate table from G2;
- existing Branch B evidence ledger;
- correction/calibration outputs where available;
- region-risk, refmap, mask, high-repeat, acrocentric, centromere/telomere
  proximity evidence;
- reference-cohort or N1 context if available;
- LOH / UPD / mosaic evidence sources when implemented.

### Required Output

One Branch B evidence/disposition row per Branch A candidate:

```text
candidate_id
sample
chrom
start
end
state
a_zscore
a_support_level
branch_b_direction_support
copy_number_like_amplitude
mosaic_proxy
loh_evidence
upd_evidence
background_source
background_status
region_risk_context
sample_noise_context
cnvpro_consistency_status
disposition
disposition_reason
report_impact
```

### G3 Pass Criteria

G3 passes only if:

- every Branch A candidate has one Branch B evidence row;
- no B-only final event is created;
- missing evidence remains `UNKNOWN` or review/no-call;
- LOH / UPD / mosaic fields are present as available or explicitly `not_available`;
- known-positive truth-overlap candidates are not hard-suppressed;
- hard thresholds are not introduced without ablation showing zero FN impact;
- review burden is summarized per sample.

### If G3 Fails

Keep Branch B output review/development-only. Do not use it to suppress or
confirm report events.

## P4: CNVpro-Inspired Consistency Comparison

### Goal

Use CNVpro methods as an independent consistency/reference track after Branch A
is fixed.

### Inputs

- fixed Branch A candidates;
- CNVpro or CNVpro-like outputs where available;
- GC/RC/waviness evidence;
- p/q-arm and acrocentric partition evidence;
- sex-aware X/Y and PAR/non-PAR evidence;
- large-event tiering and copy-number/mosaic proxy.

### Required Output

For each Branch A candidate:

```text
candidate_id
cnvpro_direction
cnvpro_size_tier
cnvpro_z_or_score
cnvpro_copy_number_proxy
gc_rc_waviness_context
arm_partition_context
acrocentric_context
same_direction_support
discordance_reason
usable_for_report_context
```

### G4 Pass Criteria

G4 passes only if:

- CNVpro-inspired checks are clearly labeled as comparison/context;
- discordance does not automatically suppress Branch A positives;
- agreement can raise review confidence but does not replace WisecondorX;
- small-region and high-repeat signals are interpreted with PGT-A resolution
  policy, usually prioritizing stable 2-4 Mb or larger signals unless evidence
  is strong.

### If G4 Fails

Do not use CNVpro-inspired fields in report disposition. Keep them as benchmark
appendices only.

## P5: Branch S / SCA Workflow Finalization

### Goal

Define how SCA evidence enters report output.

### Inputs

- sex call and sex-routing output;
- sex-specific WisecondorX predict context;
- Branch A sex-chromosome candidates;
- Branch B calibrated sex-chromosome evidence;
- X non-PAR, X PAR, Y non-PAR, and Y/PAR or homology-region annotations;
- known SCA truth cases and clean XX/XY negatives where available.

### Required Output

One SCA evidence summary per sample:

```text
sample
sex_call
expected_x_ploidy
expected_y_ploidy
x_nonpar_direction
x_par_context
y_nonpar_direction
y_par_or_homology_context
branch_a_x_support
branch_a_y_support
sca_candidate_state
sca_confidence_tier
sca_output_mode
sca_uncertainty_reason
report_text_status
```

### G5 Pass Criteria

G5 has two levels:

**G5-review passes** if:

- every SCA statement has region source and uncertainty reason;
- X non-PAR and Y non-PAR are separated from PAR/homology context;
- current unvalidated outputs are marked review/development-only;
- SCA rows do not silently add or remove final CNV events.

**G5-final passes** only if:

- locked SCA positives cover X loss, X gain, XXY/XXX-like, XYY-like, Y loss,
  and mosaic SCA cases as available;
- clean XX and XY negatives across batches do not become false SCA calls;
- sex-call concordance does not regress;
- final report schema is validated remotely.

### If G5 Final Fails

Report SCA as review/development evidence only, while still including a clear
SCA section if the report format requires it.

## P6: Report Package

### Goal

Generate report output for 0615 or later samples only after the upstream
contracts are fixed.

### Inputs

- G2 fixed Branch A reference/config;
- G3 Branch B evidence/disposition table;
- G4 CNVpro-inspired comparison where available;
- G5 Branch S/SCA output mode;
- benchmark and evaluation summaries.

### Required Output

Report package must include:

- reference ID and WisecondorX predict command;
- Branch A candidate summary and top signals per sample;
- Branch B disposition/evidence counts;
- LOH / UPD / mosaic annotation status;
- CNVpro-inspired agreement/discordance context;
- SCA section with final/review/development status;
- known limitations and background source labels;
- evidence that Y1-Y8 and H1-H6 no-FN gates still pass under the same
  reference/config.

### G6 Pass Criteria

G6 passes only if:

- report is generated by workflow outputs, not hand-edited tables;
- all components use the same reference ID and config;
- every Branch A candidate has a Branch B disposition/evidence row;
- SCA status is explicitly final/reportable, review-only, or development-only;
- missing evidence is visible and not treated as benign;
- report does not claim locked-validation performance from exploratory 0615
  samples.

### If G6 Fails

Do not hand out the report as the forward workflow result. Use it only as an
internal exploratory comparison.

## Immediate Execution Order

1. Implement P1 reference-sample Branch A audit.
2. Run remote unit tests and Snakemake dry-run on `fengxian`.
3. Materialize the audit for current H-augmented shadow reference and candidate
   expanded reference sets.
4. Decide R0/R1/R2 and freeze a named Branch A/reference contract.
5. Rerun WisecondorX predict and Branch A for Y1-Y8, H1-H6, H7-H16, and 0615
   under the fixed reference.
6. Only then implement P3 Branch B evidence refinements and P5 SCA report schema.

## Current Do-Not-Do List

- Do not treat the current H-augmented reference as final.
- Do not use old Branch B kept counts as the N0 benchmark.
- Do not make complete N0/cross-fit background a blocker for all development.
- Do not introduce new hard Branch B thresholds before no-FN ablation.
- Do not report CNVpro-inspired results as replacement calls.
- Do not leave SCA permanently undefined; define report schema even if current
  evidence remains review/development-only.
- Do not release 0615 as an old current-workflow report.
