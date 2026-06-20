# Branch B V2 Reference-Cohort Background And SCA Design

Date: 2026-06-20

## Scope

This document updates the Branch B V2 design after the H R0 shadow-reference
discussion.

The goal is to make Branch B V2 useful for candidate triage without violating
the current project constraints:

- WisecondorX predict remains the primary CNV caller.
- Branch A remains the sensitivity gate.
- Branch B V2 refines Branch A candidates and must not create B-only final
  candidates.
- Branch S remains sex-chromosome/SCA evidence and must not replace sex calling
  or final report logic before locked SCA validation.
- Any reference change first requires rerunning WisecondorX predict and Branch A
  candidates for audit. Fixed Branch B, Branch S, evaluation, benchmark, and
  report must be rerun before any report-facing release or promotion decision.

This is a design/feasibility note only. It does not claim that the proposed
rules have already been validated.

## Context Read

Primary files reviewed before this design:

- `docs/handoff/2026-06-20_0105_branch_b_v2_n1_shadow_background_handoff.md`
- `AGENTS.md`
- `skills/conversation_handoff/SKILL.md`
- `skills/pgta_reference_modeling_analysis/SKILL.md`
- `docs/reports/branch_ab_v2_rnd_constraints_2026-06-18.md`
- `docs/reports/h7_h16_reference_rebuild_eligibility_2026-06-19.tsv`
- `docs/reports/h_r0_shadow_post_rebuild_predict_2026-06-19.md`
- `docs/reports/branch_b_v2_n1_shadow_background_validation_2026-06-20.md`
- `docs/reports/branch_s_sca_validation_plan_2026-06-19.md`
- `docs/reports/branch_s_materialization_2026-06-19.md`
- `pgta/predict/branch_b/negative_bank.py`
- `pgta/predict/branch_b/matched_negative.py`
- `pgta/predict/branch_b/v2_classifier.py`
- `pgta/predict/branch_b/artifact_rules.py`
- `pgta/predict/branch_s.py`
- `pgta/predict/sex_routing.py`
- `rules/reference_layout.smk`
- `rules/predict_layout.smk`
- `rules/predict_workflow.smk`

## Current Problem

Current Branch B V2 is too conservative to be useful as a decision layer:

- Formal `N0` matched-negative background is empty.
- The latest implementation allows only configured `N1` rows as
  `SHADOW_BACKGROUND`.
- `SHADOW_BACKGROUND` remains `none_shadow_only`.
- Most truth-overlap positive candidates stay
  `UNKNOWN_BACKGROUND_POSITIVE_SUPPORT`, so V2 has not demonstrated a final-level
  improvement over the legacy Branch B report path.

The design issue is not only a threshold issue. The background source is too
narrow and is tied to legacy negative-bank labels that were derived under an
older reference/Branch B context.

The reference issue is related:

- H7-H16 are presumed negative by sample identity, but the old reference
  produces multiple Branch A signals in them.
- Old-reference Branch A signal alone should not reject a negative sample from
  reference-rebuild exploration.
- However, if a sample has a real broad sample-specific abnormality, modeling it
  into the normal reference would teach the reference the wrong distribution.

Therefore the next design must separate:

1. reference-rebuild eligibility (`R0/R1/R2`);
2. formal matched-negative calibration (`N0/N1/N2`);
3. reference-cohort background evidence for V2 review;
4. sex-aware Branch S/SCA review.

## Approach Options

### Option A: Keep Strict N0-Only Background

This is safest statistically, because only locked clean negatives become
empirical null samples.

Problem: current `N0=0`, so Branch B V2 remains mostly
`UNKNOWN_BACKGROUND` and cannot materially reduce FP or improve report triage.

Use only as the final promotion gate, not as the only development path.

### Option B: Use Reference-Cohort Limited Background

Use all selected normal reference inliers from the same reference ID as a
limited background source, while recording whether each background sample was
included in the reference and excluding the query sample when applicable.

This is less pure than held-out N0, but it answers the practical development
question: whether a candidate is unusual relative to the current normal model.

This is the recommended next step.

It must be named explicitly as reference-cohort background, not as formal N0:

```text
REF_COHORT_BACKGROUND
REF_COHORT_LIMITED_BACKGROUND
REF_COHORT_OUTLIER_POSITIVE_SUPPORT
REF_COHORT_BACKGROUND_COMPATIBLE
```

These tiers can support review/triage decisions and threshold development, but
not final PASS/artifact promotion until ablation is complete.

### Option C: Cross-Fit Or Leave-One-Out Background

Build folds where each candidate background sample is predicted against a
reference that did not include itself.

This is the cleanest way to use reference samples as negative background. It is
also more expensive because it requires multiple reference builds and full
predict reruns.

Use this as the promotion path after Option B identifies promising rules.

## Recommended Design

Use the reference-cohort background idea only after the Branch A/reference
contract has been audited and fixed. It is a limited context source, not the
main development objective and not a substitute for a formal N0/cross-fit
design.

The immediate design order is:

1. audit the current H-augmented shadow reference and candidate expanded
   reference sets;
2. freeze a named Branch A/reference contract;
3. rerun WisecondorX predict and Branch A candidates under that contract;
4. then evaluate Branch B candidate refinement, including LOH / UPD / mosaic and
   CNVpro-inspired evidence fields;
5. keep Option C as the final promotion path if reference-cohort context proves
   useful.

The revised V2 background stack should be:

| source | meaning | allowed decision impact now |
|---|---|---|
| `N0` / cross-fit held-out background | Formal empirical null. Query was not used to train its own background. | Eligible for future final promotion after validation. |
| `REF_COHORT_BACKGROUND` | Normal reference inliers from the same reference ID, excluding query sample if present. | Review/triage evidence only. |
| `REF_COHORT_LIMITED_BACKGROUND` | Same as above, but effective background count is small. | Positive/review support only; no hard artifact. |
| `N1_SHADOW_BACKGROUND` | Presumed negative shadow context from manually configured N1 labels. | Review context only. |
| `UNKNOWN_BACKGROUND` | Not enough comparable rows. | Review/no-call; never clean support. |

This makes the design wider without making it unsafe.

## Reference Sample Abnormality Audit

Before promoting a new reference, every candidate reference sample should be
audited under the reference version being evaluated.

### Required Audit Inputs

- QC PASS / WARN / FAIL.
- Sex call and BAM X/Y depth consistency.
- Reference cohort membership and sex group.
- WisecondorX predict output using the candidate reference.
- Branch A candidate table from the same predict run.
- Branch B evidence ledger and artifact-rule output from the same predict run.
- Region-risk/refmap/mask evidence for recurrent unstable regions.

Existing BAMs can be reused if the alignment FASTA and mapping contract have not
changed. BAM regeneration is required only if the alignment reference, trimming,
mapping command, or BAM contract changes.

### Branch A Audit Metrics For Reference Samples

For each reference candidate sample, record:

- total Branch A candidate count;
- strong Branch A count (`abs(a_zscore) >= 10`);
- sensitive Branch A count (`abs(a_zscore) >= 5`);
- total altered genome fraction from Branch A candidates;
- broad event count, especially events spanning a large chromosome fraction;
- recurrent-region burden in chr13/14/15 short arms, high-repeat / low
  mappability regions, centromere/telomere-proximal bins, and unstable
  chromosomes such as chr6/9/11/16;
- whether the signal is shared across many candidate reference samples from the
  same batch;
- whether the signal is sample-specific and inconsistent with batch/reference
  artifact behavior.

### R0/R1/R2 Interpretation

Use Branch A audit as a reference-stability screen, not as a direct truth call.

| label | interpretation | action |
|---|---|---|
| `R0` | QC PASS, sex concordant, no broad sample-specific abnormality after post-rebuild audit. | Eligible for primary shadow reference. |
| `R1` | QC PASS, likely negative, but has broad/recurrent/batch-risk signal that needs ablation. | Eligible for expanded shadow reference only. |
| `R2` | QC fail, identity/sex conflict, or strong sample-specific abnormality. | Exclude from reference rebuild. |

Old-reference Branch A signal alone must not force `R2`. The key question is
whether the signal persists after a named shadow reference rebuild and whether
it is sample-specific.

If a sample is used inside the same reference being tested, its post-rebuild
predict result is a stability check, not a clean independent abnormality test.
For stronger evidence, use leave-one-out or cross-fit.

### Rebuild Loop

The reference rebuild loop should be:

1. Start with a documented candidate table of `R0/R1/R2`.
2. Build a named shadow reference, for example:
   - conservative `R0-preferred`;
   - expanded `R0+R1`.
3. Rerun WisecondorX predict for all validation samples.
4. Regenerate Branch A candidates.
5. Optionally regenerate current-code Branch B/S/report outputs as diagnostic
   snapshots only.
6. Evaluate Y1-Y8 and H1-H6 Branch A truth recall.
7. Audit H7-H16 and all reference-included samples for residual Branch A
   signals.
8. If a candidate reference sample shows a persistent sample-specific broad
   positive signal, move it to `R2`, rebuild, and rerun the P1/P2 audit chain;
   current-code Branch B/S outputs remain optional diagnostic snapshots until
   their fixed contracts pass later gates.

This prevents abnormal samples from being modeled into the reference while also
avoiding the previous bias of rejecting new-batch negatives solely because the
old 32-sample reference generated Branch A calls.

## Branch B V2 Classification Update

Branch B V2 should classify only Branch A candidates. It should not add B-only
events to final reports.

This branch should not be reduced to a percentile/background hard-filter. It
must also carry candidate-level biological and technical annotations that are
needed for final review:

- LOH / UPD evidence when available;
- mosaic-like amplitude and dilution-sensitive support;
- copy-number-like amplitude relative to sex/autosome context;
- GC/RC/waviness and sample-noise context;
- CNVpro-inspired p/q-arm, acrocentric, PAR/non-PAR, and large-event tier
  consistency;
- explicit unknown/missing evidence flags.

### Candidate Evidence Gates

Proposed evidence gates:

| gate | condition | V2 output |
|---|---|---|
| `A_STRONG_B_CONSISTENT` | `abs(a_zscore) >= 10` and Branch B signal is same direction. | High-priority positive review; future confirmation candidate after validation. |
| `A_SENSITIVE_REF_OUTLIER` | `abs(a_zscore) >= 5`, same direction support exists, and reference-cohort background is limited/outlier. | Positive support review. |
| `A_SENSITIVE_UNKNOWN_BACKGROUND` | `abs(a_zscore) >= 5`, but background is unknown. | Review/no-call, not artifact. |
| `BACKGROUND_COMPATIBLE_WEAK_A` | weak A support and reference/matched background compatible. | Artifact/reject candidate after validation. |
| `TECHNICAL_CONTRACT_RISK` | same-chrom reference leakage, corrupt input, QC fail, or unusable region. | No-call/artifact depending on severity. |
| `DIRECTION_CONFLICT` | Branch A direction conflicts with reliable Branch B/ref-cohort direction. | No-call/review; hard artifact only after validation. |

### Wide But Safe Threshold Policy

The next V2 design should be less strict than the current matched-negative
`N0-only` logic:

- `abs(a_zscore) >= 5` remains the sensitivity floor inherited from Branch A.
- `abs(a_zscore) >= 7` is an important weak-positive review tier, because H6
  chr21 gain has `a_zscore=7.11` and is truth-positive.
- `abs(a_zscore) >= 10` is strong Branch A support, not a standalone PASS rule.
- If background count is below 5, do not call the event clean. Keep it
  `UNKNOWN_BACKGROUND` or `REF_COHORT_LIMITED_BACKGROUND`.
- If background count is 5-19, use the background only as limited review
  evidence. Being above all comparable background rows can support
  `REF_COHORT_OUTLIER_POSITIVE_SUPPORT`, but cannot by itself produce final
  `CONFIRMED`.
- If background count is at least 20 and the query is above p95/p99, the
  evidence is stronger, but still requires no-FN validation before final
  promotion.
- Background-compatible weak A candidates can be downgraded, but only after
  truth-positive recall and held-out negative burden are reported.

This policy is deliberately not result-backfit. The thresholds are derived from
Branch A sensitivity requirements, the known H6 weak-positive case, and finite
background sample-size limits.

### H6 Chr21 Example

Under this design, H6 chr21 gain should not be hard-filtered merely because its
Branch B z support is weak.

Expected V2 behavior:

```text
Branch A z = 7.11
background unknown or limited -> REVIEW_REQUIRED
background limited and query above comparable ref-cohort values -> REF_COHORT_OUTLIER_POSITIVE_SUPPORT
background compatible but no reliable contradiction -> REVIEW_REQUIRED, not artifact
```

This preserves the no-FN goal while still letting Branch B reduce weak,
background-compatible FP later.

## Branch S / SCA Update

Current workflow already supports sex-specific WisecondorX predict when
`reference_output_by_sex.XX` and `reference_output_by_sex.XY` are configured:

```text
WisecondorX gender
-> choose XX or XY reference through predict_gender
-> WisecondorX predict --gender <F/M>
-> Branch A candidates
-> Branch B / Branch S evidence
```

Current Branch S does not build a separate SCA reference. It consumes:

- Branch B calibrated sex-chromosome bins;
- Branch A sex-chromosome candidates;
- gender/sex-call TSV.

This is acceptable for shadow evidence, but not enough for final SCA calling.

### SCA Copy-Number Interpretation Model

The first Branch S draft treated PAR / XY shared or homologous dosage as a
possible layer for whole sex-chromosome state review. That statement was too
strong. The current evidence supports a more conservative model:

1. **Sex-status / karyotype-level SCA review**
   - Use sex-aware whole-X and whole-Y dosage evidence, Y fraction / X-Y depth
     ratios, and WisecondorX sex-specific predict context to review states such
     as `XO`, `XXY`, `XYY`, and broad sex-chromosome dosage imbalance.
   - PAR / XY shared or homologous regions must be treated as annotation,
     consistency, and artifact-risk context, not as the primary whole-SCA call
     seed. These regions are biologically pseudo-autosomal and technically
     sensitive to reference masking and ambiguous alignment.

2. **Chromosome-specific copy-number direction**
   - Use X non-PAR bins as the main X copy-number direction evidence.
   - Use Y non-PAR bins as the main Y copy-number direction evidence.
   - Interpret X and Y non-PAR evidence against sex-specific expected ploidy:
     `XX`: X ~= 2, Y ~= 0; `XY`: X ~= 1, Y ~= 1.
   - Treat X PAR as pseudo-autosomal / diploid-like context, especially in XY
     samples, rather than merging it with X non-PAR dosage evidence.

Autosomal bins remain useful for global depth normalization and sample-level
noise scale, but autosome-normalized X/Y z alone is not enough for final SCA CN.
It must be combined with sex-specific reference context, sex call evidence, and
PAR/non-PAR annotations.

Evidence review:

- CNVpro `CNVcalling.R` infers sex using Y read count over autosomal read count
  and computes sex-aware X/Y references before CNV calling.
- CNVpro segments chrX as `par1`, `pter`, `qter`, `par2`, writes
  `*.withpar.tsv`, then removes chrX PAR rows from final `*_bin.tsv` and
  `*_seg.tsv`.
- CNVpro `filterCNV.py` uses sex-specific CN baselines: male X non-PAR and Y
  are evaluated around one copy, while autosomes and PAR-like X evidence use a
  diploid baseline.
- Local cnvpro output inspection confirmed that `Y2` and `Y3` lose 262 chrX PAR
  bins from `*_bin.withpar.tsv.gz` to final `*_bin.tsv.gz`, while X non-PAR and
  Y bins remain.
- Public tool documentation is consistent with this direction: WisecondorX
  requires gender handling for X/Y prediction; DRAGEN CNV handles sex
  chromosomes and PAR based on sample sex; CNVkit documents PAR reference
  masking as a source of sex-chromosome copy-number bias.

External references reviewed:

- WisecondorX README:
  <https://github.com/CenterForMedicalGeneticsGhent/WisecondorX>
- DRAGEN CNV calling documentation:
  <https://help.dragen.illumina.com/dragen-v4.3/product-guide/dragen-v4.3/dragen-dna-pipeline/cnv-calling>
- CNVkit chromosomal sex and PAR handling:
  <https://cnvkit.readthedocs.io/en/master/sex.html>
- GATK reference genome components and PAR masking:
  <https://gatk.broadinstitute.org/hc/en-us/articles/360041155232-Reference-Genome-Components>
- SCA/NIPT caution context:
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC10491423/>

Therefore Branch S V2 should not implement a simplistic "PAR layer decides
overall SCA" rule. PAR/XY-homology should be a region-risk and consistency
feature, while final-oriented SCA review must be based on sex-aware X/Y dosage
with non-PAR evidence kept separate.

### PGT-A Resolution Policy

CNVpro targets high-resolution CNV-seq calling and can operate around
approximately 100 kb evidence units. The PGT-A workflow does not need to keep
that same resolution target for final reporting.

For this project, the SCA/Branch B V2 design should prioritize stable broad
signals:

- final-oriented review should focus on events at least `2-4 Mb`, unless a
  smaller event has very strong multi-evidence support;
- 100 kb / sub-Mb signals can remain in internal evidence tracks, plots, and
  review diagnostics;
- high-repeat, PAR, XY-homology, acrocentric, centromere/telomere-proximal, and
  low-refbin regions should increase review/artifact caution rather than
  silently remove strong broad signals;
- bin-level sensitivity should not be optimized at the cost of recurrent FP in
  unstable regions.

This keeps the algorithmic direction close to CNVpro while matching the PGT-A
use case: tolerate lower breakpoint precision to reduce unstable small-region
FP, while preserving broad and mosaic-positive sensitivity.

### Required SCA V2 Direction

Branch S should be upgraded as a sex-aware review model:

1. Keep sex-specific WisecondorX predict as the primary input path.
2. Build sex-specific reference-cohort background summaries:
   - XX background for X dosage in XX samples;
   - XY background for X/Y dosage in XY samples;
   - separate X non-PAR, X PAR, and Y non-PAR regions;
   - treat Y PAR / XY-homology according to the active reference/mask asset
     rather than assuming it is directly countable.
3. For suspicious or discordant SCA cases, optionally run a dual-reference
   shadow comparison:
   - query with XX reference;
   - query with XY reference;
   - report the difference as review evidence only.
4. Keep Branch S labels shadow-only until locked SCA truth coverage exists.

The dual-reference comparison should not replace the primary sex call. It is an
explainability tool for SCA review and for detecting cases where current sex
routing hides or exaggerates X/Y dosage.

### SCA Promotion Gate

Branch S cannot become final until validation includes:

- clean XX negatives across batches;
- clean XY negatives across batches;
- X loss / XO-like positives;
- X gain / XXX-like positives;
- XXY-like positives;
- XYY-like positives;
- Y loss or structurally abnormal Y cases;
- mosaic SCA series;
- PAR / XY-homology edge cases.

Until then, Branch S can improve review prioritization, not final PASS/FAIL.

## Feasibility

The design is feasible in the current architecture.

Existing code already has most required boundaries:

- Branch A comes from WisecondorX predict and `pgta/predict/branch_a.py`.
- Branch B evidence is candidate-level.
- `matched_negative.py` already supports background ledgers and excludes the
  query sample.
- `v2_classifier.py` already separates evidence tier from final impact.
- `negative_bank.py` already records N0 readiness.
- `branch_s.py` already outputs shadow-only SCA evidence.
- `predict_workflow.smk` already supports sex-specific reference selection.

The main implementation gap is not architecture. It is report-contract and
background semantics:

- introduce reference-cohort background as a named source;
- keep it distinct from formal N0/cross-fit background;
- update V2 tiers so limited background can support review but not hard
  artifact/PASS;
- connect `cnv_report` to Branch B V2 classifier output so future report
  packages are V2-derived rather than legacy Branch B kept-event reports;
- add reference-sample Branch A audit outputs tied to a reference ID.

## Proposed Implementation Order

1. Add a reference-sample Branch A audit report.
   - Per-sample A candidate count, strong count, broad burden, recurrent-region
     burden, and sample-specific versus shared-batch signal flags.
   - Output `R0/R1/R2` recommendations for the next shadow rebuild.
   - This must run before freezing a new Branch A/reference contract.

2. Freeze and document the Branch A/reference contract.
   - Input cohort and exclusions.
   - Reference ID, binsize, preprocessing/mask version, WisecondorX predict
     command, sex-specific reference behavior, and blacklist status.
   - Full Y1-Y8 / H1-H6 recall and H7-H16/reference-sample burden summaries.

3. Add a reference-cohort background manifest as a limited context source.
   - Input: selected/inlier reference sample IDs, sex group, reference ID,
     whether the sample was used in the reference, and evidence ledger path.
   - Output: a stable TSV consumed by Branch B V2.
   - Do not label this as formal N0.

4. Extend matched-negative/background logic.
   - Add `background_source` values:
     `N0`, `REF_COHORT`, `REF_COHORT_LIMITED`, `N1_SHADOW`, `UNKNOWN`.
   - Preserve `N0` as the only formal empirical null.
   - Exclude query sample from background when present.
   - Treat missing or limited background as review/no-call context, not as clean
     support.

5. Extend Branch B evidence and classifier tiers.
   - Add `REF_COHORT_OUTLIER_POSITIVE_SUPPORT`.
   - Add `REF_COHORT_BACKGROUND_COMPATIBLE`.
   - Add `REF_COHORT_LIMITED_REVIEW`.
   - Add LOH / UPD / mosaic and CNVpro-inspired consistency fields where the
     upstream evidence exists.
   - Keep hard PASS/artifact promotion blocked until validation passes, but
     allow V2 review classes to drive the V2 report contract.

6. Connect the agreed Branch B evidence/disposition table to report generation.
   - Add the evidence/disposition file as a `cnv_report` input for the new
     report config.
   - Report per-sample V2 class counts, evidence-tier counts, review priorities,
     and background source (`N0`, cross-fit, reference-cohort, N1 shadow, or
     unknown).
   - Preserve Branch A candidate visibility so zero V2-reportable events does
     not silently imply the sample is a clean validation negative.
   - Do not allow unvalidated Branch S output to add or remove final report
     events.

7. Finalize Branch S/SCA workflow design.
   - Keep current sex-specific predict path.
   - Add XX/XY reference-cohort SCA background tables.
   - Define which SCA outputs can enter report as final calls versus review
     labels.
   - Optional dual-reference shadow comparison only after the main V2
     background design is stable.

8. Run remote validation on `fengxian`.
   - Unit tests first.
   - Snakemake dry-runs for Y/H/0615 configs.
   - Full WisecondorX predict rerun after any reference change.
   - Evaluate Y1-Y8 and H1-H6 recall before looking at FP reduction.

## Validation Requirements

Any implementation of this design must report:

- reference ID and input cohort;
- whether predict used XX/XY sex-specific references;
- whether WisecondorX predict was rerun after reference change;
- Branch A recall on Y1-Y8 and H1-H6;
- H6 chr21 status;
- per-sample Branch A burden for H7-H16 and all reference-included samples;
- Branch B V2 tier counts by sample;
- whether the report was generated from V2 classifier output;
- number of candidates rescued from hard artifact to review;
- number of candidates downgraded by background compatibility;
- effect on 2026-06-15 five-sample report burden;
- SCA/Branch S status and whether it remains shadow-only.

Promotion is not allowed if any known-positive truth event becomes missing or
hard artifact.

## Current Recommendation

Proceed with document-to-implementation planning for:

```text
reference/Branch A audit and freeze
-> Branch B candidate evidence refinement
-> CNVpro-inspired consistency comparison
-> SCA workflow finalization
```

Do not release legacy current-workflow reports as the forward report path.
Generate future reports only after the Branch A/reference contract and the
Branch B evidence/disposition table are aligned.

Do not promote V2 evidence into unvalidated hard PASS/artifact calls yet.

Do not use legacy Branch B kept count as the final N0 benchmark.

Do not treat complete N0/cross-fit background as an immediate prerequisite for
all development. It remains the clean promotion target; until then,
reference-cohort and N1 evidence are limited context and must be labeled as
such.

Do not reject H7-H16 from reference-rebuild exploration solely because the old
reference produced Branch A signals.

Do require a two-level downstream rerun:

```text
P1/P2 audit:
WisecondorX predict
-> Branch A
-> optional current-code Branch B/S/report diagnostic snapshots

Report-facing release:
WisecondorX predict
-> Branch A
-> fixed Branch B evidence/classification
-> fixed Branch S
-> evaluation / benchmark / report
```

The next concrete step should be an implementation plan that adds
reference-sample Branch A audit as a first-class workflow output, then evaluates
whether reference-cohort context and CNVpro-inspired fields improve Branch B
triage under the fixed Branch A/reference contract. Remote tests and dry-runs
must run on `fengxian`.
