# Branch S / SCA Validation Plan

Date: 2026-06-19

## Scope

This document defines how Branch S should be validated before any sex-chromosome abnormality evidence is allowed to influence reportable CNV calls.

Branch S remains a shadow-only model. It must not replace:

- current sex calling;
- WisecondorX predict;
- Branch A candidate generation;
- Branch B final-event/report logic;
- `cnv_report`.

## Current Evidence

Remote evidence source:

```text
ssh fengxian
cd /data/project/CNV/PGT-A/refactor_validation_20260419
```

Current materialized Branch S output exists for:

```text
results_build_ref_v2_mask_only/wisecondorx/cnv/postprocess/branch_s
results_h_20260608_mask_only/wisecondorx/cnv/postprocess/branch_s
results_20260615_mask_only/wisecondorx/cnv/postprocess/branch_s
```

Existing Branch S files:

- Y1-Y8: 8 sex-chromosome evidence tables, 8 SCA state score tables, 8 summary JSON files;
- H1-H16: 16 sex-chromosome evidence tables, 16 SCA state score tables, 16 summary JSON files;
- 2026-06-15: 5 sex-chromosome evidence tables, 5 SCA state score tables, 5 summary JSON files;
- every row has `branch_s_action=SHADOW_ONLY`;
- every row has `final_report_impact=none_shadow_only`.

This proves the workflow contract is report-safe for the materialized Y1-Y8, H1-H16, and 2026-06-15 runs, but it does not prove SCA clinical/reporting performance.

## Current SCA Truth Coverage

Current truth events with sex-chromosome involvement:

| sample | truth | current limitation |
|---|---|---|
| Y3 | `chrX loss`, `45,XO` | Branch S output exists, but only one old-batch XO-like case. |
| H5 | three `chrX loss` segments, including one mosaic segment | Branch S output exists, but score direction is not validated. |
| H6 | whole-`chrX loss`, `45,XO` inference | Branch S output exists, but score direction is not validated. |

Current Branch S score-direction caveat:

| run | sample | truth context | X non-PAR mean calibrated z | X Branch A candidate count | current X scores |
|---|---|---|---:|---:|---|
| Y1-Y8 | Y3 | chrX loss / XO-like | 2.4326 | 2 | `X_GAIN=2.433`, `X_LOSS=-2.433` |
| H1-H16 | H5 | chrX loss segments, including mosaic | 2.1368 | 2 | `X_GAIN=2.137`, `X_LOSS=-2.137` |
| H1-H16 | H6 | whole-chrX loss / XO-like | 1.9344 | 2 | `X_GAIN=1.934`, `X_LOSS=-1.934` |

The corresponding Branch A chrX candidates are all loss calls with strong z support. Therefore current Branch S state scores must not be interpreted as gain/loss classification. At this stage Branch S can provide region-level evidence and review context only.

Missing locked truth classes:

- `chrX gain` / XXX-like samples;
- XY samples with extra X / XXY-like samples;
- Y gain / XYY-like samples;
- Y loss or structurally abnormal Y samples;
- mosaic SCA with known fraction series;
- clean XX and XY negatives across multiple library batches;
- PAR / XY-homology edge cases with orthogonal truth.

Conclusion: the current SCA truth set is too narrow and loss-heavy. It can guide review design, but it is not sufficient for promotion.

## Review Rubric

Branch S may assign review labels in shadow output, but these labels must not be consumed by `cnv_report` until promotion gates are met.

Allowed shadow review labels:

| label | meaning | report impact |
|---|---|---|
| `SCA_REVIEW_STRONG` | Broad non-PAR X/Y dosage signal, consistent direction, Branch A/group support present, and no hard contradiction from current sex call. | None. |
| `SCA_REVIEW_WEAK` | Partial, mosaic-like, or lower-amplitude signal that may still be biologically meaningful. | None. |
| `SCA_DISCORDANT_REVIEW` | Sex-chromosome dosage signal conflicts with current sex call or expected sex-aware baseline. | None. |
| `SCA_ARTIFACT_REVIEW` | Signal is limited to PAR/edge/high-risk regions or has technical contradiction, but evidence is not strong enough for automatic rejection. | None. |
| `SCA_NO_CALL` | Missing bins, missing Branch B calibrated bins, missing gender evidence, QC failure, or other unanalyzable condition. | None. |

Do not add a `PASS` or `CONFIRMED` Branch S label until a locked SCA validation set exists.

## Branch S Metrics Required Before Promotion

At minimum, promotion review must report:

- sample-level SCA recall and specificity;
- event-level overlap for X gain, X loss, Y gain, and Y loss;
- sex-call concordance before and after Branch S review;
- X non-PAR, X PAR, Y non-PAR, and Y PAR signal separation;
- PAR/XY-homology artifact burden;
- batch-stratified FP burden for clean XX and XY negatives;
- mosaic fraction sensitivity when dilution or orthogonal truth is available;
- whether Branch S changes any current final report event.

Promotion can only be considered if:

1. Y1-Y8 and H1-H6 known-positive recall does not regress.
2. H5 and H6 chrX loss signals are detected in the formal workflow with score direction that agrees with the locked truth definition.
3. No clean XX/XY negative is converted into a reportable SCA call by Branch S.
4. Branch S decisions are reproducible across batches and not tuned on the 2026-06-15 exploratory samples.
5. Branch S remains explainable at the region level: X non-PAR, X PAR, Y non-PAR, and Y PAR evidence must be visible.

## Report Release Boundary For The 2026-06-15 Five Samples

The 2026-06-15 five-sample run can be reported only as a current-workflow result if all of the following are true:

1. The report is generated by the frozen production/report workflow, not by a hand-edited event table.
2. WisecondorX predict, Branch A candidates, legacy Branch B final events, and `cnv_report` use the same reference ID and config.
3. Every Branch A candidate has a recorded Branch B or Branch B V2 shadow disposition.
4. Missing refmap, missing calibration null, or missing Branch S evidence is shown as `UNKNOWN` / review context, not clean support.
5. Branch S is explicitly marked as shadow-only and is not used to add or remove final report events.

The five-sample run must not be used as a locked validation set because it has already been used during rule discussion and exploratory review.

## Next Minimum Workflow Actions

1. Rework or explicitly redefine Branch S state-score direction so chrX-loss truth is not represented as positive `X_GAIN` evidence.
2. Keep H7-H16 as N1/N2 review labels; do not use them as clean SCA negatives.
3. Collect locked clean XX and XY negatives before estimating SCA specificity.
4. Collect additional SCA positives before promotion, especially X gain / XXY / XXX / XYY / Y-loss and mosaic cases.
5. Repeat formal materialization and truth summary after the score-direction contract is fixed.
