# Branch B Remaining FP Root-Cause Review 2026-06-14

## Scope

This review examines the remaining kept events after the Branch B SCA sex-consistency update.

Remote project path:

`/data/project/CNV/PGT-A/refactor_validation_20260419`

Primary current result files:

- `reports/branch_b_fp_tightening_20260614/all_kept_events.tsv`
- `reports/branch_b_threshold_sensitive_20260614/truth_hit_summary.tsv`
- `reports/branch_b_threshold_sensitive_20260614/kept_event_summary.tsv`

## Current State

Current truth recall:

| cohort | truth events | truth hits | FN | recall |
| --- | ---: | ---: | ---: | ---: |
| H1-H6 | 10 | 10 | 0 | 1.0 |
| Y1-Y8 | 10 | 10 | 0 | 1.0 |

Current H7-H16 negative kept events after SCA sex-consistency:

| sample | event | bins | A abs z | Branch B adjusted z | context |
| --- | --- | ---: | ---: | ---: | --- |
| H8 | chr21 gain 15-42 Mb | 27 | 7.87 | 0.20 | A-sensitive review |
| H13 | chr21 gain 15-42 Mb | 27 | 22.69 | 0.45 | A-strong review |
| H13 | chr12 gain 123-132 Mb | 9 | 18.21 | 0.50 | A-strong edge review |
| H14 | chr21 gain 15-42 Mb | 27 | 8.80 | 0.31 | A-sensitive review |
| H14 | chr15 gain 24-102 Mb | 78 | 22.16 | 0.16 | A-strong broad/edge review |

H7-H16 CNV input QC is PASS for H8/H13/H14:

- H8: total_counts 33,062,516; nonzero_fraction 0.935; mad_log1p 0.166; status PASS
- H13/H14 also have PASS CNV QC files under `results_h_20260608_mask_only/wisecondorx/cnv/qc/`.

## Root-Cause Pattern

The remaining FP are not explained by basic CNV input QC failure. They are Branch A-supported events that Branch B keeps only as review candidates because the current design protects sensitivity.

Key pattern:

1. `chr21:15-42 Mb gain` appears in H8/H13/H14 negatives, but the same event interval is also the best Branch B hit for real positive events:
   - H4 +21 truth
   - Y2 +21 truth
   - H6 focal chr21 gain truth, where the focal truth is contained inside the broad chr21 Branch B event

2. Low clean-bin support is not discriminative:
   - chr21 15-42 Mb has `clean_bin_fraction = 0.0` for both FP and true-positive hits.
   - This comes from region annotation / small chromosome behavior and cannot be used as a hard artifact rule.

3. Low Branch B adjusted z is not discriminative:
   - H8 FP: 0.20
   - H14 FP: 0.31
   - H6 truth chr21 focal hit: 0.33
   - H4/Y2 true +21 hits are also low by Branch B adjusted z, despite strong Branch A evidence.

4. A-branch candidate burden is not enough for hard filtering:
   - H13/H14 have many Branch A candidates, but H1-H6/Y1-Y8 positives also have many Branch A candidates.
   - Median A-branch event burden is similar across H1-H6, H7-H16, and Y1-Y8.

5. Blacklist/gap/repeat overlap is not sufficient for hard filtering:
   - Several true events still carry `blacklist_overlap` or `gap_centromere_telomere_overlap` annotations.
   - These annotations are useful review evidence, but not a safe standalone reject criterion.

## Unsafe Hard Filters

The following possible filters would reduce FP but are unsafe because they would also affect truth hits:

| Candidate filter | Why unsafe |
| --- | --- |
| Hard reject chr21 15-42 Mb gain | Directly overlaps H4/Y2 +21 and H6 focal chr21 truth. |
| Raise A-sensitive review threshold above 9 | Would lose H6 chr21 focal truth with A abs z around 7.39. |
| Require Branch B adjusted z above 0.5 | Would lose H6 chr21, H4 +21, Y2 +21 and other weak/mosaic/broad events. |
| Hard reject clean_bin_fraction = 0 | Would reject true chr21/chr15 events in current truth panels. |
| Hard reject blacklist/gap overlap | Would reject multiple current truth hits. |
| Sample-level A-candidate burden hard reject | Positives have comparable A-candidate burden. |

## Interpretation

The current remaining FP are a reference/batch/recurrent-review problem more than an event-level threshold problem.

Branch B is currently doing the correct conservative thing for the stated clinical strategy:

- preserve Branch A sensitivity
- keep weak but plausible positive regions as review
- reject high-confidence technical artifacts only when evidence is specific enough

Further event-level hard filtering would likely reintroduce FN.

## Recommended Next Step

Do not add another hard Branch B filter for the remaining five events yet.

Recommended direction:

1. Add more true negative samples and measure recurrent review burden by region.
2. Rebuild/evaluate reference with selected H7-H16 normal samples, but exclude samples with unexplained kept review events from the first reference expansion trial.
3. Keep H7, H9, H10, H11, H12, H15, H16 as the first H-batch candidate normals for ref evaluation.
4. Hold H8, H13, H14 out initially because they have retained chr21/chr12/chr15 review events.
5. Continue to report chr21 15-42 Mb gain as candidate_review unless a new reference or larger negative panel proves it is recurrent technical noise without truth-panel overlap.

## Validation Evidence

Commands were run on remote `fengxian` with absolute executables.

Unit test:

```bash
cd /data/project/CNV/PGT-A/refactor_validation_20260419 &&
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python \
-m pytest tests/unit/test_branch_b_artifact_rules.py -q
```

Result:

`38 passed in 0.57s`

Workflow reruns:

- H1-H16: `36 of 36 steps (100%) done`
- Y1-Y8: `22 of 22 steps (100%) done`

Summary evaluation:

- H1-H6 recall: 10/10, FN=0
- Y1-Y8 recall: 10/10, FN=0
- H7-H16 kept events: 5

