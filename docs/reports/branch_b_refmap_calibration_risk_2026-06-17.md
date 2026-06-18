# Branch B Refmap And Calibration Risk Review

Date: 2026-06-17

Context:

- Local repo: `D:\Pipeline\PGT-A`
- Remote mirror used for evidence: `/data/project/CNV/PGT-A/refactor_validation_20260419`
- Active 0615 run: `/data/project/CNV/PGT-A/refactor_validation_20260419/results_20260615_mask_only`
- Current handoff source: `docs/handoff/2026-06-15_0125_truth_eval_scope_handoff.md`

## Summary

The current Branch B report compression on the 2026-06-15 five-sample run is effective, but some hard filters are not yet supported by a complete stability evidence chain. In particular, `cnvseq_subchrom_boundary_weak_support` and direction-discordance filtering are using weak Branch B calibrated evidence to suppress events that can still have strong Branch A evidence.

This creates a possible mosaic or weak-positive false-negative risk.

The main issue is not only threshold strictness. The underlying problem is that event-level WisecondorX ref-bin stability evidence is currently missing, and Branch B calibration falls back to a broad null pool because no bins are marked as `calibration_null_eligible`.

## Observed 0615 Report Compression

The 0615 run contains five samples:

| sample | QC | sex | Branch A candidates | Branch B total | Branch B kept | PASS | review |
|---|---|---:|---:|---:|---:|---:|---:|
| JZ26125843-56-56 | PASS | XX | 29 | 29 | 2 | 0 | 2 |
| JZ26125844-59-59 | PASS | XX | 29 | 29 | 2 | 0 | 2 |
| JZ26125845-60-60 | PASS | XY | 47 | 47 | 20 | 0 | 20 |
| JZ26125846-61-61 | PASS | XY | 38 | 38 | 2 | 0 | 2 |
| JZ26125847-62-62 | PASS | XY | 27 | 27 | 2 | 0 | 2 |

For 56, 59, 61, and 62, Branch B strongly reduces the reportable/review list. The dominant `filter_reason` values are:

| filter reason | count across 56/59/61/62 |
|---|---:|
| `cnvseq_subchrom_boundary_weak_support` | 74 |
| `signal_support_below_minimum` | 73 |
| `chromosome_fraction_too_large` | 42 |
| `branch_b_direction_discordant_with_candidate_state` | 31 |
| `recurrent_artifact_region_weak_support` | 11 |
| `segmental_duplication_overlap_with_limited_clean_support` | 8 |
| `recurrent_edge_lowcal_weak_support` | 7 |
| `bin_count_below_minimum` | 3 |

## Recent Rule Changes Most Responsible

Recent commits show the likely sources of the reduced report count:

| rule | introduced in | status |
|---|---|---|
| `cnvseq_subchrom_boundary_weak_support` | `95bc2d5 Add CNVseq-style Branch B artifact refinement` | new hard filter |
| `branch_b_direction_discordant_with_candidate_state` | `0357957 Refine Branch B direction consistency` | new hard filter |
| `recurrent_artifact_region_weak_support` | `bbe3540 Refine Branch B CNV-seq reporting tiers` | new hard filter |
| `recurrent_edge_lowcal_weak_support` | `bbe3540 Refine Branch B CNV-seq reporting tiers` | new hard filter |
| `signal_support_below_minimum` | present before these commits | older rule, still highly active |
| `chromosome_fraction_too_large` | present before these commits | older rule, still highly active |

The most important new filter is `cnvseq_subchrom_boundary_weak_support`.

## Evidence Chain Gap

### Refmap evidence is missing in the active report

The 0615 final artifact tables contain the expected refmap columns, but the actual count fields are empty:

| column | current state in 0615 artifact events |
|---|---|
| `wisecondorx_ref_bin_count` | all null |
| `reference_bin_count` | all null |
| `ref_bin_count` | all null |
| `ref_size_after_cutoff` | all null |
| `low_refbin_fraction` | present but default 0 |
| `same_chrom_ref_bin_count` | present but default 0 |

The active configs do not provide `reference_assets.refmap_tsvs`, and no current result file under the build-ref or 0615 result tree provides a usable refmap metric table.

This means the report currently has refmap-shaped columns, but it does not have real event-level ref-bin stability evidence.

### Calibration null is not clean

In every 0615 calling bins table:

| field | value |
|---|---:|
| total bins | 3113 |
| `calibration_null_eligible=1` | 0 |
| `region_risk_class=moderate` | 2658 |
| `region_risk_class=high` | 455 |
| `region_risk_class=clean` | 0 |

Because there are no eligible calibration-null bins, `pgta/predict/branch_b/calibration.py` falls back to a broad robust-z quantile pool and then to all bins if needed.

The 0615 calibration summaries show large null scales in the samples where signals are strongly compressed:

| sample | null_bin_count | null_scale |
|---|---:|---:|
| JZ26125843-56-56 | 2489 | 13.4037 |
| JZ26125844-59-59 | 2489 | 15.3189 |
| JZ26125845-60-60 | 2489 | 2.1047 |
| JZ26125846-61-61 | 2489 | 13.6173 |
| JZ26125847-62-62 | 2489 | 11.8160 |

When `null_scale` is this high, Branch B calibrated z can become close to zero even for events with strong Branch A evidence.

## Examples Of Potential Over-Filtering

Several `cnvseq_subchrom_boundary_weak_support` events have strong Branch A z but near-zero Branch B median z:

| sample | event | A z | Branch B median z | filter |
|---|---|---:|---:|---|
| JZ26125846-61-61 | chr12:37.0-132.0Mb gain | 43.54 | -0.019 | `cnvseq_subchrom_boundary_weak_support` |
| JZ26125844-59-59 | chr8:46.0-145.0Mb loss | -42.97 | -0.015 | `cnvseq_subchrom_boundary_weak_support;branch_b_direction_discordant_with_candidate_state` |
| JZ26125846-61-61 | chr13:54.0-112.0Mb loss | -35.40 | -0.034 | `cnvseq_subchrom_boundary_weak_support` |
| JZ26125847-62-62 | chr7:63.0-76.0Mb gain | 30.53 | -0.115 | `segmental_duplication_overlap_with_limited_clean_support;cnvseq_subchrom_boundary_weak_support` |
| JZ26125844-59-59 | chr9:70.0-141.0Mb gain | 27.19 | 0.009 | `cnvseq_subchrom_boundary_weak_support` |

These rows do not prove true positives. They do show that the current hard filter can suppress strong Branch A candidates using a Branch B signal that may be unstable or over-normalized.

## Risk Assessment

The current behavior is acceptable as an FP compression experiment, but it should not yet be treated as a final clinical-style rule set.

Primary risk:

- weak positive or mosaic events can be filtered if they are large subchromosomal events near a high-risk boundary and Branch B calibrated z is compressed by a large null scale.

The risk is highest when:

- A has moderate or strong z;
- Branch B median or calibrated z is near zero;
- event spans a centromere/telomere/gap-adjacent transition;
- refmap stability evidence is missing;
- calibration null has no clean eligible bins;
- event is low-fraction mosaic.

## Recommendations

### Immediate behavior change to evaluate

Do not keep `cnvseq_subchrom_boundary_weak_support` as a hard filter while refmap and calibration-null evidence are missing.

Recommended temporary behavior:

- downgrade to review-low;
- keep the flag and downgrade reason;
- do not use it alone to set `keep_event=0` when Branch A support is strong or sensitive-review eligible.

The same caution applies to direction-discordance when Branch A is strong and Branch B is calibrated against a large fallback null scale.

### Refmap contract fix

Missing refmap metrics should remain unknown, not silently become benign.

Recommended changes:

- keep missing `wisecondorx_ref_bin_count` as NaN;
- keep missing `low_refbin_fraction` as NaN or `unknown`, not 0;
- keep missing `same_chrom_ref_bin_count` as NaN or `unknown`, not 0;
- add a report-level flag such as `refmap_metrics_available=0`;
- do not use absent refmap columns as evidence that a region is stable.

### Calibration-null fix

The current `calibration_null_eligible=0` for all bins means Branch B null calibration is not following the intended clean-bin contract.

Recommended changes:

- define a true clean autosomal null set;
- do not block all pass/moderate bins from calibration if no clean bins exist;
- expose `null_selection_mode` in each summary:
  - `clean_neutral`
  - `eligible_quantile`
  - `all_bins_fallback`
- if fallback is used, report-level filters should be less aggressive.

### Validation matrix before promoting hard filters

Run an ablation matrix on:

- Y1-Y8 truth positives;
- H1-H6 truth positives;
- H7-H16 negatives;
- 2026-06-15 five-sample run.

Suggested variants:

| variant | change |
|---|---|
| baseline_current | current rules |
| boundary_review_only | `cnvseq_subchrom_boundary_weak_support` becomes review-only |
| direction_protect_30 | lower A-protect threshold for direction discordance from 50 to 30 |
| recurrent_review_only | recurrent artifact filters become review-only |
| repaired_calibration_null | fix null eligibility and rerun current filters |
| repaired_refmap_plus_null | add real refmap metrics and repaired null calibration |

Primary acceptance criteria:

- no reduction in known truth recall;
- no known mosaic/weak positive loss;
- report events increase is quantified per sample;
- 0615 suspected negatives do not explode into unmanageable reports;
- every hard filter has an explicit evidence source in the output table.

## Evidence Commands

All result evidence was collected on `ssh fengxian`.

Commands used:

```bash
cd /data/project/CNV/PGT-A/refactor_validation_20260419

/biosoftware/miniconda/envs/snakemake_env/bin/python - <<'PY'
from pathlib import Path
import pandas as pd
base=Path('results_20260615_mask_only/wisecondorx/cnv/postprocess/artifact_rules')
for f in sorted(base.glob('*.candidate_events.tsv')):
    df=pd.read_csv(f, sep='\t')
    cols=[c for c in [
        'wisecondorx_ref_bin_count',
        'reference_bin_count',
        'ref_bin_count',
        'ref_size_after_cutoff',
        'low_refbin_fraction',
        'same_chrom_ref_bin_count',
    ] if c in df.columns]
    print(f.name, {c:int(pd.to_numeric(df[c], errors='coerce').notna().sum()) for c in cols})
PY
```

Result:

- all event-level ref-bin count fields had 0 non-null values;
- `low_refbin_fraction` and `same_chrom_ref_bin_count` were present, but defaulted to 0.

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/python - <<'PY'
from pathlib import Path
import pandas as pd
base=Path('results_20260615_mask_only/wisecondorx/cnv/postprocess/calling')
for f in sorted(base.glob('*.bins.tsv')):
    df=pd.read_csv(f, sep='\t')
    print(f.name, df['calibration_null_eligible'].value_counts(dropna=False).to_dict())
PY
```

Result:

- every 0615 sample had `calibration_null_eligible={0: 3113}`.

Calibration summary files reviewed:

- `results_20260615_mask_only/wisecondorx/cnv/postprocess/calibration/*.summary.json`

## Conclusion

The current Branch B filters are useful for identifying likely artifacts, but several hard-filter decisions are currently ahead of the evidence chain. Before using them as final report-suppression rules, the workflow should repair refmap propagation and calibration-null selection, then rerun an ablation matrix against known positives and the 0615 run.
