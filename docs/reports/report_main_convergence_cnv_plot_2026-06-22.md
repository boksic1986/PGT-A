# Report Main Convergence And CNV Plot Contract

Date: 2026-06-22

## Scope

This loop updates the report layer and CNV plot output only.

Unchanged:

- active 1Mb reference: `h_r0_shadow_ref_20260619`
- Branch A caller: WisecondorX/CBS with explicit `merge_gap_bp=2_000_000` overlay in the active validation configs
- Branch B V2: Branch-A-anchored report-layer evidence, no B-only events
- Branch S: review/development-only, no final SCA promotion
- reference build and BWA/mapping

## Implemented Contract

### Branch B V2 report visibility

Branch B V2 now exposes report visibility labels:

- `report_strong_event`
- `report_weak_event`
- `internal_review_event`
- `filtered_event`
- `branch_s_event`

Final autosomal CNV report event input is the V2 benchmark `report_events.tsv`,
which contains only `report_strong_event` and `report_weak_event` rows.
Internal review rows remain in V2/sample summaries and audit ledgers.
Filtered rows remain audit-only. Branch S rows remain in the SCA/sex-chromosome
section and do not enter the autosomal CNV main table.

The weak-report rule is conservative: a sensitive A signal is report-weak only
when it has same-direction B-side support, clean support, and same-direction
lowres support. This protects H6 chr21 while avoiding promotion of every
background-unknown review row to the main table.

### Report sample universe

`cnv_report_summary` now consumes Branch B V2 benchmark summary and sample
summary whenever `branch_b_v2_benchmark` is available in the workflow. It no
longer depends on `branch_b_v2_benchmark` being explicitly requested by the
user.

The report sample table is expanded from the V2 sample summary, so samples with
zero autosomal final report events remain visible with event counts set to zero.
This fixed the prior report-layer gap where Y3 and several H samples were absent
from `cnv_summary.tsv` when refreshing the report alone.

### WisecondorX-style CNV plot

`pgta/predict/branch_b/plot.py` now uses only `calibrated_z` as the bin-level
plot signal.

Rules:

- missing `calibrated_z` column is an error
- non-numeric or non-finite `calibrated_z` bins are skipped
- no fallback to `robust_z`, `raw_robust_z`, or normalized signal
- per-sample plot table is written to
  `wisecondorx/cnv/plots/{sample}.plot_bins.tsv`
- plot TSV columns include `chrom`, `start`, `end`, `genome_pos`, `z`,
  `report_state`
- `report_state` is limited to `dup`, `del`, `neutral`
- only final autosomal report events are highlighted
- SVG legend is limited to `dup`, `del`, `neutral bin`, `report z trend`
- `dup` bins use WisecondorX-style yellow
- `del` bins use blue
- no genome-wide smooth polyline is drawn
- each final autosomal report event gets one red horizontal report-z trend line,
  spanning only that event interval and using the event median `calibrated_z`

### Event-level copy-number plot supplement

Each sample now also gets a copy-number visualization beside the calibrated-z
plot.

New outputs:

- `wisecondorx/cnv/plots/{sample}.final_cnv_cn.svg`
- `wisecondorx/cnv/plots/{sample}.plot_bins_cn.tsv`

CN plot contract:

- neutral bins are displayed as `CN=2`
- bins overlapping final autosomal report events inherit that event's
  `copy_number_estimate`
- if `copy_number_estimate` is absent, the plotter can fall back to
  `sex_adjusted_copy_number`, then `a_ratio` using `CN = 2 * (1 + a_ratio)`
- if all event-level CN sources are absent, CN plotting fails instead of
  fabricating CN from `calibrated_z` or `normalized_signal`
- the y axis is `Copy number`
- reference lines are `CN=1`, `CN=2`, and `CN=3`
- `dup` is yellow, `del` is blue, neutral bins are grey, and event CN trend
  lines are red
- no genome-wide smooth CN line is drawn
- the SVG legend is limited to `dup`, `del`, `neutral bin`, and
  `report CN trend`

The CN plot is an event-level report visualization. It is not a bin-level CN
caller and does not change Branch A, Branch B V2, Branch S, report visibility,
or filtering.

## Remote Validation

Remote path:
`/data/project/CNV/PGT-A/refactor_validation_20260419`

Executables:

- python: `/biosoftware/miniconda/envs/snakemake_env/bin/python`
- snakemake: `/biosoftware/miniconda/envs/snakemake_env/bin/snakemake`

Unit tests:

```text
tests/unit/test_cnv_report.py
tests/unit/test_branch_ab_phase12_workflow_contract.py
tests/unit/test_branch_b_plot.py
tests/unit/test_branch_b_v2_classifier.py
tests/unit/test_branch_b_v2_benchmark.py
```

Result:

```text
79 passed
```

Dry-run:

```text
branch_b_v2_benchmark branch_s_review cnv_report -n
```

All four active lowres-enabled gap2m configs parsed successfully:

- `config_predict_y_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`
- `config_predict_h_20260608_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`
- `config_predict_g1_g8_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`
- `config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`

Report refresh:

```text
snakemake --forcerun cnv_report_summary cnv_report
```

This was required because `cnv_report` is an aggregate target; forcing only the
aggregate target does not rerun `report.py`.

## Materialized Result Summary

| cohort | samples | candidates | truth | preserved | FN | truth filtered | report | strong | weak | internal | filtered | Branch S | SVG | plot TSV | report rows |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Y1-Y8 | 8 | 97 | 10 | 10 | 0 | 0 | 40 | 21 | 19 | 31 | 13 | 13 | 8 | 8 | 8 |
| H1-H16 | 16 | 105 | 10 | 10 | 0 | 0 | 23 | 6 | 17 | 32 | 8 | 42 | 16 | 16 | 16 |
| G1-G8 | 8 | 75 | 10 | 10 | 0 | 0 | 26 | 15 | 11 | 29 | 7 | 13 | 8 | 8 | 8 |
| 2026-06-15 | 5 | 165 | 0 | 0 | 0 | 0 | 71 | 52 | 19 | 57 | 23 | 14 | 5 | 5 | 5 |

Key locked truth checks:

- Y1-Y8: truth preserved 10/10, FN=0, truth filtered=0.
- H1-H6: truth preserved 10/10, FN=0, truth filtered=0.
- H6 chr21 gain remains `report_weak_event`,
  `keep_report_layer_event`, A z = 7.1135.
- H6 chrX remains routed to Branch S review.
- G1-G8: truth preserved 10/10, FN=0, truth filtered=0.
- G2 chr8 gain remains visible as `internal_review_event`; it is not filtered.
- G2 chr11 loss remains `report_strong_event`.
- 2026-06-15 has no locked truth and remains burden/context only.

Plot contract checks:

- all cohorts have matching `.final_cnv.svg` and `.plot_bins.tsv` per sample
- all plot TSV files contain required columns
- `report_state` values are restricted to `dup`, `del`, `neutral`
- SVG legend contains only `dup`, `del`, `neutral bin`, `report z trend`
- SVG does not contain Branch A/B, internal, filtered, mask, or rejected legend
  terms

## 2026-06-22 Plot Style Refinement

This follow-up changed only CNV plot rendering. Branch A, Branch B V2, Branch S,
report classification, and reference outputs were unchanged.

Updated rendering:

- removed chromosome/global smoothed polylines
- added red horizontal `report-z-trend` lines only over final autosomal report
  event intervals
- changed duplication highlight color to WisecondorX-style yellow
- kept deletion highlight color blue
- kept neutral bins grey

Remote validation:

```text
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_plot.py \
  tests/unit/test_cnv_report.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py -q
```

Result:

```text
32 passed
```

0615 plot rematerialization:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml \
  --cores 4 --forcerun cnv_branch_ab_plot cnv_report
```

Result:

```text
9 of 9 steps (100%) done
Complete log: .snakemake/log/2026-06-22T093707.721304.snakemake.log
```

0615 SVG content check after local sync:

| sample | report trend lines | polyline present | smooth legend present | dup yellow | del blue | red trend |
|---|---:|---|---|---|---|---|
| JZ26125843-56-56 | 7 | false | false | true | true | true |
| JZ26125844-59-59 | 11 | false | false | true | true | true |
| JZ26125845-60-60 | 23 | false | false | true | true | true |
| JZ26125846-61-61 | 18 | false | false | true | true | true |
| JZ26125847-62-62 | 12 | false | false | true | true | true |

Local synced SVG path:

`D:\Pipeline\PGT-A\reports\0615_cnv_plots\*.final_cnv.svg`

## 2026-06-22 Copy Number Plot Supplement

This follow-up added event-level copy-number plots while preserving the existing
calibrated-z plot.

Remote unit tests:

```text
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_plot.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py \
  tests/unit/test_cnv_report.py -q
```

Result:

```text
33 passed in 1.07s
```

0615 dry-run:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml \
  --cores 1 -n --forcerun cnv_branch_ab_plot cnv_report
```

Result:

```text
DAG planned only 5 cnv_branch_ab_plot jobs, cnv_report_summary, cnv_report,
collect_runtime_tracking, and all. No mapping or reference build jobs were
requested.
```

0615 materialization:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml \
  --cores 4 --forcerun cnv_branch_ab_plot cnv_report
```

Result:

```text
9 of 9 steps (100%) done
Complete log: .snakemake/log/2026-06-22T102013.293346.snakemake.log
```

0615 CN SVG content check:

| sample | CN trend lines | polyline present | copy-number label | calibrated-z label | dup yellow | del blue | red trend |
|---|---:|---|---|---|---|---|---|
| JZ26125843-56-56 | 7 | false | true | false | true | true | true |
| JZ26125844-59-59 | 11 | false | true | false | true | true | true |
| JZ26125845-60-60 | 23 | false | true | false | true | true | true |
| JZ26125846-61-61 | 18 | false | true | false | true | true | true |
| JZ26125847-62-62 | 12 | false | true | false | true | true | true |

0615 CN bin table checks:

| sample | CN bins | non-neutral CN bins | CN source |
|---|---:|---:|---|
| JZ26125843-56-56 | 4144 | 503 | `copy_number_estimate` plus neutral baseline |
| JZ26125844-59-59 | 4144 | 923 | `copy_number_estimate` plus neutral baseline |
| JZ26125845-60-60 | 4144 | 1274 | `copy_number_estimate` plus neutral baseline |
| JZ26125846-61-61 | 4144 | 1193 | `copy_number_estimate` plus neutral baseline |
| JZ26125847-62-62 | 4144 | 1016 | `copy_number_estimate` plus neutral baseline |

The regenerated `cnv_summary.tsv` contains `copy_number_plot_svg` and links all
5/5 0615 samples to `.final_cnv_cn.svg`.

Local synced CN outputs:

- `D:\Pipeline\PGT-A\reports\0615_cnv_plots\*.final_cnv_cn.svg`
- `D:\Pipeline\PGT-A\reports\0615_cnv_plots\*.plot_bins_cn.tsv`

## Interpretation

This loop makes the report package more usable:

- autosomal final report main table now contains only V2 report-layer events
- strong and weak report events are both visible
- internal review and filtered events are separated from the main report table
- Branch S remains separate from autosomal CNV
- zero-event samples remain present in the sample report summary
- each sample has a WisecondorX-style calibrated-z CNV plot and matching bin
  table

This is still a development report contract. It does not promote the shadow
reference, Branch B V2, Branch S, or lowres evidence to production-final.

## Remaining Work

The 2026-06-15 burden remains high, especially in samples 60 and 61. Further
reduction should not start with another broad threshold. The next gate should
inspect candidate-level evidence among remaining `report_strong_event` and
`report_weak_event` rows:

- same-direction lowres support
- ref-MAD stability
- clean support and region-risk
- B-side signal context
- SCA/sex-chromosome separation
- truth-protecting ablation across Y/H/G

Any future rule must keep Y/H/G FN=0 and must not filter H6 chr21 or G2 locked
truth.

## 2026-06-22 CN Plot V2 Update

This update changes only the copy-number visualization. Branch A, Branch B V2,
Branch S, report-event classification, filtering, mapping, and reference build
are unchanged.

### CN V2 contract

- z plot remains unchanged and still uses only `calibrated_z`.
- CN plot width is increased to `2560`.
- CN plot uses a dark plotting background, wider chromosome gaps, and 50Mb
  intra-chromosome ticks.
- structural gap / centromere-telomere bins are blanked when
  `is_gap_centromere_telomere` is true or
  `gap_centromere_telomere_overlap_fraction >= 0.5`.
- report CN trend remains a red horizontal event-level line and is split across
  non-gap contiguous chunks.
- CN scatter points are drawn per bin:
  - neutral outside final events: `CN=2`
  - final `dup`: dark blue
  - final `del`: red
- event-bin scatter uses an event-anchored CN proxy:
  `2 + calibrated_z * (event_cn - 2) / median_event_z` when median z is
  informative and direction-consistent.
- if event median z is uninformative, event bins fall back to uniform event CN.
- the CN proxy is for visual review only and is not an independent CN caller.

### Remote validation

TDD red run:

```text
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_plot.py -q
```

Result on old implementation:

```text
2 failed, 2 passed
```

Final unit tests:

```text
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_plot.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py \
  tests/unit/test_cnv_report.py -q
```

Result:

```text
34 passed in 0.89s
```

0615 dry-run:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml \
  --cores 1 -n --forcerun cnv_branch_ab_plot cnv_report
```

Result:

```text
Only 5 cnv_branch_ab_plot jobs, cnv_report_summary, cnv_report,
collect_runtime_tracking, and all were planned. No mapping or reference rebuild
jobs were requested.
```

0615 materialization:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml \
  --cores 4 --forcerun cnv_branch_ab_plot cnv_report
```

Result:

```text
9 of 9 steps (100%) done
Complete log: .snakemake/log/2026-06-22T105631.231123.snakemake.log
```

### 0615 CN V2 output check

| sample | width 2560 | trend chunks | gap blanks | 50Mb ticks | dup/del-only legend | polyline |
|---|---|---:|---:|---|---|---|
| JZ26125843-56-56 | true | 15 | 544 | true | true | false |
| JZ26125844-59-59 | true | 26 | 544 | true | true | false |
| JZ26125845-60-60 | true | 47 | 544 | true | true | false |
| JZ26125846-61-61 | true | 42 | 544 | true | true | false |
| JZ26125847-62-62 | true | 34 | 544 | true | true | false |

| sample | CN TSV rows | structural-gap blank bins | event-scaled proxy bins |
|---|---:|---:|---:|
| JZ26125843-56-56 | 4144 | 544 | 0 |
| JZ26125844-59-59 | 4144 | 544 | 0 |
| JZ26125845-60-60 | 4144 | 544 | 1003 |
| JZ26125846-61-61 | 4144 | 544 | 0 |
| JZ26125847-62-62 | 4144 | 544 | 0 |

Local synced CN V2 outputs:

- `D:\Pipeline\PGT-A\reports\0615_cnv_plots\*.final_cnv_cn.svg`
- `D:\Pipeline\PGT-A\reports\0615_cnv_plots\*.plot_bins_cn.tsv`
