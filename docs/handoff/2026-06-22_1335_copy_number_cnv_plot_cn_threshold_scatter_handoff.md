# Copy Number CNV Plot CN-Threshold Scatter Handoff

Date: 2026-06-22 13:35

## Scope

This handoff supersedes
`docs/handoff/2026-06-22_1305_copy_number_cnv_plot_bin_z_scatter_handoff.md`
for copy-number plot interpretation.

The change is visualization-only:

- no Branch A change;
- no Branch B V2 filtering or report-event classification change;
- no Branch S/SCA change;
- no reference, mapping, or WisecondorX predict change.

## Updated CN Scatter Contract

The CN plot now draws a scatter point for every non-centromere 1Mb bin.

Scatter color is based only on the bin-level visual CN proxy, not on the merged
report-event interval:

```text
CN_proxy = clip(2 + calibrated_z * 0.05, 0, 4)
CN < 1.7          -> del, red #ef4444
1.7 <= CN <= 2.3 -> neutral, grey #64748b
CN > 2.3          -> dup, blue #1d4ed8
```

The 1.7/2.3 thresholds correspond to a 30% mosaic review band around diploid:

```text
30% loss CN = 1.7
30% gain CN = 2.3
```

The output CN TSV now includes:

```text
chrom
start
end
genome_pos
z
copy_number
report_state
event_report_state
copy_number_source
```

`report_state` is the bin-level color state. `event_report_state` is only an
audit field that records overlap with final autosomal report events.

The horizontal CN trend lines still represent event-level copy-number estimates
from the final report events and remain separate from bin scatter. They do not
control scatter color.

## Large CN Value Explanation

The previously observed very large CN values were produced by the now-retired
event-anchored visualization formula:

```text
2 + calibrated_z * abs(event_cn - 2) / max(median(abs(event_z)), 0.25)
```

Example from sample `JZ26125843-56-56`:

```text
chr11:96.00-96.75Mb
calibrated_z = 37.920977
old visual CN proxy = 14.923469
```

That value was not a real copy-number estimate. It reflected an extreme
calibrated-z bin being amplified by the event-level scaling denominator.

Current output clips the visual proxy to `0..4` and records
`copy_number_source=calibrated_z_mosaic30_cn_proxy`.

## Files Modified

- `pgta/predict/branch_b/plot.py`
- `tests/unit/test_branch_b_plot.py`
- `tests/unit/test_current_context_index.py`
- `docs/CURRENT_CONTEXT_INDEX.md`
- `CURRENT_STATE.md`
- `PLANS.md`
- `docs/reports/report_main_convergence_cnv_plot_2026-06-22.md`
- this handoff

## Remote Validation

Remote mirror:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419
```

Executables:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/python
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake
```

TDD red run on the old implementation:

```text
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_plot.py -q
```

Result:

```text
5 failed, 2 passed
```

Expected failures confirmed old behavior:

- missing `event_report_state`;
- scatter radius still `2.00`;
- no-event scatter not colored by CN proxy;
- event bins still used event-anchored CN scaling.

Final remote unit tests:

```text
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_plot.py \
  tests/unit/test_cnv_report.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py \
  tests/unit/test_current_context_index.py -q
```

Result:

```text
40 passed in 1.19s
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
collect_runtime_tracking, and all were planned. No mapping or reference build
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
Complete log: .snakemake/log/2026-06-22T132746.835251.snakemake.log
```

## Materialized 0615 Check

All five 0615 samples have updated CN plots and CN TSVs.

Per sample:

```text
CN SVG scatter points: 4024
CN scatter radius: 1.35
centromere blanks: 120
chromosome separators: 24
chrom-background rectangles: absent
```

CN TSV sources:

```text
calibrated_z_mosaic30_cn_proxy: 4024
structure_gap_blank: 120
```

Current CN proxy ranges and bin color counts:

| sample | CN range | bin states |
|---|---|---|
| JZ26125843-56-56 | 1.894-4.000 | neutral=3678, dup=346 |
| JZ26125844-59-59 | 1.911-4.000 | neutral=3679, dup=345 |
| JZ26125845-60-60 | 0.000-2.188 | neutral=3766, del=258 |
| JZ26125846-61-61 | 1.823-4.000 | neutral=3682, dup=342 |
| JZ26125847-62-62 | 1.891-4.000 | neutral=3679, dup=345 |

Local synced outputs:

```text
D:\Pipeline\PGT-A\reports\0615_cnv_plots\*.final_cnv_cn.svg
D:\Pipeline\PGT-A\reports\0615_cnv_plots\*.plot_bins_cn.tsv
```

## Current Interpretation

Use the CN plot for visual review only.

- A blue/red scatter cluster means many bins cross the visual 30% mosaic CN
  proxy threshold.
- A horizontal CN trend line means a final report event has an event-level CN
  estimate.
- Scatter color and trend line are intentionally separate.
- Neither the scatter proxy nor the trend line changes Branch B V2 filtering or
  final report-event status.

## Next Recommended Step

Continue manual review of 0615 samples, starting with sample 56, using:

- the calibrated-z plot for signed z evidence;
- the CN threshold scatter for bin-level coherence;
- the event-level horizontal CN trend for merged report-event amplitude.

Any future candidate-level filtering proposal must still be tested against
Y1-Y8, H1-H6, and G1-G8 locked truth with FN=0, H6 chr21 retained, and G2 truth
visible.
