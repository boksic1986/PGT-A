# Handoff: Copy Number CNV Plot Supplement

Date: 2026-06-22 10:25

## 1. Goal

Add a copy-number CNV plot for each sample beside the existing calibrated-z
WisecondorX-style plot.

This is visualization-only. It does not change Branch A, Branch B V2, Branch S,
report-event classification, filtering, reference build, or mapping.

## 2. Confirmed Project State

- Active reference remains `h_r0_shadow_ref_20260619`.
- Active 0615 config remains
  `config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`.
- 0615 has no locked truth labels, so no TP/FP/FN conclusion is allowed.
- Existing calibrated-z plot remains the primary signal plot.

## 3. Completed Changes

Modified:

- `pgta/predict/branch_b/plot.py`
- `pgta/predict/report.py`
- `rules/predict_layout.smk`
- `rules/predict_workflow.smk`
- `rules/target_assembly.smk`
- `tests/unit/test_branch_b_plot.py`
- `tests/unit/test_branch_ab_phase12_workflow_contract.py`
- `tests/unit/test_cnv_report.py`
- `docs/reports/report_main_convergence_cnv_plot_2026-06-22.md`
- `CURRENT_STATE.md`
- `PLANS.md`
- `docs/CURRENT_CONTEXT_INDEX.md`

New plot outputs:

- `{sample}.final_cnv_cn.svg`
- `{sample}.plot_bins_cn.tsv`

CN contract:

- neutral bins are displayed as `CN=2`
- final autosomal report-event bins use event-level `copy_number_estimate`
- fallback order is `copy_number_estimate`, `sex_adjusted_copy_number`,
  then `CN = 2 * (1 + a_ratio)`
- if all event-level CN sources are absent, plotting fails
- `calibrated_z` and `normalized_signal` are not used to compute CN
- CN plot legend is limited to `dup`, `del`, `neutral bin`, and
  `report CN trend`
- no genome-wide smooth CN line is drawn

Report integration:

- `cnv_branch_ab_plot` now emits z SVG/TSV and CN SVG/TSV.
- `cnv_report_summary` receives CN SVGs and writes
  `copy_number_plot_svg` in `cnv_summary.tsv/json`.
- Markdown/HTML sample tables include a `CN Plot` link.

## 4. Remote Validation

Remote path:

`/data/project/CNV/PGT-A/refactor_validation_20260419`

Executables:

- python: `/biosoftware/miniconda/envs/snakemake_env/bin/python`
- snakemake: `/biosoftware/miniconda/envs/snakemake_env/bin/snakemake`

Unit tests:

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

Dry-run:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml \
  --cores 1 -n --forcerun cnv_branch_ab_plot cnv_report
```

Result:

```text
exit 0; planned 5 cnv_branch_ab_plot jobs, cnv_report_summary, cnv_report,
collect_runtime_tracking, and all. No mapping or reference rebuild jobs were
requested.
```

Materialization:

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

## 5. Materialized Output Check

0615 outputs:

- CN SVG count: `5`
- CN bin table count: `5`
- z SVG count: `5`
- `cnv_summary.tsv` contains `copy_number_plot_svg`
- 5/5 report rows link to `.final_cnv_cn.svg`

SVG checks:

| sample | CN trend lines | polyline present | copy-number label | calibrated-z label | dup yellow | del blue | red trend |
|---|---:|---|---|---|---|---|---|
| JZ26125843-56-56 | 7 | false | true | false | true | true | true |
| JZ26125844-59-59 | 11 | false | true | false | true | true | true |
| JZ26125845-60-60 | 23 | false | true | false | true | true | true |
| JZ26125846-61-61 | 18 | false | true | false | true | true | true |
| JZ26125847-62-62 | 12 | false | true | false | true | true | true |

CN bin table checks:

| sample | bins | non-neutral CN bins | CN source |
|---|---:|---:|---|
| JZ26125843-56-56 | 4144 | 503 | `copy_number_estimate` plus neutral baseline |
| JZ26125844-59-59 | 4144 | 923 | `copy_number_estimate` plus neutral baseline |
| JZ26125845-60-60 | 4144 | 1274 | `copy_number_estimate` plus neutral baseline |
| JZ26125846-61-61 | 4144 | 1193 | `copy_number_estimate` plus neutral baseline |
| JZ26125847-62-62 | 4144 | 1016 | `copy_number_estimate` plus neutral baseline |

Local synced outputs:

- `D:\Pipeline\PGT-A\reports\0615_cnv_plots\*.final_cnv_cn.svg`
- `D:\Pipeline\PGT-A\reports\0615_cnv_plots\*.plot_bins_cn.tsv`

## 6. Current Conclusion

The report workflow now generates both:

- calibrated-z plot: bin-level `calibrated_z` points plus final-event red
  horizontal z trends
- copy-number plot: event-level CN blocks plus final-event red horizontal CN
  trends

The CN plot is easier for report review, but it is derived from final report
event CN estimates. It must not be used as independent evidence to create,
filter, or promote CNV calls.

## 7. Suggested Next Step

Use the synced 0615 z/CN plot pairs for manual sample-by-sample review,
starting with sample 56 if continuing the prior review sequence.

If report burden is further reduced, candidate-level evidence must still be
ablated against Y/H/G locked truth with FN=0, H6 chr21 retained, and G2 truth
visible.

## 8. Core File Sync

- `CURRENT_STATE.md`: updated with CN plot contract and validation.
- `PLANS.md`: updated to include CN plot in manual review while keeping it out
  of filter logic.
- `docs/CURRENT_CONTEXT_INDEX.md`: updated to point to this handoff.
- `REPO_MAP.md`: not updated; no stable entrypoint or directory structure
  changed.
- `AGENTS.md`: not updated; no hard repository constraint changed.
