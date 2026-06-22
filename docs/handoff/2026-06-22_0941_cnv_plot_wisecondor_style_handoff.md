# Handoff: CNV Plot Wisecondor-Style Refinement

Date: 2026-06-22 09:41

## 1. Goal

Optimize the per-sample CNV SVG plot style without changing Branch A, Branch B
V2, Branch S, report classification, or reference outputs.

User-requested plot behavior:

- do not connect z-trend lines across the genome or across chromosomes
- draw trend lines only over final duplication/deletion report intervals
- draw those trend lines as horizontal red lines
- use Wisecondor-like yellow for duplication bins
- keep deletion bins blue

## 2. Confirmed Project State

- Active reference remains `h_r0_shadow_ref_20260619`.
- Active 0615 config remains
  `config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`.
- 0615 still has no locked truth labels, so no TP/FP/FN conclusion is allowed.
- This loop did not modify candidate classification, Branch S, report
  visibility, or any threshold.

## 3. Completed Changes

Modified:

- `pgta/predict/branch_b/plot.py`
- `tests/unit/test_branch_b_plot.py`
- `docs/reports/report_main_convergence_cnv_plot_2026-06-22.md`
- `CURRENT_STATE.md`
- `PLANS.md`
- `docs/CURRENT_CONTEXT_INDEX.md`

Plot implementation changes:

- Removed the prior genome/chromosome smooth polyline rendering.
- Added event-local `report-z-trend` lines.
- Each `report-z-trend` line is horizontal and spans only the final autosomal
  report event interval.
- The y value is the event-bin median `calibrated_z`.
- Duplication bins are `#facc15`.
- Deletion bins remain `#2563eb`.
- Trend lines are `#dc2626`.
- Legend now says `report z trend`, not `smooth z trend`.

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
  tests/unit/test_cnv_report.py \
  tests/unit/test_branch_ab_phase12_workflow_contract.py -q
```

Result:

```text
32 passed in 0.89s
```

Dry-run:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml \
  --cores 1 -n --forcerun cnv_branch_ab_plot cnv_report
```

Result:

```text
exit 0; planned cnv_branch_ab_plot for 5 samples plus cnv_report refresh;
no mapping or reference rebuild jobs were requested
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

Local synced plot path:

`D:\Pipeline\PGT-A\reports\0615_cnv_plots\*.final_cnv.svg`

## 5. Current Conclusion

The CNV plot style now matches the requested behavior for the 0615 materialized
SVGs. This is a plot/rendering-only change. It does not reduce report burden or
change any candidate classification.

## 6. Next Step

Use the synced 0615 SVGs for manual sample-by-sample interpretation, starting
from sample 56 if continuing the previous review. Any future report-burden
reduction still needs candidate-level evidence review and locked-truth
ablation.

## 7. Core File Sync

- `CURRENT_STATE.md`: updated with the refined CNV plot behavior.
- `PLANS.md`: updated to make the visual review criterion the red horizontal
  report-event trend.
- `docs/CURRENT_CONTEXT_INDEX.md`: updated to point to this handoff as the
  latest plot-refinement context.
- `REPO_MAP.md`: not updated; no structure or entrypoint changed.
- `AGENTS.md`: not updated; no repository hard constraint changed.
