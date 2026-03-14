# Changelog

## 2026-03-14

- Refactored workflow organization into modular rule files:
  - `rules/common_preprocess.smk`
  - `rules/reference_workflow.smk`
  - `rules/predict_workflow.smk`
- Slimmed top-level `Snakefile` to focus on config loading, target assembly, and includes.
- Added mode-specific config entry files:
  - `config_reference.yaml`
  - `config_predict.yaml`
- Updated reference group generation logic in `scripts/init_wes_project.py`:
  - Fixed Excel columns to `batch`, `sample_ids`, `sample_type`, `sex`
  - Batch2 includes only `A-H` and `XY`
  - Batch3 ID extraction uses the first numeric ID from names like `JZ26040817-1-1_combined`
  - Output group IDs normalized for batch2 (`A-H`) and batch3 (`1-26`)
- Added BAM-level QC utility:
  - `scripts/bam_uniformity_qc.py`
  - Outputs `qc_metrics.tsv`, `target_bin_profile.tsv`, `reference_summary.tsv`, and `target_vs_ref_profile.png`
- Added tmux remote workflow memo:
  - `tmux_tut.md`
