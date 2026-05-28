# Build Ref V2 Current Contract

## WisecondorX Executable

- newref executable path: `/biosoftware/miniconda/envs/wise_env/bin/WisecondorX`
- predict executable path: `/biosoftware/miniconda/envs/wise_env/bin/WisecondorX`
- remote project path: `/data/project/CNV/PGT-A/refactor_validation_20260419`

Remote commands recorded:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && /biosoftware/miniconda/envs/wise_env/bin/WisecondorX newref --help > docs/reports/wisecondorx_newref_help_2026-05-27.txt 2>&1'
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && /biosoftware/miniconda/envs/wise_env/bin/WisecondorX predict --help > docs/reports/wisecondorx_predict_help_2026-05-27.txt 2>&1'
```

Observed CLI contract:

- `newref` supports `--nipt`, `--yfrac`, `--plotyfrac`, `--refsize`, `--binsize`, `--cpus`.
- `newref` does not expose `--blacklist`.
- `predict` supports `--blacklist`.
- No `newref` PCA component CLI parameter is exposed.

Remote evidence files:

- `/data/project/CNV/PGT-A/refactor_validation_20260419/docs/reports/wisecondorx_newref_help_2026-05-27.txt`
- `/data/project/CNV/PGT-A/refactor_validation_20260419/docs/reports/wisecondorx_predict_help_2026-05-27.txt`
- `/data/project/CNV/PGT-A/refactor_validation_20260419/docs/reports/wisecondorx_source_inspect_2026-05-27.txt`

## Current Reference Build

Current primary path remains WisecondorX-first:

```text
WisecondorX convert
-> project-level reference QC/tuning
-> WisecondorX newref
```

Current `pgta/reference/tune.py` reads WisecondorX NPZ files as `sample` chromosome count vectors, builds a project-level chr1-22 matrix, and computes diagnostic PCA/QC metrics. This project PCA is not a verified WisecondorX `newref` parameter.

`pgta/reference/build.py` calls `WisecondorX newref` directly with NPZ inputs, output reference path, `--binsize`, `--cpus`, and optional `--yfrac`.

## Current Predict

Current primary predict path remains:

```text
WisecondorX predict
-> WisecondorX BED/plot/CBS-derived output
-> project postprocess/report
```

Build Ref V2 first-round changes make `WisecondorX predict` pass `--blacklist <path>` when configured and present as a non-empty file. If the path is empty, missing, or empty, the predict log records why no blacklist was passed.

## Source Inspect Findings

Remote inspect command recorded WisecondorX source from:

- `/biosoftware/miniconda/envs/wise_env/lib/python3.8/site-packages/wisecondorx/newref_tools.py`
- `/biosoftware/miniconda/envs/wise_env/lib/python3.8/site-packages/wisecondorx/predict_tools.py`

Confirmed from `newref_tools.get_mask`:

- WisecondorX builds per-bin mask from normalized count matrix.
- Bins are retained where `sum_per_bin > 0`.
- Therefore Build Ref V2 mask-only preprocessing uses hard-mask-to-zero for `newref`, not per-chromosome median neutralization.

Confirmed from `newref_tools.get_reference`:

- For a target chromosome, reference-bin selection concatenates bins before and after that chromosome and excludes bins on the same chromosome.
- Source-inspect status: same-chromosome exclusion verified by source inspect.

Confirmed from `newref_control.py` grep:

- `train_pca(masked_data)` is called without a CLI-controlled PCA component argument.
- `train_pca` default is `pcacomp=5` in inspected source.
- Project `best_pca_components` / `project_qc_pca_components` must not be reported as a WisecondorX runtime parameter unless a future patched executable proves otherwise.

## Current Risks

- `newref` lacks `--blacklist`; reference masking must be implemented through count-like NPZ preprocessing or a validated WisecondorX patch.
- Hard-mask median neutralization may pollute WisecondorX reference-bin selection and is not default.
- Corrected NPZ files must remain positive, finite, count-scale vectors; logR/z-score/centered signals are prohibited as `newref` inputs.
- Project-level tuning is diagnostic and not the WisecondorX internal reference-bin selection algorithm.
- Branch B HMM must stay sidecar by default and must not become a second primary caller.
- Gender should remain based on raw NPZ/BAM depth until validated otherwise.
- PAR and XY homology remain annotation/review labels by default, not global hard masks.
