# PGT-A Snakemake Pipeline

This pipeline is for PGT-A analysis with `fastp + bwa + WisecondorX`.

## 0. Recent updates

- Refactored NPZ parsing to follow official WisecondorX convert format (`binsize + sample + quality`).
- Added reference prefilter robustness: unusable NPZ samples are logged and excluded before model QC.
- Improved `WisecondorX convert` throughput with parallel execution (`--threads`) and per-sample logs.
- Added stricter pipeline reproducibility logs (git branch/commit + tool versions).
- Added repository cleanup guidance and ignore rules for temporary scratch files.

Supported targets:
- `mapping`: fastp + bwa + samtools sort/index
- `reference_qc`: prefilter + bin/pca tuning (when tuning is enabled)
- `reference`: build final reference model
- `cnv_qc`: CNV input QC (when CNV is enabled)
- `cnv`: WisecondorX predict (when CNV is enabled)

## 1. Initialize config and update reference groups

Entry script: `scripts/init_wes_project.py`

Available actions:
- `init_config`: initialize project structure and `config.yaml`
- `build_reference_groups`: update `build_reference.groups` from Excel

### 1.1 Init config

```bash
python scripts/init_wes_project.py \
  --action init_config \
  --project /path/to/project \
  --fq_dir /path/to/fastq \
  --output_config /path/to/project/config.yaml \
  --template_snakefile /home/jiucheng/pipelines/PGT_A/Snakefile
```

Optional: init config and fill XX/XY groups from Excel in one step.

```bash
python scripts/init_wes_project.py \
  --action init_config \
  --project /path/to/project \
  --fq_dir /path/to/fastq \
  --output_config /path/to/project/config.yaml \
  --build_reference_mode excel \
  --sample_info_xlsx /home/jiucheng/pipelines/PGT_A/sample_info/CNV-seq.xlsx \
  --batch2_fq_dir /data/project/CNV/PGT-A/rawdata/lib_test/2026-02-05 \
  --batch3_fq_dir /data/project/CNV/PGT-A/rawdata/lib_test/2026-03-03
```

### 1.2 Update only build_reference groups

```bash
python scripts/init_wes_project.py \
  --action build_reference_groups \
  --project /path/to/project \
  --output_config /path/to/project/config.yaml \
  --sample_info_xlsx /home/jiucheng/pipelines/PGT_A/sample_info/CNV-seq.xlsx \
  --batch2_fq_dir /data/project/CNV/PGT-A/rawdata/lib_test/2026-02-05 \
  --batch3_fq_dir /data/project/CNV/PGT-A/rawdata/lib_test/2026-03-03
```

Notes:
- `--project` is required for `build_reference_groups`.
- If `--output_config` does not exist, a template config is created automatically.
- `samples` contains all batch2 + batch3 samples.
- `build_reference.groups` is generated from Excel rules.

## 2. Default key parameters (from init_wes_project)

- `core.wisecondorx.binsize: 100000`
- `core.wisecondorx.tuning.bin_sizes: [100000, 200000, 300000, 500000, 750000, 1000000]`
- `core.wisecondorx.tuning.pca.min_explained_variance: 0.0`
- `core.wisecondorx.reference_prefilter.binsize: 100000`
- `core.wisecondorx.cnv.convert_binsize: 100000`
- `core.wisecondorx.cnv.zscore: 5`
- `core.wisecondorx.cnv.alpha: 0.001`

## 3. Reference output layout

Default root: `core.wisecondorx.reference_model_root: reference`

- `reference/XX/prefilter/`
- `reference/XX/tuning/`
- `reference/XX/result/`
- `reference/XY/prefilter/`
- `reference/XY/tuning/`
- `reference/XY/result/`

## 4. Run the pipeline

Use Snakemake in your server env:
`/biosoftware/miniconda/envs/snakemake_env/bin/snakemake`

### 4.1 Mapping

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake \
  -s /home/jiucheng/pipelines/PGT_A/Snakefile \
  --configfile /path/to/project/config.yaml \
  --cores 32 -p mapping
```

### 4.2 Build reference (XX and XY in one run)

Run reference QC stage first:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake \
  -s /home/jiucheng/pipelines/PGT_A/Snakefile \
  --configfile /path/to/project/config.yaml \
  --cores 32 -p reference_qc
```

Then build final reference:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake \
  -s /home/jiucheng/pipelines/PGT_A/Snakefile \
  --configfile /path/to/project/config.yaml \
  --cores 32 --rerun-incomplete -p reference
```

Or run all together:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake \
  -s /home/jiucheng/pipelines/PGT_A/Snakefile \
  --configfile /path/to/project/config.yaml \
  --cores 32 --rerun-incomplete -p mapping reference_qc reference
```

### 4.3 CNV

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake \
  -s /home/jiucheng/pipelines/PGT_A/Snakefile \
  --configfile /path/to/project/config.yaml \
  --cores 32 -p mapping cnv_qc cnv
```

## 5. Main outputs

Reference:
- `reference/{SEX}/prefilter/reference_sample_qc.tsv`
- `reference/{SEX}/prefilter/reference_sample_qc.svg`
- `reference/{SEX}/prefilter/reference_inlier_samples.txt`
- `reference/{SEX}/prefilter/prefilter_summary.yaml`
- `reference/{SEX}/tuning/bin_pca_grid.tsv`
- `reference/{SEX}/tuning/best_params.yaml`
- `reference/{SEX}/tuning/best_bin_pca_elbow.svg`
- `reference/{SEX}/tuning/reference_qc_metrics.svg`
- `reference/{SEX}/tuning/reference_inlier_samples.txt`
- `reference/{SEX}/result/*.npz`

CNV:
- `wisecondorx/cnv/qc/{sample}.qc.tsv`
- `wisecondorx/cnv/qc/{sample}.qc.svg`
- `wisecondorx/cnv/predict/{sample}*`

## 6. Logs and audit trail

Logs are written to:
- `logs/fastp/`
- `logs/bwa/`
- `logs/metadata/`
- `logs/wisecondorx/`
- `logs/cnv/`

Python scripts now use the `logging` module for pipeline step tracing
(start/end, command execution, sample filtering, QC summaries, and failure reasons).
For key rules, script logging is persisted via `--log` in addition to rule-level shell logs.

Audit metadata in logs includes:
- `git_branch`
- `git_commit`
- tool versions (`fastp`, `bwa`, `samtools`, `WisecondorX`, `python`)

## 7. Performance notes

- `scripts/reference_prefilter_qc.py` and `scripts/tune_wisecondorx_bin_pca.py` now run `WisecondorX convert` in parallel, controlled by `--threads`.
- Each sample convert job writes its own log under `.../converted/convert_logs/`.
- Main pipeline log keeps a compact progress summary for convert start/done events.
- Some repeated list/set traversals were reduced in sample filtering and outlier set construction.
- NPZ parsing now follows official WisecondorX layout: `binsize + sample(dict chr1..24) + quality`, avoiding expensive generic object scans.

## 8. Troubleshooting

### 8.1 IndentationError in Snakefile/pipeline.smk

```bash
cd /home/jiucheng/pipelines/PGT_A
sed -i 's/\r$//' Snakefile rules/pipeline.smk
nl -ba rules/pipeline.smk | sed -n '70,90p;240,260p'
```

### 8.2 WisecondorX convert argument error

`WisecondorX convert` in this pipeline no longer uses `--cpus` because some versions do not support it.

### 8.3 `No usable numeric arrays found in .../*.npz`

`reference_prefilter_qc.py` now auto-drops unusable NPZ samples and records them in log/summary.
NPZ parsing follows official WisecondorX convert output layout:
- `binsize`
- `sample` (dict with chromosome keys `1..24`)
- `quality`

If remaining usable samples are fewer than `--min-reference-samples`, the rule still fails by design.
In that case:
- check sample-level convert logs under `reference/{SEX}/prefilter/converted/convert_logs/`
- remove problematic BAM/NPZ from reference groups or regenerate those BAM files
- rerun `reference_qc`

### 8.4 Repository temporary files

Temporary files (for ad-hoc debugging/source lookups) are not part of the pipeline and should not remain in repo root.
Ignored by default:
- `tmp_*`
- `*.tmp`

## 9. Workflow structure (modularized)

The workflow is now organized into modular rule files:

- `Snakefile`: config loading, target assembly, top-level entry rules, and includes
- `rules/common_preprocess.smk`: shared preprocessing (`collect_run_metadata`, `fastp_bwa`)
- `rules/reference_workflow.smk`: reference-only stages (prefilter, tuning, reference build)
- `rules/predict_workflow.smk`: predict-only stages (`convert`, CNV QC, CNV predict)

This keeps common preprocessing and downstream analysis modes separated while preserving existing rule names, outputs, and DAG behavior.

## 10. Config entry points

In addition to `config.yaml`, two mode-specific config entry files are available:

- `config_reference.yaml`
  - `base_config: "config.yaml"`
  - `pipeline.targets: ["mapping", "metadata", "reference_qc", "reference"]`
- `config_predict.yaml`
  - `base_config: "config.yaml"`
  - `pipeline.targets: ["mapping", "metadata", "cnv_qc", "cnv"]`

Run examples:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake \
  -s /home/jiucheng/pipelines/PGT_A/Snakefile \
  --configfile /path/to/project/config_reference.yaml \
  --cores 32 -p
```

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake \
  -s /home/jiucheng/pipelines/PGT_A/Snakefile \
  --configfile /path/to/project/config_predict.yaml \
  --cores 32 -p
```

## 11. PGT-A algorithm notes

This section explains algorithm concepts only. Runtime behavior, thresholds, QC logic, CNV logic, mosaic logic, sex logic, and output schema are still defined by code/config.

### 11.1 End-to-end stages

1. Raw data preprocessing
- FASTQ input
- FASTQ QC (`fastp`)
- mapping (`bwa`) + BAM sort/index (`samtools`)

2. Bin-level signal construction
- convert aligned reads into fixed-size genomic bin signals

3. Reference building
- estimate normal-background distributions on reference samples
- optional prefilter and tuning before final reference model generation

4. Predict
- process test sample with reference-consistent settings
- quantify deviation from reference
- call CNV results

### 11.2 Notation

- `X_ij`: raw count/signal of sample `i` at bin `j`
- `X~_ij`: normalized/corrected bin signal
- `mu_j`, `sigma_j`: reference mean/std at bin `j`
- `Z_j`: per-bin z-score for test sample

### 11.3 Typical formulas

Library-size normalization:

`X_norm_ij = X_ij / sum_j X_ij`

z-score:

`Z_j = (X~_tj - mu_j) / sigma_j`

Optional floor-stabilized form:

`Z_j = (X~_tj - mu_j) / max(sigma_j, sigma_min)`

If log-ratio is used by a specific implementation path:

`R_j = log2((X~_tj + eps) / (mu_j + eps))`

### 11.4 PCA / denoising intuition

Reference matrix (sample x bin) can be centered and projected to principal components to capture systematic technical variance (for example GC/batch/amplification effects). Tuning determines suitable bin size / PCA dimension combinations. In this pipeline, this behavior is implemented by tuning scripts and WisecondorX-driven workflows.

### 11.5 Segmentation note

CNV inference commonly performs segmentation over ordered bin signals. Exact segmentation internals are tool-implementation specific; this pipeline delegates core prediction behavior to `WisecondorX predict`.

## 12. Implementation mapping

| Module | Rule / Script | Purpose | Main input | Main output |
|---|---|---|---|---|
| metadata | `collect_run_metadata` / `scripts/collect_run_metadata.py` | collect run/tool metadata | project + tool paths | `logs/run_metadata.tsv` |
| preprocess | `fastp_bwa` | fastq QC + mapping | sample R1/R2 | cleaned FASTQ, sorted BAM/BAI |
| reference prefilter | `reference_prefilter_*` / `scripts/reference_prefilter_qc.py` | remove unusable/outlier reference samples | reference BAMs | prefilter QC/report/inlier list |
| reference tuning | `tune_wisecondorx_reference_qc*` / `scripts/tune_wisecondorx_bin_pca.py` | tune bin/PCA and sample QC | prefilter inliers + BAMs | tuning summary/best params/inliers |
| reference build | `build_wisecondorx_reference_*` / `scripts/build_reference_from_tuning.py` | generate final reference | tuning outputs | reference `.npz` |
| predict convert | `wisecondorx_convert_for_cnv` | convert BAM to CNV NPZ | sample BAM | CNV NPZ |
| predict QC | `wisecondorx_qc_for_predict` / `scripts/cnv_qc.py` | per-sample CNV input QC | CNV NPZ | QC TSV/SVG/pass marker |
| predict | `wisecondorx_predict_cnv` | CNV calling | sample NPZ + reference NPZ | predict outputs + done marker |
