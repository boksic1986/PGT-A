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
- `logs/wisecondorx/`
- `logs/cnv/`

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
