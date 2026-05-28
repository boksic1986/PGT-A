# Build Ref V2 Validation

## Scope

First-round Build Ref V2 implementation only:

- WisecondorX-first contract report.
- Static reference mask classification.
- Mask-only WisecondorX-compatible NPZ preprocessing.
- `WisecondorX predict --blacklist` passthrough.
- Branch B HMM sidecar default.
- Reference build final-newref mask-only preprocessing hook.
- Predict pre-`predict` mask-only NPZ route.

GC/RC correction remains offline-only design work and was not promoted into `newref`.

## Remote Commands

All validation commands below were executed on `fengxian` under:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419
```

### WisecondorX Help And Source Inspect

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && mkdir -p docs/reports && /biosoftware/miniconda/envs/wise_env/bin/WisecondorX newref --help > docs/reports/wisecondorx_newref_help_2026-05-27.txt 2>&1 && /biosoftware/miniconda/envs/wise_env/bin/WisecondorX predict --help > docs/reports/wisecondorx_predict_help_2026-05-27.txt 2>&1'
```

Result: success.

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && /biosoftware/miniconda/envs/wise_env/bin/python /tmp/wisecondorx_source_inspect_20260527.py > docs/reports/wisecondorx_source_inspect_2026-05-27.txt 2>&1'
```

Result: success.

Key findings:

- `newref` has no `--blacklist`.
- `predict` has `--blacklist`.
- `newref_tools.get_mask` retains bins where normalized cohort count sum is `> 0`.
- `newref_tools.get_reference` excludes same-chromosome bins during reference-bin selection.
- `train_pca(masked_data)` is called without exposed CLI PCA component control.

## Unit Tests

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 /biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest tests/unit/test_reference_preprocess.py tests/unit/test_wisecondorx_predict_blacklist.py tests/unit/test_branch_b_calling.py -q'
```

Result:

```text
19 passed in 0.79s
```

Validated:

- chrM hard mask classification.
- PAR is annotation-only, not global hard mask.
- gap/centromere/telomere hard mask classification.
- GC/RC fitting mask uses clean autosomal pass bins only.
- mask-only NPZ hard bins become zero and count-like contract is preserved.
- `WisecondorX predict` rule includes blacklist passthrough and explicit skip logging.
- `WisecondorX predict` consumes `CNV_PREDICT_INPUT_NPZ`; default `mask_only` adds `cnv_mask_npz_for_predict`.
- Branch B HMM sidecar/disabled modes do not emit standalone segment candidates.
- `legacy_candidate` remains available only for comparison behavior.

## Snakemake Dry-Run

### Predict

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && /biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile /data/project/CNV/PGT-A/test_project/config_predict.yaml --cores 1 -n cnv cnv_report'
```

Result: success dry-run.

Summary:

```text
total jobs: 88
included rules include:
- build_reference_bin_annotations
- build_reference_masks
- wisecondorx_convert_for_cnv
- cnv_mask_npz_for_predict
- wisecondorx_gender_for_predict
- wisecondorx_predict_cnv
- cnv_correction_branch_b
- cnv_calling_branch_b
- cnv_report
```

### Reference

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && /biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile /tmp/build_ref_v2_config_reference_with_baseline.yaml --cores 1 -n reference'
```

Result: success dry-run.

```text
total jobs: 55
included rules include:
- build_reference_bin_annotations
- build_reference_masks
- tune_wisecondorx_reference_qc_by_cohort
- build_wisecondorx_reference_from_tuning_by_cohort_sex
- build_wisecondorx_gender_reference_from_tuning_by_cohort
```

Interpretation:

- `build_reference_bin_annotations` no longer depends on `TUNING_BEST`.
- `build_reference_masks` is upstream of final WisecondorX `newref` rules.
- final reference rules now receive `combined_mask.tsv` and can use `mask_only` preprocessing.

## Mask-Only NPZ Compatibility

Single-sample transform check:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 /biosoftware/miniconda/envs/snakemake_env/bin/python -m pgta.reference.preprocess mask_npz --input-npz /data/project/CNV/PGT-A/test_project/wisecondorx/tuning/pass_warn/bin_100000/converted/3.npz --combined-mask /data/project/CNV/PGT-A/test_project/reference/assets/mask/combined_mask.tsv --output-npz /tmp/build_ref_v2_mask_npz_check/3.mask_only.npz --output-summary /tmp/build_ref_v2_mask_npz_check/3.mask_summary.json'
```

Result:

```text
preprocess_strategy=mask_only
mask_policy=hard_mask_to_zero_soft_keep
hard_masked_bin_count=7518
count_like_contract_validated=true
```

Small WisecondorX `newref` compatibility check:

- Masked all existing `/data/project/CNV/PGT-A/test_project/wisecondorx/tuning/pass_warn/bin_100000/converted/*.npz` files into `/tmp/build_ref_v2_newref_check_many/masked/`.
- Ran `/biosoftware/miniconda/envs/wise_env/bin/WisecondorX newref ... --binsize 100000 --cpus 1`.
- Output `/tmp/build_ref_v2_newref_check_many/ref_mask_only.npz` was generated, size `101M`.
- `/biosoftware/miniconda/envs/snakemake_env/bin/python` successfully loaded the output with `numpy.load`; observed standard WisecondorX reference keys including `bins_per_chr`, `binsize`, `distances`, `has_female`, `has_male`, `mask`, `masked_bins_per_chr`.

Note: a 3-sample trial failed inside WisecondorX gender-model training (`IndexError`) because the validation set was too small; the larger existing cohort generated a loadable reference.

## WisecondorX-First Status

Preserved:

- `WisecondorX convert` remains the input generator.
- `WisecondorX newref` remains the reference builder.
- `WisecondorX predict` remains the primary CNV caller.
- WisecondorX CBS/segment output remains primary evidence.

Not introduced:

- No custom HMM primary caller.
- No custom segmentation replacing CBS.
- No logR/z-score/centered corrected NPZ sent to `newref`.
- No hard-mask median neutralization in default `newref` path.

## Current Limitations

- Tuning still uses current project chr1-22 matrix. It remains project QC and is not claimed to reproduce WisecondorX internal reference-bin selection.
- GC/RC correction is not promoted into reference build. It remains a Phase 3 offline evaluation task.
- The reference dry-run used `/tmp/build_ref_v2_config_reference_with_baseline.yaml`, derived from test_project config with `baseline_qc` added to requested targets, so the original `config_reference.yaml` alone still is not the complete reference validation entrypoint.

## Promotion State

Current recommendation:

- Proceed with full mask-only reference validation on the intended retained cohort.
- Do not promote GC/RC corrected NPZ to `newref` yet.
- Keep Branch B HMM sidecar as default.

## 2026-05-28 Mask-Only Reference Predict Validation

Validation target:

```text
Y1-Y8 positive samples
reference root: /data/project/CNV/PGT-A/refactor_validation_20260419/results_build_ref_v2_mask_only
predict config: /data/project/CNV/PGT-A/refactor_validation_20260419/config_predict_build_ref_v2_mask_only.yaml
```

The validation config points `WisecondorX predict` to the completed mask-only reference outputs:

```text
XX: /data/project/CNV/PGT-A/refactor_validation_20260419/results_build_ref_v2_mask_only/reference/XX/result/ref_xx_best.npz
XY: /data/project/CNV/PGT-A/refactor_validation_20260419/results_build_ref_v2_mask_only/reference/XY/result/ref_xy_best.npz
gender: /data/project/CNV/PGT-A/refactor_validation_20260419/results_build_ref_v2_mask_only/reference/gender/result/ref_gender_best.npz
common binsize: /data/project/CNV/PGT-A/refactor_validation_20260419/results_build_ref_v2_mask_only/reference/gender/common_best_binsize.txt
blacklist BED: /data/project/CNV/PGT-A/refactor_validation_20260419/results_build_ref_v2_mask_only/reference/assets/mask/hard_mask.bed
```

The blacklist path is now a BED export derived from `hard_mask.tsv`, not the TSV itself. This keeps `WisecondorX predict --blacklist` contract-compatible.

Remote unit test:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 /biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest tests/unit/test_reference_preprocess.py tests/unit/test_wisecondorx_predict_blacklist.py tests/unit/test_branch_b_calling.py -q'
```

Result:

```text
20 passed in 0.81s
```

Remote dry-run:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && /biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile /data/project/CNV/PGT-A/refactor_validation_20260419/config_predict_build_ref_v2_mask_only.yaml --cores 64 --rerun-triggers mtime -n cnv cnv_eval cnv_benchmark cnv_report'
```

Result:

```text
success dry-run
total jobs: 81
no fastp/bwa/mapping jobs selected
selected rules include:
- export_reference_hard_mask_bed
- wisecondorx_convert_for_cnv
- cnv_mask_npz_for_predict
- wisecondorx_gender_for_predict
- wisecondorx_predict_cnv
- cnv_correction_branch_b
- cnv_calling_branch_b
- cnv_calibration_branch_b
- cnv_artifact_rules_branch_b
- cnv_evaluation
- cnv_benchmark_framework
- cnv_report_summary
```

The `--rerun-triggers mtime` flag is intentional for this validation run. Without it, Snakemake provenance changes would try to rerun reference asset construction even though the mask-only reference has already been built. This run freezes the built reference assets and only derives the compatible `hard_mask.bed` from the existing hard mask TSV.

Background run:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && /data/project/CNV/PGT-A/refactor_validation_20260419/logs/build_ref_v2/predict_ref_v2_mask_only_20260528_101126.cmd.sh'
```

Runtime metadata:

```text
PID: 100936
Snakemake child PID: 100938
log: /data/project/CNV/PGT-A/refactor_validation_20260419/logs/build_ref_v2/predict_ref_v2_mask_only_20260528_101126.log
pid file: /data/project/CNV/PGT-A/refactor_validation_20260419/logs/build_ref_v2/predict_ref_v2_mask_only_20260528_101126.pid
command script: /data/project/CNV/PGT-A/refactor_validation_20260419/logs/build_ref_v2/predict_ref_v2_mask_only_20260528_101126.cmd.sh
```

Initial health check:

```text
status: running
active stage: WisecondorX convert for Y1-Y8
active WisecondorX convert processes observed: 8
CPU: about 8 active cores during convert
memory: low
disk: /data 56% used
monitor/runtime.db: not found, so no historical runtime estimate is available
```

Concurrency interpretation:

```text
The workflow was launched with --cores 64.
The first active stage is limited by the eight Y samples and by WisecondorX convert behavior.
No BWA/remapping jobs were selected because existing sorted BAM/BAI files under the mask-only output root are reused.
```

### Predict Validation Completion

The first background run failed before CNV prediction because the validation config pointed `mosaic_fraction_truth_tsv` at:

```text
/data/project/CNV/PGT-A/test_project/sample_info/mosaic_fraction_validation_y7_y8.tsv
```

That file does not exist. The correct file is:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419/sample_info/mosaic_fraction_validation_y7_y8.tsv
```

The second run then failed in `WisecondorX predict --blacklist` because the first `hard_mask.bed` export was BED4 (`chrom/start/end/reason`), while this WisecondorX executable imports blacklist files with strict three-column parsing:

```text
chr_name, s, e = line.strip().split("\t")
```

Root-cause fix:

```text
pgta/reference/assets.py write_hard_mask_bed now emits strict BED3 only.
tests/unit/test_reference_preprocess.py now asserts BED3 output.
```

Remote verification:

```text
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 /biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest tests/unit/test_reference_preprocess.py tests/unit/test_wisecondorx_predict_blacklist.py -q
8 passed in 0.55s

hard_mask.bed column validation:
bad_lines=0
```

Final successful run:

```text
PID: 113818
log: /data/project/CNV/PGT-A/refactor_validation_20260419/logs/build_ref_v2/predict_ref_v2_mask_only_bed3_20260528_102633.log
command script: /data/project/CNV/PGT-A/refactor_validation_20260419/logs/build_ref_v2/predict_ref_v2_mask_only_bed3_20260528_102633.cmd.sh
result: 39 of 39 steps (100%) done
```

Final targets:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419/results_build_ref_v2_mask_only/wisecondorx/cnv/evaluation/summary.json
/data/project/CNV/PGT-A/refactor_validation_20260419/results_build_ref_v2_mask_only/wisecondorx/cnv/benchmarking/summary.json
/data/project/CNV/PGT-A/refactor_validation_20260419/results_build_ref_v2_mask_only/wisecondorx/cnv/benchmarking/mosaic_truth_validation.json
/data/project/CNV/PGT-A/refactor_validation_20260419/results_build_ref_v2_mask_only/wisecondorx/cnv/report/cnv_summary.tsv
/data/project/CNV/PGT-A/refactor_validation_20260419/results_build_ref_v2_mask_only/wisecondorx/cnv/report/cnv_summary.html
```

Key metrics:

```text
sample_count: 8
truth_event_count: 10
A-branch detected truth events: 10 / 10
A-branch detection rate: 1.0
Branch B detected truth events: 3 / 10
Branch B detection rate: 0.3
Branch B event count: 25
Branch B kept event count: 22
Branch B PASS event count: 0
Branch B REVIEW event count: 22
Branch B artifact event count: 3
```

Truth-hit breakdown:

```text
Detected by Branch B:
- Y1 chr21 loss
- Y2 chr21 gain
- Y4 chr13 gain

Detected by Branch A but not retained by Branch B:
- Y2 chr14 gain
- Y3 chrX loss
- Y5 chr15 loss
- Y6 chr7 loss
- Y6 chr7 gain
- Y7 chr8 loss
- Y8 chr13 gain
```

Interpretation:

```text
The mask-only reference did not reduce Branch A truth sensitivity in this Y1-Y8 validation set.
Branch B remains too conservative or misaligned with its intended role as refinement of Branch A evidence.
The next iteration should focus on Branch B retaining/classifying A-supported true events instead of acting as an independent gate that drops them.
```

## 2026-05-28 Branch A Candidate Refinement Iteration

Implementation scope:

```text
Branch A:
  WisecondorX aberrations BED -> normalized A candidate table
  same chrom/state adjacent intervals are merged
  A candidate table keeps WisecondorX ratio/zscore/rank/support metadata

Branch B:
  consumes Branch A candidate table as its primary candidate input
  refines A candidates with corrected-bin, calibration, artifact, fraction, and sidecar HMM evidence
  HMM remains sidecar by default and does not emit standalone final candidates
  candidate annotation / artifact / fraction labeling remains in Branch B
```

Remote tests:

```text
cd /data/project/CNV/PGT-A/refactor_validation_20260419
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
  /biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_b_artifact_rules.py \
  tests/unit/test_predict_truth.py \
  tests/unit/test_benchmark.py \
  tests/unit/test_branch_a_candidates.py \
  tests/unit/test_branch_b_calling.py -q

result:
34 passed in 0.88s
```

Remote dry-run:

```text
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake \
  -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile \
  --configfile /data/project/CNV/PGT-A/refactor_validation_20260419/config_predict_build_ref_v2_mask_only.yaml \
  --cores 64 --rerun-triggers mtime --rerun-incomplete -n \
  cnv cnv_eval cnv_benchmark cnv_report

result:
dry-run succeeded; 39 jobs selected.
No BWA, mapping, WisecondorX convert/predict, or reference rebuild jobs were selected.
```

Remote execution:

```text
command:
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake \
  -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile \
  --configfile /data/project/CNV/PGT-A/refactor_validation_20260419/config_predict_build_ref_v2_mask_only.yaml \
  --cores 64 --rerun-triggers mtime --rerun-incomplete --latency-wait 60 \
  -R cnv_artifact_rules_branch_b cnv_evaluation cnv_benchmark_framework cnv_report_summary \
  cnv cnv_eval cnv_benchmark cnv_report

PID: 12657
log: /data/project/CNV/PGT-A/refactor_validation_20260419/logs/build_ref_v2/predict_ref_v2_a_refine_b_iter2_force_20260528_130548.log
result: 17 of 17 steps (100%) done
```

Key metrics after A-refined Branch B:

```text
truth_event_count: 10
A-branch detected truth events: 10 / 10
A-branch detection rate: 1.0
Branch B detected truth events: 10 / 10
Branch B detection rate: 1.0
truth_match_count: 10 / 10
truth_recall: 1.0

event_count: 86
kept_event_count: 44
PASS event count: 0
REVIEW event count: 44
artifact event count: 42
```

Truth-hit breakdown:

```text
Y1 chr21 loss: detected by A and A-refined B
Y2 chr14 gain: detected by A and A-refined B
Y2 chr21 gain: detected by A and A-refined B
Y3 chrX loss: detected by A and A-refined B
Y4 chr13 gain: detected by A and A-refined B
Y5 chr15 loss: detected by A and A-refined B
Y6 chr7 loss: detected by A and A-refined B
Y6 chr7 gain: detected by A and A-refined B
Y7 chr8 loss: detected by A and A-refined B
Y8 chr13 gain: detected by A and A-refined B
```

Interpretation:

```text
This iteration fixes the main FN failure mode: Branch B no longer drops A-supported truth events by behaving as an independent primary gate.
The remaining limitation is classification, not recall: all retained events are still REVIEW and no event is promoted to PASS.
The next iteration should tune Branch B PASS/REVIEW/artifact stratification and reduce review burden while preserving the 10/10 recall.
```
