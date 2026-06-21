# Branch B V2 Low-Resolution Materialization

status: in_progress_reference_build_started
date: 2026-06-22
active_reference_id: h_r0_shadow_ref_20260619
branch_a_input: WisecondorX/CBS with explicit merge_gap_bp=2,000,000 overlay

## Scope

This gate materializes the already implemented 2Mb/3Mb low-resolution evidence
interface. The 1Mb Branch A chain remains the primary WisecondorX/CBS discovery
layer. The 2Mb/3Mb outputs are auxiliary evidence only: they can strengthen
confidence when same-direction support is present, but absence of support must
not suppress a candidate by itself.

## Configs Added

Reference configs:

- `config_reference_h_r0_shadow_2mb_20260622.yaml`
- `config_reference_h_r0_shadow_3mb_20260622.yaml`

Low-resolution Branch A predict configs:

- `config_predict_y_h_r0_shadow_2mb_branch_a_gap2m_20260622.yaml`
- `config_predict_y_h_r0_shadow_3mb_branch_a_gap2m_20260622.yaml`
- `config_predict_h_20260608_h_r0_shadow_2mb_branch_a_gap2m_20260622.yaml`
- `config_predict_h_20260608_h_r0_shadow_3mb_branch_a_gap2m_20260622.yaml`
- `config_predict_g1_g8_h_r0_shadow_2mb_branch_a_gap2m_20260622.yaml`
- `config_predict_g1_g8_h_r0_shadow_3mb_branch_a_gap2m_20260622.yaml`
- `config_predict_20260615_h_r0_shadow_2mb_branch_a_gap2m_20260622.yaml`
- `config_predict_20260615_h_r0_shadow_3mb_branch_a_gap2m_20260622.yaml`

Lowres-enabled primary benchmark/report configs:

- `config_predict_y_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`
- `config_predict_h_20260608_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`
- `config_predict_g1_g8_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`
- `config_predict_20260615_h_r0_shadow_branch_b_v2_gap2m_lowres_evidence_20260622.yaml`

## Workflow Contract Update

`rules/predict_layout.smk` now accepts
`core.wisecondorx.cnv.lowres_evidence.reference_npz_glob` in addition to the
existing explicit `reference_npz` list. The glob is resolved during workflow
layout construction and deduplicated before use. This keeps the primary ref
NPZ source explicit while avoiding repeated 38-path lists in each config.

The active glob points to the current primary ref inlier mask-only NPZ set:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419/results_h_r0_shadow_ref_20260619/wisecondorx/tuning/pass_warn_h_r0/bin_750000/mask_only_for_newref/*.npz
```

Remote inspection showed this glob resolves to 38 files at the time of this
handoff.

## Remote Validation So Far

Unit tests on `ssh fengxian`:

```bash
cd /data/project/CNV/PGT-A/refactor_validation_20260419
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
  tests/unit/test_branch_ab_phase12_workflow_contract.py \
  tests/unit/test_lowres_evidence.py \
  tests/unit/test_ref_stability.py \
  tests/unit/test_branch_b_v2_classifier.py -q
```

Result:

```text
59 passed in 1.25s
```

Reference dry-runs on `ssh fengxian`:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile config_reference_h_r0_shadow_2mb_20260622.yaml \
  --cores 1 -n reference

/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s Snakefile \
  --configfile config_reference_h_r0_shadow_3mb_20260622.yaml \
  --cores 1 -n reference
```

Result:

```text
Both dry-runs exited 0. Each DAG contains 19 reference/QC/tuning/newref jobs.
No BWA/mapping jobs were scheduled.
```

Low-resolution predict dry-run before reference completion was intentionally
blocked by missing low-res reference output:

```text
Missing input:
results_h_r0_shadow_ref_20260619_2mb/reference/gender/common_best_binsize.txt
results_h_r0_shadow_ref_20260619_3mb/reference/gender/common_best_binsize.txt
```

This is the expected dependency order. Predict/benchmark/report materialization
must wait for 2Mb/3Mb reference completion.

## Long Task Started

No historical runtime database was available:

```text
monitor/runtime.db: not present
```

The 2Mb and 3Mb reference builds were started as one sequential background
driver to cap concurrent usage at 60 Snakemake cores:

```text
PID: 71117
PID file: /data/project/CNV/PGT-A/refactor_validation_20260419/logs/driver/lowres_ref_2mb_3mb_20260622.pid
Log: /data/project/CNV/PGT-A/refactor_validation_20260419/logs/driver/lowres_ref_2mb_3mb_20260622.log
Command:
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile /data/project/CNV/PGT-A/refactor_validation_20260419/config_reference_h_r0_shadow_2mb_20260622.yaml --cores 60 reference &&
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile /data/project/CNV/PGT-A/refactor_validation_20260419/config_reference_h_r0_shadow_3mb_20260622.yaml --cores 60 reference
```

Initial status check:

```text
Driver PID 71117 was running.
Child Snakemake PID 71123 was running the 2Mb reference DAG.
The 2Mb run had entered collect_run_metadata, build_reference_bin_annotations,
and baseline_bam_uniformity_qc.
No mapping/BWA step was triggered.
```

Expected reference targets:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419/results_h_r0_shadow_ref_20260619_2mb/reference/gender/common_best_binsize.txt
/data/project/CNV/PGT-A/refactor_validation_20260419/results_h_r0_shadow_ref_20260619_2mb/reference/XX/ref.npz
/data/project/CNV/PGT-A/refactor_validation_20260419/results_h_r0_shadow_ref_20260619_2mb/reference/XY/ref.npz
/data/project/CNV/PGT-A/refactor_validation_20260419/results_h_r0_shadow_ref_20260619_3mb/reference/gender/common_best_binsize.txt
/data/project/CNV/PGT-A/refactor_validation_20260419/results_h_r0_shadow_ref_20260619_3mb/reference/XX/ref.npz
/data/project/CNV/PGT-A/refactor_validation_20260419/results_h_r0_shadow_ref_20260619_3mb/reference/XY/ref.npz
```

## Pending

After the 2Mb/3Mb reference build completes:

1. Dry-run and materialize low-res Branch A predict for Y1-Y8, H1-H16,
   G1-G8, and 2026-06-15 at 2Mb/3Mb.
2. Dry-run and materialize the lowres-enabled primary benchmark/report configs.
3. Verify:
   - Y1-Y8 truth preserved 10/10, FN=0.
   - H1-H6 truth preserved 10/10, FN=0, H6 chr21 visible.
   - G1-G8 truth preserved 10/10, FN=0, G2 truth visible.
   - truth filtered count=0.
   - 2026-06-15 remains burden/context only.
4. Output `>=4Mb` truth same-direction support rates for 2Mb/3Mb.
5. Treat 3-4Mb events primarily by 2Mb evidence; 3Mb absence is
   `boundary_diluted_context`, not negative evidence.

## Current Boundary

This report does not claim that low-res evidence has passed performance
validation. It records that configs, workflow contract support, unit tests, and
reference dry-runs are in place, and that the long-running reference build has
started.
