# Handoff: Branch B V2 Low-Resolution Materialization

timestamp: 2026-06-21_2339
status: in_progress_reference_build_started
active_reference_id: h_r0_shadow_ref_20260619

## Context Entry

Use this handoff after `docs/CURRENT_CONTEXT_INDEX.md`.

Required read order:

1. `docs/CURRENT_CONTEXT_INDEX.md`
2. `docs/handoff/2026-06-21_2339_branch_b_v2_lowres_materialization_handoff.md`
3. `AGENTS.md`
4. `skills/conversation_handoff/SKILL.md`
5. `skills/pgta_reference_modeling_analysis/SKILL.md`
6. Current remote result files and active configs

## Completed In This Loop

- Confirmed PR #26 merged the lowres/ref-MAD interface into `main`.
- Created branch `codex/lowres-ref-materialization` from current `main`.
- Added 2Mb/3Mb shadow reference configs.
- Added 2Mb/3Mb low-res Branch A predict configs for Y, H, G, and 2026-06-15.
- Added lowres-enabled primary benchmark/report configs for Y, H, G, and
  2026-06-15.
- Added local G1-G8 gap2m configs that had existed in the remote mirror but not
  in local main.
- Extended `rules/predict_layout.smk` so lowres evidence can use either
  explicit `reference_npz` paths or `reference_npz_glob`.
- Updated workflow-contract tests for the new glob option.
- Synced the staged code/config snapshot to the remote mirror for validation.

## Remote Validation

Remote unit test command:

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

Reference dry-runs:

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
Both exited 0. Each DAG contained 19 jobs and did not schedule BWA/mapping.
```

Predict dry-run status:

```text
Y/H 2Mb and 3Mb predict dry-runs currently fail with MissingInput for the
low-res reference common_best_binsize.txt files. This is expected before the
2Mb/3Mb reference builds finish.
```

## Long Task Running

Sequential 2Mb then 3Mb reference build driver:

```text
PID: 71117
PID file: /data/project/CNV/PGT-A/refactor_validation_20260419/logs/driver/lowres_ref_2mb_3mb_20260622.pid
Log: /data/project/CNV/PGT-A/refactor_validation_20260419/logs/driver/lowres_ref_2mb_3mb_20260622.log
```

Command:

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile /data/project/CNV/PGT-A/refactor_validation_20260419/config_reference_h_r0_shadow_2mb_20260622.yaml --cores 60 reference &&
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile /data/project/CNV/PGT-A/refactor_validation_20260419/config_reference_h_r0_shadow_3mb_20260622.yaml --cores 60 reference
```

Initial check:

```text
Driver PID running.
2Mb child Snakemake running.
2Mb DAG entered metadata, annotation, and baseline QC jobs.
No mapping/BWA was triggered.
No runtime.db was present, so no historical runtime estimate is available.
```

## Next Steps

1. Poll reference build by checking target outputs and log, not PID alone.
2. After 2Mb/3Mb refs complete, run low-res predict dry-runs for all four
   cohorts and confirm no mapping/BWA dependency.
3. Materialize low-res Branch A events for all four cohorts at 2Mb/3Mb.
4. Dry-run and materialize lowres-enabled primary benchmark/report configs.
5. Verify truth preservation:
   - Y1-Y8: 10/10, FN=0.
   - H1-H6: 10/10, FN=0, H6 chr21 visible.
   - G1-G8: 10/10, FN=0, G2 truth visible.
   - truth filtered count=0.
6. Summarize 2Mb/3Mb same-direction support, especially for `>=4Mb` truth and
   3-4Mb boundary-sensitive events.

## Constraints To Preserve

- Do not modify Branch A thresholds in this gate.
- Do not rebuild or promote the 1Mb primary shadow ref.
- Do not treat low-res absence as a standalone filter.
- Do not let low-res evidence generate B-only events.
- Keep Branch S as `review_reportable_with_limitations`.
- Keep 2026-06-15 as burden/context only unless locked truth labels are added.
