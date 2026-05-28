# Build Ref V2 WisecondorX-First Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Rework reference build so masking and GC/RC correction are applied before reference modeling, while keeping WisecondorX as the primary CNV method and using Branch B only for refinement and review.

**Architecture:** Add a reference preprocessing layer that creates a single set of mask/annotation assets, applies mask-aware GC/RC correction to reference and predict NPZ inputs, and runs WisecondorX `newref`/`predict` on consistently transformed inputs. Keep HMM as an optional sidecar consistency signal, not a main caller. Re-tune binsize and PCA after preprocessing instead of assuming 1 Mb and PCA20 are optimal.

**Tech Stack:** Snakemake, Python, WisecondorX, pandas/numpy, pytest/unittest, remote validation on `ssh fengxian` with `/biosoftware/miniconda/envs/snakemake_env/bin/python` and `/biosoftware/miniconda/envs/wise_env/bin/WisecondorX`.

---

## Scope And Non-Goals

This plan changes the build reference architecture and predict input contract. It does not change final result schema unless a later task explicitly proposes and validates a schema extension.

In scope:
- Move unified blacklist/mask assets into the reference build path.
- Add mask-aware GC/RC correction before reference modeling.
- Ensure predict samples undergo the same mask/correction transform as reference samples.
- Re-run binsize/PCA tuning on transformed data.
- Keep WisecondorX as the primary caller.
- Add `--blacklist` to WisecondorX predict.
- Convert Branch B HMM from main-caller behavior to optional sidecar evidence.
- Validate with current Y/G truth and low-fraction sensitivity checks.

Out of scope:
- Copy-neutral UPD calling.
- Replacing WisecondorX CBS segmentation with a custom HMM caller.
- Changing mosaic/sex/result schema before validation.
- Treating local test output as workflow validation; all real tests run on `fengxian`.

## Current Findings To Preserve

- Remote WisecondorX confirms:
  - `WisecondorX predict` supports `--blacklist`.
  - `WisecondorX newref` does not support `--blacklist`.
  - `WisecondorX newref` exposes `--refsize`, `--binsize`, `--yfrac`, and `--cpus`.
- Current `pgta/reference/build.py` calls `WisecondorX newref` without blacklist support because the executable has no such option.
- Current `pgta/reference/tune.py` builds a custom chr1-22 matrix for QC/tuning. It does not explicitly implement target-bin same-chromosome exclusion for correlation-style reference-bin selection.
- Current GC/mappability LOESS exists in `pgta/predict/branch_b/correction.py`, but it is downstream postprocessing, not a reference-build preprocessing contract.
- Current predict workflow uses `WisecondorX predict`, but blacklist assets are mainly used by Branch B annotation/artifact rules. The `--blacklist` passthrough must be made explicit.

## Target Data Flow

```text
reference BAMs
-> WisecondorX convert raw NPZ
-> project bin annotation and combined mask assets
-> mask-aware GC/RC correction on clean autosomal bins
-> corrected/masked NPZ
-> tuning on corrected/masked matrix
-> WisecondorX newref on corrected/masked NPZ

predict BAM
-> WisecondorX convert raw NPZ
-> same mask-aware GC/RC correction
-> corrected/masked NPZ
-> WisecondorX predict with corrected ref and --blacklist
-> Branch B refinement/report
```

## File Structure

Modify:
- `pgta/reference/assets.py`: ensure combined blacklist/mask BEDs explicitly encode chrM exclusion and hard/soft mask classes.
- `pgta/reference/tune.py`: make tuning consume transformed NPZ/matrix inputs and record preprocessing strategy in outputs.
- `pgta/reference/build.py`: build references from transformed NPZ paths and record preprocessing metadata.
- `pgta/predict/branch_b/correction.py`: extract reusable GC/RC correction functions into a shared module.
- `pgta/predict/branch_b/calling.py`: keep HMM output optional and avoid treating it as the primary caller in default workflow.
- `rules/reference_layout.smk`: add paths for transformed NPZs, correction summaries, mask assets, and refbuild-v2 outputs.
- `rules/reference_workflow.smk`: insert mask/correction rules before tuning and reference build.
- `rules/predict_layout.smk`: add predict corrected NPZ path and blacklist-predict configuration.
- `rules/predict_workflow.smk`: apply the same correction to predict NPZ and pass `--blacklist` to WisecondorX predict.
- `README.md`: document WisecondorX-first refbuild-v2 and Branch B sidecar role.
- `docs/handoff/*.md` or a new handoff at the end: record decisions and validation state.

Create:
- `pgta/reference/preprocess.py`: shared mask-aware matrix/NPZ transform and GC/RC correction entrypoint.
- `tests/unit/test_reference_preprocess.py`: unit coverage for mask/correction contracts.
- `tests/unit/test_wisecondorx_predict_blacklist.py`: rule/command construction checks for `--blacklist`.
- `docs/reports/build_ref_v2_validation_YYYY-MM-DD.md`: final validation summary after remote experiments.

## Task 1: Audit And Freeze Current Reference Contract

**Files:**
- Create: `docs/reports/build_ref_v2_current_contract_2026-05-27.md`

- [ ] **Step 1: Write current contract report**

Document current behavior:

```markdown
# Build Ref Current Contract

## WisecondorX executable
- `newref --help`: no `--blacklist` option.
- `predict --help`: supports `--blacklist`.

## Current reference build
- Converts BAM to NPZ for candidate binsizes.
- Builds a custom chr1-22 matrix in `pgta/reference/tune.py`.
- Runs WisecondorX `newref` on selected inlier NPZ files.
- No explicit pre-newref blacklist/mask transform.

## Current predict
- Runs WisecondorX `predict`.
- Must be checked/changed to pass `--blacklist`.

## Current risk
- Reference samples and predict samples may not pass through the same mask/correction transform.
- Branch B correction currently compensates after predict instead of defining the comparison space before modeling.
```

- [ ] **Step 2: Verify remote executable help**

Run:

```bash
ssh fengxian '/biosoftware/miniconda/envs/wise_env/bin/WisecondorX newref --help | head -120'
ssh fengxian '/biosoftware/miniconda/envs/wise_env/bin/WisecondorX predict --help | head -160'
```

Expected:
- `newref` shows no `--blacklist`.
- `predict` shows `--blacklist`.

- [ ] **Step 3: Record command results**

Append executable, command, and result to the contract report.

## Task 2: Define Unified Mask Asset Contract

**Files:**
- Modify: `pgta/reference/assets.py`
- Modify: `rules/reference_assets.smk`
- Test: `tests/unit/test_reference_preprocess.py`

- [ ] **Step 1: Write failing tests for mask classes**

Add tests that construct small fake annotation rows and assert:

```python
def test_chrM_is_hard_excluded_from_reference_mask():
    row = build_mask_row(chrom="chrM", start=1, end=1000)
    assert row["mask_label"] == "hard"
    assert "chrM" in row["mask_reason"]

def test_par_is_label_not_global_hard_mask():
    row = build_mask_row(chrom="chrX", start=60001, end=2699520, par_overlap_fraction=1.0)
    assert row["is_PAR"] == 1
    assert row["mask_label"] != "hard"

def test_gap_blacklist_is_hard_mask():
    row = build_mask_row(chrom="chr1", start=1, end=1000, gap_centromere_telomere_overlap_fraction=1.0)
    assert row["mask_label"] == "hard"
```

- [ ] **Step 2: Run failing test remotely**

Run:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 /biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest tests/unit/test_reference_preprocess.py -q'
```

Expected: FAIL before implementation.

- [ ] **Step 3: Implement mask classification helper**

Add a helper in `pgta/reference/assets.py`:

```python
def classify_reference_mask(row):
    chrom = str(row.get("chrom", ""))
    if chrom in {"chrM", "MT", "M"}:
        return "hard", "chrM_excluded"
    if float(row.get("blacklist_overlap_fraction", 0.0) or 0.0) > 0.0:
        return "hard", "blacklist"
    if float(row.get("gap_centromere_telomere_overlap_fraction", 0.0) or 0.0) >= 0.50:
        return "hard", "gap_centromere_telomere"
    if float(row.get("low_mappability_overlap_fraction", 0.0) or 0.0) >= 0.25:
        return "soft", "low_mappability"
    if float(row.get("segmental_duplication_overlap_fraction", 0.0) or 0.0) >= 0.25:
        return "soft", "segmental_duplication"
    if float(row.get("repeat_rich_overlap_fraction", 0.0) or 0.0) >= 0.25:
        return "soft", "repeat_rich"
    return "pass", ""
```

Then use it when writing combined mask assets. Keep PAR and sex homology as annotation labels, not default hard masks.

- [ ] **Step 4: Run unit tests remotely**

Run:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 /biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest tests/unit/test_reference_preprocess.py -q'
```

Expected: PASS.

## Task 3: Extract GC/RC Correction Into Reference Preprocess Module

**Files:**
- Create: `pgta/reference/preprocess.py`
- Modify: `pgta/predict/branch_b/correction.py`
- Test: `tests/unit/test_reference_preprocess.py`

- [ ] **Step 1: Write tests for correction fit scope**

Add tests that assert:

```python
def test_gc_rc_fit_uses_only_clean_autosomal_bins():
    bins = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chrX", "chr2"],
            "normalized_signal": [10.0, 11.0, 30.0, 12.0],
            "gc_fraction": [0.40, 0.45, 0.50, 0.60],
            "mappability_score": [1.0, 1.0, 1.0, 1.0],
            "mask_label": ["pass", "soft", "pass", "pass"],
            "is_autosome": [1, 1, 0, 1],
        }
    )
    fit_mask = build_gc_rc_fit_mask(bins, include_mask_labels={"pass"})
    assert fit_mask.tolist() == [True, False, False, True]
```

Also test that correction returns both:
- `signal_for_reference`
- correction summary fields with `fit_bin_count`, `model`, `loess_frac`.

- [ ] **Step 2: Implement shared functions**

Create `pgta/reference/preprocess.py` with:

```python
def build_gc_rc_fit_mask(bins_df, include_mask_labels):
    labels = bins_df["mask_label"].fillna("pass").astype(str)
    autosome = bins_df["is_autosome"].fillna(0).astype(int).eq(1)
    finite = (
        pd.to_numeric(bins_df["normalized_signal"], errors="coerce").notna()
        & pd.to_numeric(bins_df["gc_fraction"], errors="coerce").notna()
    )
    return autosome & labels.isin(include_mask_labels) & finite
```

Move reusable local regression helpers from Branch B correction into this module:
- `tricube`
- `bisquare_robust_weights`
- `local_linear_predict`
- `fit_loess_model`
- `apply_gc_rc_correction`

Keep the implementation deterministic and single-sample scoped.

- [ ] **Step 3: Keep Branch B compatibility**

Update `pgta/predict/branch_b/correction.py` to import shared helpers instead of duplicating them. Its outputs must remain unchanged.

- [ ] **Step 4: Run tests remotely**

Run:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 /biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest tests/unit/test_reference_preprocess.py tests/unit/test_branch_b_correction.py -q'
```

Expected: PASS.

## Task 4: Generate Corrected/Masked NPZ For Reference Tuning

**Files:**
- Modify: `pgta/reference/preprocess.py`
- Modify: `rules/reference_layout.smk`
- Modify: `rules/reference_workflow.smk`
- Test: `tests/unit/test_reference_preprocess.py`

- [ ] **Step 1: Define transformed NPZ contract**

Use WisecondorX-compatible NPZ layout:

```text
binsize: original binsize
sample: dict of chromosome numeric vectors
quality: original quality if available
pgta_preprocess: metadata object
```

Metadata must include:

```python
{
    "preprocess_version": "build_ref_v2",
    "mask_policy": "hard_mask_to_zero_soft_keep_weighted",
    "correction_model": "loess_gc_rc_mappability",
    "fit_bin_count": int,
    "hard_masked_bin_count": int,
}
```

- [ ] **Step 2: Implement corrected NPZ writer**

Add a CLI action in `pgta/reference/preprocess.py`:

```bash
python -m pgta.reference.preprocess correct_npz \
  --input-npz sample.raw.npz \
  --annotations reference/assets/bins/analysis_bin_annotations.tsv \
  --combined-mask reference/assets/masks/analysis_combined_mask.tsv \
  --output-npz sample.corrected.npz \
  --output-summary sample.correction_summary.json
```

Hard-masked bins should be neutralized consistently. For first implementation, set hard-masked bins to the per-chromosome median of clean corrected signal rather than zero, to avoid artificial CNV breakpoints.

- [ ] **Step 3: Add Snakemake paths**

In `rules/reference_layout.smk`, add:

```python
REF_CORRECTED_NPZ = project_path("wisecondorx", "tuning", "{cohort}", "bin_{binsize}", "corrected", "{sample}.npz")
REF_CORRECTION_SUMMARY = project_path("wisecondorx", "tuning", "{cohort}", "bin_{binsize}", "correction", "{sample}.summary.json")
```

- [ ] **Step 4: Insert correction after convert**

In `rules/reference_workflow.smk`, after WisecondorX convert for each candidate binsize, run correction for every selected reference sample before `run_tune_wisecondorx` loads NPZ arrays.

- [ ] **Step 5: Run unit tests remotely**

Run:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 /biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest tests/unit/test_reference_preprocess.py -q'
```

Expected: PASS.

## Task 5: Retune Binsize/PCA On Transformed Matrix

**Files:**
- Modify: `pgta/reference/tune.py`
- Modify: `rules/reference_workflow.smk`
- Test: `tests/unit/test_reference_preprocess.py`
- Test: existing `tests/unit/test_reference_cohort_resequencing.py`

- [ ] **Step 1: Add preprocessing metadata columns**

Extend `bin_pca_grid.tsv` and `best_params.yaml` with:

```text
preprocess_strategy
mask_policy
correction_model
fit_bin_median
hard_masked_bin_median
```

- [ ] **Step 2: Extend candidate grid**

Use config-controlled binsize list. For the validation run, include:

```yaml
core:
  wisecondorx:
    tuning:
      bin_sizes: [100000, 250000, 500000, 1000000]
      pca:
        min_components: 2
        max_components: 30
        min_explained_variance: 0.80
```

Do not hard-code PCA20 as default.

- [ ] **Step 3: Add same-chrom exclusion audit metric**

Do not claim to reproduce WisecondorX internal reference-bin selection unless verified. Add a tuning diagnostic that reports whether project-level correlation/CV metrics are:

```text
global_chr1_22_matrix
leave_one_chromosome_out_available=false
```

Then add a later optional metric:

```text
leave_one_chromosome_out_reconstruction_mse
```

The first implementation must at least make this limitation explicit in `best_params.yaml`.

- [ ] **Step 4: Run reference tuning dry-run remotely**

Run:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && /biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile <refbuild_v2_config.yaml> --cores 1 -n -- reference'
```

Expected: DAG parses and corrected NPZ/tuning outputs are dependencies before reference build.

## Task 6: Build WisecondorX References From Corrected NPZ

**Files:**
- Modify: `pgta/reference/build.py`
- Modify: `rules/reference_workflow.smk`
- Test: add command-level tests if repository already has rule construction tests.

- [ ] **Step 1: Preserve WisecondorX newref as reference builder**

`build_reference_from_npz_paths()` should still call:

```bash
/biosoftware/miniconda/envs/wise_env/bin/WisecondorX newref <corrected_npz...> <ref.npz> --binsize <best_binsize> --cpus <threads>
```

Do not attempt to pass `--blacklist` to `newref`.

- [ ] **Step 2: Add metadata sidecar**

For every ref output, write:

```text
reference/.../result/ref_*.metadata.json
```

Include:

```python
{
    "source_npz_type": "pgta_gc_rc_corrected_npz",
    "best_binsize": best_binsize,
    "best_pca_components": best_pca,
    "preprocess_strategy": "mask_gc_rc_v2",
    "wisecondorx_newref_supports_blacklist": False,
}
```

- [ ] **Step 3: Rework reference groups**

For v2 validation:
- Build autosome/mixed reference from all retained XX+XY samples.
- Preserve XX/XY sex-specific references for sex chromosomes if current WisecondorX output format does not support chromosome-partitioned refs.
- Build gender ref from mixed samples, retaining manual `--yfrac` logic.

If chromosome-partitioned refs are not feasible in WisecondorX, document that current implementation remains whole-genome references but all autosome tuning metrics are mixed-sex.

## Task 7: Apply Same Transform To Predict Inputs And Pass Blacklist

**Files:**
- Modify: `rules/predict_layout.smk`
- Modify: `rules/predict_workflow.smk`
- Modify: `pgta/reference/preprocess.py`
- Test: `tests/unit/test_wisecondorx_predict_blacklist.py`

- [ ] **Step 1: Write test for `--blacklist` predict command**

Add a unit or static test asserting that when `CNV_POSTPROCESS_BLACKLIST_BED` is non-empty, the `wisecondorx_predict_cnv` rule command includes:

```bash
--blacklist <path>
```

- [ ] **Step 2: Add predict corrected NPZ**

Add:

```python
CNV_CORRECTED_NPZ = str(Path(CNV_DIR) / "npz_corrected" / "{sample}.npz")
CNV_CORRECTION_SUMMARY = str(Path(CNV_DIR) / "npz_corrected" / "{sample}.correction_summary.json")
```

- [ ] **Step 3: Insert correction before predict**

New predict flow:

```text
wisecondorx_convert_for_cnv -> cnv_correct_npz_for_predict -> wisecondorx_gender_for_predict / wisecondorx_predict_cnv
```

Use the same reference assets and correction settings as build ref.

- [ ] **Step 4: Keep gender robust**

Evaluate whether gender should use raw NPZ or corrected NPZ. First validation should run both and compare sex calls. Do not switch gender to corrected NPZ unless calls match expected labels.

## Task 8: Demote Branch B HMM To Sidecar Evidence

**Files:**
- Modify: `pgta/predict/branch_b/calling.py`
- Modify: `rules/predict_workflow.smk`
- Test: `tests/unit/test_branch_b_calling.py`
- Documentation: `README.md`

- [ ] **Step 1: Add config flag**

Add config option:

```yaml
core:
  wisecondorx:
    cnv:
      postprocess:
        branch_b:
          hmm_role: sidecar
```

Accepted values:

```text
sidecar
disabled
legacy_candidate
```

Default: `sidecar`.

- [ ] **Step 2: Make HMM non-primary by default**

In sidecar mode:
- Keep `hmm_state` in bins.
- Compute consistency metrics.
- Do not create final candidate events solely from HMM segment calls unless they overlap WisecondorX candidates or are explicitly marked as rescue review.

- [ ] **Step 3: Preserve old behavior behind flag**

If `hmm_role=legacy_candidate`, preserve current behavior for comparison runs only.

- [ ] **Step 4: Unit test**

Assert:

```python
def test_sidecar_hmm_does_not_emit_standalone_candidate():
    events = build_candidate_events_with_hmm_role(hmm_role="sidecar", wisecondorx_candidates=[])
    assert events == []
```

## Task 9: Remote Validation Matrix

**Files:**
- Create: `docs/reports/build_ref_v2_validation_2026-05-27.md`

- [ ] **Step 1: Prepare configs**

Create remote configs under:

```text
/data/project/CNV/PGT-A/refbuild_v2_validation_20260527/
```

Variants:

```text
baseline_current
mask_only
mask_gc_rc
mask_gc_rc_pca80
mask_gc_rc_pca_grid
```

- [ ] **Step 2: Run dry-runs**

Run each dry-run on remote:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && /biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile <config> --cores 1 -n -- reference cnv cnv_report'
```

Expected: DAG parses.

- [ ] **Step 3: Launch long tasks with PID/log handoff**

For each long Snakemake run, follow project rule:
- return PID
- log path
- target paths
- command
- runtime estimate from `monitor/runtime.db` if available
- do not block indefinitely

- [ ] **Step 4: Compare metrics**

Collect:

```text
reference inlier count
best_binsize
best_pca_components
cumulative explained variance
CV reconstruction MSE
baseline bin CV
chr-level shift
Y1-Y8 truth recall
G1-G8 truth recall if truth table is active
Branch A merged TP/FP/FN
Branch B TP/FP/FN
10/20/30/50% admixture sensitivity
sex call concordance
```

- [ ] **Step 5: Promotion criteria**

Promote v2 default only if:
- truth recall is not lower than current PCA20 baseline
- low-fraction sensitivity does not regress
- FP or review burden improves
- sex calls remain correct
- correction does not flatten known broad gains/losses

## Task 10: Documentation And Handoff

**Files:**
- Modify: `README.md`
- Create: `docs/reports/build_ref_v2_validation_2026-05-27.md`
- Create: `docs/handoff/YYYY-MM-DD_HHMM_build_ref_v2_handoff.md`

- [ ] **Step 1: Update README**

Document:
- WisecondorX is the primary CNV method.
- WisecondorX CBS remains primary segmentation.
- Project mask/GC-RC preprocessing defines the comparison space.
- Branch B HMM is sidecar by default.
- Autosome modeling is mixed-sex; sex chromosome interpretation remains sex-aware.

- [ ] **Step 2: Write validation report**

Include final table:

```text
variant  ref_samples  best_binsize  pca  explained_var  Y_recall  G_recall  FP  FN  low_fraction_recall  decision
```

- [ ] **Step 3: Write handoff**

Use `skills/conversation_handoff/SKILL.md` format. Include:
- completed changes
- remote commands
- paths
- validation state
- unresolved risks

## Validation Commands

All validation commands must run on remote `fengxian`.

Unit tests:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 /biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest tests/unit/test_reference_preprocess.py tests/unit/test_branch_b_correction.py tests/unit/test_wisecondorx_predict_blacklist.py -q'
```

Reference dry-run:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && /biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile <refbuild_v2_config.yaml> --cores 1 -n -- reference'
```

Predict dry-run:

```bash
ssh fengxian 'cd /data/project/CNV/PGT-A/refactor_validation_20260419 && /biosoftware/miniconda/envs/snakemake_env/bin/snakemake -s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile --configfile <predict_v2_config.yaml> --cores 1 -n -- cnv cnv_report'
```

Full validation should be launched as a long task with PID/log handoff, not by blocking the conversation.

## Risks And Review Questions

- WisecondorX `newref` does not accept blacklist; transformed NPZ must remain compatible with its expected layout.
- Neutralizing hard-masked bins may affect WisecondorX internal reference-bin selection; validation must compare mask-only versus mask+GC-RC.
- GC/RC correction can flatten real broad events if fit bins are not clean enough; fit must use clean autosomal bins only.
- Gender calling may need raw NPZ rather than corrected NPZ; validate both before switching.
- WisecondorX internal CBS/reference-bin selection should remain the primary method. Project-level HMM must not silently become the default caller.
- Same-chrom exclusion in WisecondorX internals is not directly controlled by this code; project-level tuning should not overclaim equivalence until verified.

## GPT Pro Review Checklist

Ask the reviewer to focus on:
- Whether transformed NPZ before `newref` is a valid way to apply mask/GC-RC given WisecondorX lacks `newref --blacklist`.
- Whether hard-mask neutralization by per-chrom clean median is safer than zeroing bins.
- Whether GC/RC correction should be sample-specific, cohort-level, or hybrid.
- Whether autosome mixed-sex reference plus sex-aware interpretation is correct for PGT-A/CNV-seq.
- Whether HMM sidecar role is appropriately constrained.
- Whether validation metrics are sufficient to prove no recall regression.
