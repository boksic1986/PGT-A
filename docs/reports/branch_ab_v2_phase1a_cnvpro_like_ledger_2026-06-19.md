# Branch A/B V2 Phase 1A CNVpro-like Ledger

Date: 2026-06-19

## Scope

This change implements Phase 1A from `branch_ab_v2_rnd_constraints_2026-06-18.md`.

The implementation extends the existing Branch B candidate evidence ledger with CNVpro/CNVseq-inspired shadow evidence fields. It does not change:

- WisecondorX predict;
- Branch A candidate inclusion;
- legacy Branch B final report decisions;
- Branch S or sex calling;
- `cnv_report` output schema.

## Added Shadow Fields

The ledger now records:

- GC/RC background availability;
- dynamic-reference/refmap status derived from existing ref-bin evidence;
- matched-negative placeholder status;
- copy-number and mosaic proxy values;
- arm/PAR/centromere context;
- whole-chromosome and large-segment tier;
- waviness/sample-noise context;
- a shadow-only CNVpro-like evidence status.

Missing background remains `UNKNOWN_BACKGROUND` or `NaN`. Missing evidence is not converted into clean support.

## Design Notes

The copy-number and mosaic fields are proxies for ablation review, not production thresholds. The current proxy uses Branch A ratio when no event-level copy-number field exists:

```text
copy_number_estimate = 2 * (1 + a_ratio)
```

The large-segment tier is also shadow-only. A candidate is not labeled `whole_chromosome` solely because it covers most of the bins present in a partial input table; the event must be at least 20 Mb and cover most available chromosome bins.

## Remote Validation

All validation ran on `fengxian` under:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419
```

Commands:

```bash
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
tests/unit/test_branch_b_evidence_ledger.py -q
```

Result:

```text
4 passed in 0.58s
```

```bash
PYTHONPATH=/data/project/CNV/PGT-A/refactor_validation_20260419 \
/biosoftware/miniconda/envs/snakemake_env/bin/python -m pytest \
tests/unit/test_branch_ab_phase12_workflow_contract.py \
tests/unit/test_branch_b_evidence_ledger.py \
tests/unit/test_branch_b_matched_negative.py \
tests/unit/test_branch_b_v2_classifier.py \
tests/unit/test_negative_bank.py \
tests/unit/test_branch_b_calling.py \
tests/unit/test_branch_b_artifact_rules.py \
tests/unit/test_branch_b_correction.py \
tests/unit/test_branch_s_shadow.py -q
```

Result:

```text
83 passed in 1.07s
```

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake \
-s /data/project/CNV/PGT-A/refactor_validation_20260419/Snakefile \
--configfile /data/project/CNV/PGT-A/refactor_validation_20260419/config_predict_build_ref_v2_mask_only.yaml \
--cores 1 -n cnv cnv_report
```

Result:

```text
Nothing to be done (all requested files are present and up to date).
```

## Remaining Work

- Materialize refreshed Phase 1A ledgers for Y1-Y8, H1-H16, and 2026-06-15 if the current outputs need the new columns immediately.
- Add per-field missingness summaries after materialization.
- Use Phase 3 matched-negative evidence when N0/N1/N2 labels are finalized; do not treat current `UNKNOWN_BACKGROUND` as clean support.
- Run locked-data ablation before promoting any CNVpro-like evidence into report decisions.
