# Branch A/B V2 Phase 1B Asset Decision Matrix

Date: 2026-06-19

## Scope

Phase 1B is an asset decision step, not a caller change.

It consolidates current hg19 mask, mappability, PAR/maskpar, CNVseq unique-bin, and WisecondorX ref-bin evidence into a workflow decision matrix. It does not change WisecondorX predict, Branch A candidate inclusion, legacy Branch B final events, Branch S, or `cnv_report`.

Current default remains WisecondorX-first:

- WisecondorX predict/CBS is the primary CNV evidence.
- Branch A is the high-sensitivity candidate discovery layer.
- Branch B V2 evidence remains candidate-level shadow evidence until locked ablation proves no recall loss.
- Branch S remains sex-chromosome/SCA shadow evidence and does not replace sex calling or final reporting.

## Evidence Sources

Repository files:

- `docs/reports/branch_ab_v2_rnd_constraints_2026-06-18.md`
- `docs/handoff/2026-06-11_2358_hg19_repeat_mask_probe_handoff.md`
- `pgta/reference/assets.py`
- `rules/reference_assets.smk`
- `pgta/predict/branch_b/evidence_ledger.py`

Remote evidence was refreshed on `fengxian` from:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419/results_build_ref_v2_mask_only/reference/assets
```

Remote command:

```bash
cd /data/project/CNV/PGT-A/refactor_validation_20260419 &&
/biosoftware/miniconda/envs/snakemake_env/bin/python -
```

## Current Asset Inventory

Current reference FASTA recorded in the asset summary:

```text
/data/Database/index/hg19/hg19.fa
```

Current bin levels:

| bin level | bin size | total bins | PAR bins | sex-homology bins | segdup bins | low-map bins | blacklist bins | near-telomere bins | near-centromere bins |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| atomic | 50 kb | 61,927 | 61 | 61 | 8,427 | 547 | 4,881 | 4,652 | 6,268 |
| analysis | 100 kb | 30,970 | 31 | 31 | 6,245 | 414 | 2,573 | 2,348 | 3,148 |
| QC | 1 Mb | 3,113 | 5 | 5 | 2,220 | 234 | 437 | 277 | 340 |

Current combined mask summary:

| label | total bins |
|---|---:|
| pass | 81,406 |
| hard | 7,892 |
| soft | 6,710 |
| dynamic | 2 |

Mask labels by level:

| bin level | hard | soft | dynamic | pass |
|---|---:|---:|---:|---:|
| atomic | 4,881 | 4,271 | 0 | 52,775 |
| analysis | 2,573 | 2,197 | 0 | 26,200 |
| QC | 438 | 242 | 2 | 2,431 |

Current hard-mask BED:

```text
/data/project/CNV/PGT-A/refactor_validation_20260419/results_build_ref_v2_mask_only/reference/assets/mask/hard_mask.bed
```

It contains 4,881 atomic intervals and is exported as strict BED3 for WisecondorX predict blacklist input.

## Current Mask Semantics

`pgta/reference/assets.py::classify_reference_mask()` currently uses:

| asset or condition | current label | current meaning |
|---|---|---|
| chrM | hard | excluded from reference mask logic. |
| blacklist overlap > 0 | hard | exported to hard-mask BED through atomic hard bins. |
| gap / centromere / telomere overlap >= 0.50 | hard | unavailable or boundary-prone region. |
| WisecondorX ref bins < 50 | hard | insufficient reference-bin support. |
| WisecondorX ref bins < 100 plus high repeat/low-map/ambiguous evidence | hard | unstable region with insufficient reference-bin support. |
| proximal telomere/centromere plus ref bins < 150 | hard | boundary-prone and low reference support. |
| WisecondorX ref bins < 150 only | dynamic | not a hard artifact by itself. |
| proximal telomere/centromere plus high dynamic noise | hard | boundary-prone and empirically unstable. |
| low mappability >= 0.25 | soft | artifact evidence only. |
| segmental duplication >= 0.25 | soft | artifact evidence only. |
| repeat-rich >= 0.25 | soft | artifact evidence only. |
| PAR / sex homology | annotation only | not global hard mask. |

At the analysis 100 kb level, current high-repeat and mapping-related annotations are common:

| annotation | bins with overlap > 0 | bins with overlap >= 0.25 | bins with overlap >= 0.50 |
|---|---:|---:|---:|
| PAR | 31 | 31 | 30 |
| sex homology | 31 | 31 | 30 |
| segmental duplication | 6,245 | 2,342 | 1,833 |
| low mappability | 414 | 148 | 100 |
| gap / centromere / telomere | 2,621 | 2,482 | 2,357 |
| repeat-rich | 28,786 | 171 | 92 |
| blacklist | 2,573 | 2,571 | 2,569 |
| ambiguous alignment | 0 | 0 | 0 |

This distribution is why repeat/segdup burden cannot become a broad hard event-level filter without locked ablation.

## CNVseq / maskpar Evidence

The 2026-06-11 probe reviewed older CNVseq resources:

- standard CNVseq uses `hg19_mask_par.fa.gz` and `uniqBin5k.bin2pos.gz`;
- AZF CNVseq uses `test_azf`, which is a custom binary mapper/index, not a normal hg19 FASTA;
- the generated AZF unique-bin complement excluded about 3.887% of the genome;
- chr6, chr9, chr11, and chr16 excluded fractions were low and had zero overlap with the then-current kept chr6/9/11/16 events.

Event-level findings from that probe remain important:

- H1 chr16 truth overlaps repeat/segdup soft mask substantially.
- Y5 chr16 false-positive-like events had a similar soft-mask burden pattern to H1 chr16 truth.
- H2 chr6 truth had low but nonzero soft/hard mask overlap.

Therefore:

- CNVseq unique-bin complement is useful as a reference asset to inspect, but current evidence does not justify promoting it to a production hard filter.
- `hg19_mask_par.fa.gz` / `maskpar` is useful for PAR and sex-homology interpretation, especially Branch S/SCA review.
- Full remapping to `hg19.mapability.fa.gz` or `maskpar` FASTA is not justified by current evidence.

## Decision Matrix

| option | decision | current use | promotion requirement |
|---|---|---|---|
| `annotation_only` | keep as default | Use PAR, sex-homology, low-map, segdup, repeat, blacklist, ref-bin support, and proximity fields as evidence annotations. | Unit/static checks plus evidence completeness checks. No downstream rerun needed if no report logic changes. |
| `reference_build_mask` | allowed only as shadow ref | Convert selected assets into reference-build masks without changing the alignment FASTA. | Build a new reference ID and rerun WisecondorX predict, Branch A candidates, Branch B, evaluation, benchmark, and report. Prove no known-positive recall regression. |
| `full_remap_reference` | not recommended now | Replace alignment reference with `hg19.mapability.fa.gz`, `hg19_mask_par.fa.gz`, or another altered FASTA. | Last resort only. Requires rebuilding aligner index, remapping BAMs, regenerating raw/mask-only NPZs, rebuilding reference, rerunning all downstream outputs, and proving no recall regression. |

Current Phase 1B decision:

```text
annotation_only remains the default.
reference_build_mask can be tested only as an explicitly named shadow reference.
full_remap_reference is not justified now.
```

## What Should Not Change Yet

Do not add any of the following based on current Phase 1B evidence:

- no Branch A artifact filter;
- no hard event-level filter based only on repeat, segmental duplication, low mappability, or blacklist overlap burden;
- no global PAR or sex-homology hard mask;
- no production reference promotion using H7-H16;
- no full hg19 FASTA replacement;
- no final-report change from Branch S or Branch B V2 shadow evidence.

## Why Event Blacklist Or Repeat Burden Cannot Be A Hard Standard

The bin size and event size make event-level overlap burden coarse. At 100 kb analysis bins, a broad event can contain a masked subregion and still have a consistent CNV signal on both sides. The previous chr16 and chr6 observations also show that known positives can carry the same soft-mask burden pattern as false-positive-like events.

Therefore, a region-risk burden should be used as context:

- if the event is mostly masked and clean-bin support is low, downgrade to artifact/review;
- if Branch A is strong and clean-bin support is not extremely low, keep as review rather than direct artifact;
- if the event is PAR/sex-homology related, route to sex-aware/Branch S review context;
- if same-chromosome ref leakage is observed, treat it as a method-contract violation and artifact-risk evidence.

This is a forward rule based on asset semantics and WisecondorX contract, not a result-backfilled threshold.

## Current Gaps

1. `mappability_score` in the current bin annotation is an ATGC-derived proxy unless an external mappability BED is configured. The low-mappability BED exists as interval evidence, but it is not a full per-base mappability track.
2. `ambiguous_alignment` is currently unpopulated in the active asset table.
3. `XTR` is currently zero under the configured region set.
4. H7-H16 are not N0 clean negatives. Phase 1B cannot use them to justify hard thresholds.
5. Any reference-build mask promotion needs a separate reference ID and full downstream rerun.

## Recommended Next Step

Keep Phase 1B as a documented decision matrix and proceed with Branch B V2 shadow evaluation. The next implementation-level step should be a shadow reference experiment only if the question is specifically whether a selected reference-build mask improves stability without recall regression.

For normal Branch B V2 iteration, use the current asset fields as evidence only and require candidate-level ablation before any field changes report disposition.
