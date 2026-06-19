# CNVpro Standard Check

Date: 2026-06-19

## Purpose

Confirm the CNVpro/CNVseq rules before using them as Branch A/B V2 evidence.
This report corrects the earlier loose wording around chr13/14/15 p-arm or
high-repeat residual signals.

## Sources Read

Remote source locations:

- `fengxian:/home/jiucheng/project/cnvpro/cnvpro_minimal_20260614/code/cnvseqpipe/CNVcalling.R`
- `fengxian:/home/jiucheng/project/cnvpro/cnvpro_minimal_20260614/code/cnvseqpipe/filterCNV.py`
- `fengxian:/home/jiucheng/project/cnvpro/cnvpro_minimal_20260614/code/cnvseqpipe/GCcorrect.py`
- `fengxian:/home/jiucheng/project/cnvpro/cnvpro_minimal_20260614/code/cnvseqpipe/cytoBand.txt.gz`
- `BS:/sg2/7.hudahui/CNV/cnvpropipe/README.md`

## Confirmed CNVpro Rules

### README-level project rules

`BS:/sg2/7.hudahui/CNV/cnvpropipe/README.md` states:

- CNV filtering is similar to CNV-seq.
- Keep CNVs above 50 kb.
- Keep 20% mosaic aneuploidy.
- Exclude mosaic CNVs below 1 Mb.
- The example config uses `hg19.mapability.fa.gz` as the reference genome and
  `block=6`, described as 5 kb times 6 for CNV calling.

### Acrocentric chromosome segmentation

`CNVcalling.R` handles acrocentric chromosomes explicitly:

- chr13/14/15/21/22 are segmented only on `qter`.
- Other autosomes are segmented as `pter` and `qter`.
- chrX is segmented into PAR1, pter, qter, and PAR2.

This means the source CNVpro logic does not merely annotate chr13/14/15 p-arm
risk after calling. It avoids standard segmentation on the acrocentric short
arms by design.

### Cytoband context

`cytoBand.txt.gz` marks chr13/14/15 short-arm regions with `gvar`, `stalk`, and
`acen` categories. These regions are used by CNVpro for cytoband context and
plotting. They should not be converted into a PGT-A hard filter unless a
separate ablation proves no known-positive recall regression.

### filterCNV thresholds

`filterCNV.py` applies:

- copy-number thresholds: default `cnupper=2.3`, `cnlower=1.7`;
- mosaic thresholds: default `mosupper=2.7`, `moslower=1.3`;
- length threshold: default `length=100,000`;
- mosaic event handling that skips mosaic CNVs below 1 Mb;
- whole-chromosome / aneuploidy detection when same-direction CNV markers cover
  at least `chroverlap=0.9` of the chromosome marker count;
- sex-aware PAR / non-PAR handling;
- AnnotSV rank routing: events with rank above 2 are primary positive; lower
  rank events are retained as other/negative annotation context.

There is no source evidence in the inspected `filterCNV.py` for a dedicated
chr13/14/15 p-arm hard artifact filter.

### GC / dynamic reference

`GCcorrect.py` applies LOESS GC correction and can optionally apply dynamic
reference normalization using nearest controls. It separates autosome, chrX,
and chrY bins for dynamic-reference normalization.

## Implication For PGT-A Branch A/B V2

Use CNVpro as evidence in these ways:

1. Treat acrocentric chr13/14/15/21/22 qter-only segmentation as a verified
   method-contract idea to evaluate.
2. Do not write “chr13/14/15 p-arm high-repeat contribution” as a confirmed
   CNVpro standard unless it is tied to the qter-only segmentation contract or
   cytoband asset evidence.
3. Do not add a Branch A hard filter for chr13/14/15 p-arms.
4. If acrocentric p-arm masking is tested, it must be a named shadow
   reference-build mask experiment with full downstream rerun.
5. Keep CNVpro thresholds as priors for ablation, not production defaults.

## Documentation Correction

The Branch A post-reference rebuild check should say:

```text
Evaluate whether residual signals fall in regions explicitly handled by source
CNVpro logic. Current source review shows CNVpro segments chr13/14/15/21/22
only on qter; this is an acrocentric segmentation-contract check, not a guessed
chr13/14/15 p-arm hard filter.
```
