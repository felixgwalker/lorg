# Cas Variant Selector

Ranks Cas nucleases and editors for a target locus and desired editing goal.

## Overview

With a rapidly expanding toolkit of Cas variants, selecting the optimal nuclease for a given locus and
editing application requires balancing PAM availability, protein size (for delivery constraints), and
editing window compatibility. This tool scores a curated panel of ten Cas variants against these
criteria and returns a ranked recommendation.

## Approach

**Inputs:** Target locus FASTA; editing goal selection (knockout, base-edit, prime-edit, activation,
repression).

**Core method:** For each Cas variant in the database (SpCas9, SaCas9, Cas9-NG, SpRY, AsCas12a,
LbCas12a, Cas12b, CasX, ABE8e, CBE4max, PE2): (1) PAM density is computed by IUPAC scanning of the
locus. (2) Goal compatibility score: for knockout — all DSB-capable Cas variants score 1.0; for
base-edit — only deaminase fusions score 1.0, DSB variants score 0.1; for prime-edit — only PE2
variants score 1.0. (3) Delivery score: AAV-compatible variants (SaCas9, LbCas12a, Cas12b, CasX)
score higher when AAV delivery is implied. Composite score is the sum of weighted components.

**Outputs:** Ranked variant TSV with per-criterion scores and recommendation notes; optional bar chart.

**How it ships:** `python run_selector.py locus.fa --goal knockout`; delegated from
`main.py → src.pipeline.main() → run_selector.py`.

## Usage

```bash
python run_selector.py locus.fa --goal knockout -o results/
python run_selector.py locus.fa --goal base-edit -o results/
python run_selector.py --demo -o results/
```

## Dependencies

- biopython>=1.81
- pandas>=2.0.0
- numpy>=1.24.0
- matplotlib>=3.7.0

## Status

Planned
