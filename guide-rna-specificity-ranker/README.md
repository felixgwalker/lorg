# Guide RNA Specificity Ranker

Ranks sgRNAs by predicted genome-wide specificity using CFD scoring against enumerated off-target sites.

## Overview

Off-target editing is a primary safety concern in therapeutic CRISPR applications. This tool scores each
guide RNA by the Cutting Frequency Determination (CFD) method — enumerating all genomic positions with
up to 3 mismatches and computing a weighted penalty sum — then ranks guides from most to least specific.

## Approach

**Inputs:** FASTA or TSV of 20 nt guide spacers; either a reference genome FASTA for de novo off-target
search, or a pre-computed off-target BED (e.g. from Cas-OFFinder, CRISPOR).

**Core method:** CFD scoring (Doench et al. 2016): each mismatched position is assigned a penalty from
the published rN:dN mismatch penalty matrix; the product of (1 − penalty) across all mismatches gives
the CFD score for each off-target site. The guide specificity score is 1 − Σ(off-target CFDs), capped
at 0. Guides are ranked by specificity score and off-target site count at ≤1, ≤2, ≤3 mismatches.

**Outputs:** Ranked TSV with specificity score, off-target counts, and band (poor/moderate/good/
excellent); optional off-target distribution plot.

**Dependencies reused:** biopython for FASTA I/O; numpy for CFD matrix operations.

**How it ships:** `python run_ranker.py guides.fa --genome genome.fa`; delegated from
`main.py → src.pipeline.main() → run_ranker.py`.

## Usage

```bash
python run_ranker.py guides.fa --genome genome.fa -o results/
python run_ranker.py guides.fa --offtargets sites.bed -o results/
python run_ranker.py --demo -o results/
```

## Dependencies

- biopython>=1.81
- pandas>=2.0.0
- numpy>=1.24.0
- matplotlib>=3.7.0

## Status

Planned
