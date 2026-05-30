# CRISPR Knockout Designer

Designs sgRNAs for gene knockout, targeting early coding exons with Doench 2016-style on-target scoring.

## Overview

For gene knockout experiments, the most effective guides target the first one to three coding exons,
where frameshift-inducing indels reliably ablate protein function. This tool scans the gene sequence
for NGG PAM sites in early exons, scores each candidate guide, and returns the top-N designs ranked
by predicted on-target efficiency and frameshift probability.

## Approach

**Inputs:** Gene or CDS FASTA (full gene or CDS only); Cas9 PAM motif (default NGG); number of top
guides to return.

**Core method:** The tool targets the first 33 % of the CDS (or exons 1–3 when annotation is provided).
Each NGG PAM site on either strand yields a 20 nt spacer scored by Doench 2016-style features: GC
content (30–70 % optimal), seed-region GC (positions 1–12 from PAM), penalise A at position 20 (weak
binding), G clamp at position 1 (strong initiation), homopolymer runs, and an approximate thermodynamic
accessibility score. Frameshift probability is estimated as the fraction of inDelphi-predicted outcomes
that cause a +1 or −2 indel (the two most common frameshift-inducing outcomes).

**Outputs:** Ranked guide TSV with on-target score, predicted frameshift rate, and position; optional
gene schematic plot showing guide locations.

**Dependencies reused:** biopython for FASTA I/O; numpy for scoring.

**How it ships:** `python run_designer.py gene.fa`; delegated from
`main.py → src.pipeline.main() → run_designer.py`.

## Usage

```bash
python run_designer.py gene.fa -o results/
python run_designer.py gene.fa --n-guides 10 -o results/
python run_designer.py --demo -o results/
```

## Dependencies

- biopython>=1.81
- pandas>=2.0.0
- numpy>=1.24.0
- matplotlib>=3.7.0

## Status

Planned
