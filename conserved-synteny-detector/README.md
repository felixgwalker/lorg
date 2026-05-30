# Conserved Synteny Detector

Detects conserved synteny blocks between two genomes by chaining ortholog anchor pairs with a collinearity algorithm.

## Overview

Given ortholog gene pairs and gene position BEDs for two species, this tool identifies regions where gene order and orientation are conserved, reporting synteny blocks with their genomic coordinates, number of anchors, and orientation relative to both genomes.

## Approach

**Inputs:** TSV of ortholog gene pairs (gene_A, gene_B); BED of gene positions in species A; BED of gene positions in species B.

**Core method:** Ortholog pairs are mapped to their genomic coordinates in both genomes to form a set of anchor points. A collinearity scoring algorithm (similar to MCScan or i-ADHoRe) scans for windows of anchors that are co-linear within a tolerance of ≤ 10 intervening non-anchor genes. Anchor chains scoring above the collinearity cutoff are merged into synteny blocks. Orientation (same or inverted) is inferred from the strand of anchor pairs. Blocks below the minimum anchor count or length are discarded.

**Outputs:** TSV of synteny blocks (`synteny_blocks.tsv`); TSV of anchor pairs per block; optional dotplot.

**How it ships:** `python run_detector.py --orthologs orthologs.tsv --positions-a speciesA.bed --positions-b speciesB.bed`; `main.py` delegates to `src.pipeline.main()` which loads `run_detector.py` via `importlib`.

## Usage

```bash
# Detect synteny blocks between two species
python run_detector.py --orthologs orthologs.tsv --positions-a human.bed --positions-b mouse.bed -o results/

# Synthetic demo (no real input required)
python run_detector.py --demo -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- matplotlib>=3.7.0

## Status

Planned
