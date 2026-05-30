# Regulatory Rewiring Analyser

Analyses regulatory rewiring between two species by classifying enhancers and promoters near orthologous genes as gained, lost, conserved, or relocated.

## Overview

Given BED files of regulatory elements (e.g., from ATAC-seq or ChIP-seq peaks) for two species and an ortholog table, this tool matches elements near orthologous genes, compares their sequence conservation via genome alignment, and classifies each element's evolutionary fate.

## Approach

**Inputs:** BED of regulatory elements in species A (e.g., open chromatin peaks); BED of regulatory elements in species B; TSV of ortholog gene pairs.

**Core method:** For each orthologous gene pair, regulatory elements within the promoter window (default 2 kb upstream + 500 bp downstream) are retrieved in both species. Each species-A element is searched for a sequence-conserved counterpart in species B using the pairwise genome alignment (or sequence liftover); elements with identity ≥ 70 % are called conserved; those below 30 % but with a synteny-supported homologous position are called relocated; absent elements are called lost in B or gained in B depending on direction. TFBS content is compared between conserved pairs to flag TF binding changes.

**Outputs:** TSV of regulatory rewiring events (`regulatory_rewirings.tsv`); optional bar chart of rewiring type counts.

**How it ships:** `python run_analyser.py --elements-a speciesA.bed --elements-b speciesB.bed --orthologs orthologs.tsv`; `main.py` delegates to `src.pipeline.main()` which loads `run_analyser.py` via `importlib`.

## Usage

```bash
# Analyse regulatory rewiring between two species
python run_analyser.py --elements-a human_atac.bed --elements-b mouse_atac.bed --orthologs orthologs.tsv -o results/

# Synthetic demo (no real input required)
python run_analyser.py --demo -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0
- matplotlib>=3.7.0

## Status

Planned
