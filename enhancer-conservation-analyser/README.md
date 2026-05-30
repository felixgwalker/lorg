# Enhancer Conservation Analyser

Analyses evolutionary conservation of enhancer elements using phastCons/PhyloP scores and multiple alignment coverage.

## Overview

Given a BED of enhancer elements and conservation data (phastCons/PhyloP bigWig or MAF alignment), this tool retrieves per-element conservation scores and classifies each enhancer as highly conserved, moderately conserved, or lineage-specific.

## Approach

**Inputs:** BED of enhancer elements; optional phastCons or PhyloP bigWig; optional multiple alignment MAF.

**Core method:** For each enhancer, the mean phastCons score is extracted from the bigWig. If a MAF alignment is provided, the fraction of alignment columns covered across vertebrate species is also computed. Enhancers with mean phastCons ≥ 0.8 (or PhyloP ≥ 2.0) and coverage in ≥ 10 species are classified as highly conserved; 0.4–0.8 as moderately conserved; below 0.4 or present in < 20 % of species as lineage-specific. Functional conservation can be optionally assessed by checking if conserved orthologs also have ATAC/H3K27ac signal in matched cell types.

**Outputs:** TSV of conservation scores per element (`enhancer_conservation.tsv`); optional bar chart of conservation class distribution.

**How it ships:** `python run_analyser.py --enhancers enhancers.bed --conservation phastcons.bw`; `main.py` delegates to `src.pipeline.main()` which loads `run_analyser.py` via `importlib`.

## Usage

```bash
# Analyse enhancer conservation
python run_analyser.py --enhancers enhancers.bed --conservation phastcons.bw -o results/

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
