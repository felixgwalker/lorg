# Coexpression Module Finder

Identifies coexpression modules from RNA-seq data using WGCNA, clique-based clustering, or k-means, with optional module-trait correlation.

## Overview

Given a normalised expression matrix, this tool constructs a gene coexpression network and clusters genes into modules of similar expression patterns, identifying hub genes and optionally correlating module eigengenes with sample phenotypes/traits.

## Approach

**Inputs:** Normalised gene expression matrix (genes × samples, log₂ TPM or VST-normalised counts); optional sample trait TSV (samples × traits) for module-trait correlation.

**Core method:** Under **WGCNA**: a signed weighted adjacency matrix is computed by raising absolute Pearson correlations to a soft-thresholding power (default 6) to approximate scale-free topology; topological overlap matrix (TOM) is derived; hierarchical clustering with dynamic tree cutting identifies modules; small modules (< `min_module_size`) are merged. Module eigengenes (first PC of each module's expression matrix) are computed and correlated with traits if provided. **Clique-based**: genes with Pearson r ≥ 0.8 to each other form dense subgraphs. **k-means**: expression vectors are k-means clustered after dimensionality reduction. Hub genes are the genes with highest intra-module connectivity (kIM score).

**Outputs:** TSV of gene module assignments (`module_assignments.tsv`); module eigengene matrix; optional heatmap and module-trait correlation plot.

**How it ships:** `python run_finder.py --expression expression.tsv`; `main.py` delegates to `src.pipeline.main()` which loads `run_finder.py` via `importlib`.

## Usage

```bash
# Find coexpression modules
python run_finder.py --expression expression.tsv -o results/

# Synthetic demo (no real input required)
python run_finder.py --demo -o results/

# Include trait correlation
python run_finder.py --expression expression.tsv --trait-data traits.tsv -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0
- matplotlib>=3.7.0

## Status

Planned
