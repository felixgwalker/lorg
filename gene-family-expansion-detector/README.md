# Gene Family Expansion Detector

Detects gene families that have expanded significantly on specific lineages using a birth-death model.

## Overview

Given a gene family size table (genes × species) and a species phylogeny, this tool fits a birth-death model genome-wide (CAFE-style), estimates the per-family rate of gain/loss, and tests each family for significant expansion on individual branches using a likelihood ratio test against the background rate.

## Approach

**Inputs:** TSV of gene family sizes per species (rows = family IDs, columns = species); Newick species tree with branch lengths.

**Core method:** A global birth-death rate (λ) is estimated by maximum likelihood across all families and branches. For each family, an independent branch-specific rate is tested; families with a significantly higher rate on a branch than the background rate (FDR-corrected LRT p < threshold) are flagged as expansions. Fold change is computed as observed size / expected size under the background rate. Functional annotations (if provided) are reported alongside each expansion.

**Outputs:** TSV of significant family expansions (`gene_family_expansions.tsv`); optional bubble chart of fold-change vs. p-value.

**How it ships:** `python run_detector.py --family-table families.tsv --phylogeny tree.nwk`; `main.py` delegates to `src.pipeline.main()` which loads `run_detector.py` via `importlib`.

## Usage

```bash
# Detect gene family expansions
python run_detector.py --family-table families.tsv --phylogeny tree.nwk -o results/

# Synthetic demo (no real input required)
python run_detector.py --demo -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0
- matplotlib>=3.7.0

## Status

Planned
