# Gene Essentiality Predictor

Predicts gene essentiality by aggregating DepMap CRISPR fitness scores, RNAi dependency data, and gnomAD constraint metrics into a composite classification.

## Overview

Given a list of genes and optional DepMap/gnomAD data, this tool combines CRISPR gene effect scores, RNAi dependency scores, and population constraint metrics (LOEUF, pLI) into a weighted composite essentiality score, classifying each gene as essential, context-dependent, or non-essential.

## Approach

**Inputs:** Text file of gene IDs or symbols; optional DepMap CRISPR gene effect TSV; optional gnomAD constraint TSV.

**Core method:** Evidence is collected per gene: (1) **CRISPR fitness score** (DepMap CERES/Chronos) — scores < −1.0 indicate essentiality in ≥ 50 % of cell lines; (2) **RNAi dependency** — DEMETER2 scores < −1.5 indicate essentiality; (3) **LOEUF** — LOEUF < 0.35 suggests haploinsufficiency/constraint; (4) **pLI** — pLI > 0.9 indicates intolerance to LoF variants. Each evidence type is weighted and summed into a composite score [0, 1]. Score ≥ 0.7 → essential; 0.4–0.7 → context-dependent; < 0.4 → non-essential. Context-specific essentiality is reported for a user-supplied cell line.

**Outputs:** TSV of per-gene essentiality predictions (`essentiality_predictions.tsv`); optional violin plot of composite score distributions.

**How it ships:** `python run_predictor.py --gene-list genes.txt --depmap-scores depmap.tsv`; `main.py` delegates to `src.pipeline.main()` which loads `run_predictor.py` via `importlib`.

## Usage

```bash
# Predict gene essentiality
python run_predictor.py --gene-list genes.txt --depmap-scores depmap.tsv --constraint-data gnomad.tsv -o results/

# Synthetic demo (no real input required)
python run_predictor.py --demo -o results/

# Context-specific for a DepMap cell line
python run_predictor.py --gene-list genes.txt --depmap-scores depmap.tsv --cell-line ACH-000001 -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0
- matplotlib>=3.7.0

## Status

Planned
