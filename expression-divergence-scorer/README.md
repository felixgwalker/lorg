# Expression Divergence Scorer

Scores expression divergence between orthologous gene pairs across two species using log fold change, Jensen-Shannon Index, tau, or Euclidean distance.

## Overview

Given matched expression matrices for two species and an ortholog table, this tool computes a divergence score for each ortholog pair across tissues or conditions, identifying genes with significantly higher or lower expression divergence than expected by chance.

## Approach

**Inputs:** Gene expression matrices (genes × tissues, TPM or normalised counts) for species A and B; TSV of ortholog gene pairs; matched tissue/condition labels.

**Core method:** For each ortholog pair, the expression vectors across matched tissues are extracted and normalised. The chosen metric is applied: **log fold change** — mean absolute log₂(TPM_A / TPM_B) across tissues; **JSI** — Jensen-Shannon Index of the expression proportion vectors; **tau** — tissue specificity index per species, then delta-tau as divergence; **Euclidean** — L2 distance of normalised expression vectors. Statistical significance is assessed by permuting tissue labels (1000 reps). FDR correction (Benjamini-Hochberg) is applied across gene pairs.

**Outputs:** TSV of divergence scores and p-values (`expression_divergence.tsv`); optional scatter plot of divergence vs. evolutionary distance.

**How it ships:** `python run_scorer.py --expression-a speciesA.tsv --expression-b speciesB.tsv --orthologs orthologs.tsv`; `main.py` delegates to `src.pipeline.main()` which loads `run_scorer.py` via `importlib`.

## Usage

```bash
# Score expression divergence between two species
python run_scorer.py --expression-a speciesA.tsv --expression-b speciesB.tsv --orthologs orthologs.tsv -o results/

# Synthetic demo (no real input required)
python run_scorer.py --demo -o results/

# Use Jensen-Shannon Index
python run_scorer.py --expression-a speciesA.tsv --expression-b speciesB.tsv --orthologs orthologs.tsv --metric JSI -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0
- matplotlib>=3.7.0

## Status

Planned
