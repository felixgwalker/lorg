# Ancestral Gene Content Reconstructor

Reconstructs the gene content of ancestral genomes at internal phylogenetic nodes using Dollo parsimony or a Bayesian gain/loss model.

## Overview

Given a gene presence/absence table and a species phylogeny, this tool infers which genes were present in each ancestral genome at every internal node, reporting the reconstructed gene set and posterior probabilities (Bayesian mode) or parsimony support (parsimony mode).

## Approach

**Inputs:** TSV of gene presence/absence per species (rows = genes, columns = species, values = 1/0); Newick species tree.

**Core method:** Under Dollo parsimony, each gene is assumed to have arisen once (on the most ancestral branch consistent with the observed presence pattern) and may have been lost multiple times. Under the Bayesian model, gain and loss rates are estimated by MCMC with flat priors; the posterior probability of presence at each node is sampled. Internal nodes are reconstructed as present (P > 0.8), absent (P < 0.2), or uncertain. Gene counts per ancestral node and branch-level gain/loss tallies are reported.

**Outputs:** TSV of ancestral gene states per node (`ancestral_gene_states.tsv`); ancestral genome size table (`ancestral_genome_sizes.tsv`); optional annotated phylogeny.

**How it ships:** `python run_reconstructor.py --presence-table genes.tsv --phylogeny tree.nwk`; `main.py` delegates to `src.pipeline.main()` which loads `run_reconstructor.py` via `importlib`.

## Usage

```bash
# Reconstruct ancestral gene content
python run_reconstructor.py --presence-table genes.tsv --phylogeny tree.nwk -o results/

# Synthetic demo (no real input required)
python run_reconstructor.py --demo -o results/

# Use Bayesian model
python run_reconstructor.py --presence-table genes.tsv --phylogeny tree.nwk --method bayesian -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0
- matplotlib>=3.7.0

## Status

Planned
