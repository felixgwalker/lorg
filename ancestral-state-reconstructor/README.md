# Ancestral State Reconstructor

Reconstructs ancestral nucleotide or amino acid sequences at internal phylogenetic nodes using marginal maximum likelihood, parsimony, or Bayesian inference.

## Overview

Given a multiple sequence alignment and a phylogenetic tree, this tool computes the most probable ancestral state at each site and each internal node, outputting full ancestral FASTA sequences and per-site posterior probabilities.

## Approach

**Inputs:** FASTA multiple sequence alignment; Newick phylogenetic tree.

**Core method:** Under the **ML** method, Felsenstein's pruning algorithm computes the marginal likelihood of each character state at each internal node given the alignment and a specified substitution model (default GTR+G with 4 gamma rate categories). The state with the highest posterior is assigned; sites with maximum posterior < 0.5 are flagged as uncertain. **Parsimony** uses Fitch parsimony to assign the most parsimonious ancestral state with no branch-length model. **Bayesian** samples ancestral states under the GTR+G model via MCMC, summarising posterior distributions per site.

**Outputs:** FASTA of ancestral sequences per internal node (`ancestral_sequences.fa`); TSV of per-site posterior probabilities (`site_posteriors.tsv`); optional phylogeny with uncertain sites highlighted.

**How it ships:** `python run_reconstructor.py alignment.fa --phylogeny tree.nwk`; `main.py` delegates to `src.pipeline.main()` which loads `run_reconstructor.py` via `importlib`.

## Usage

```bash
# Reconstruct ancestral sequences (ML method)
python run_reconstructor.py alignment.fa --phylogeny tree.nwk -o results/

# Synthetic demo (no real input required)
python run_reconstructor.py --demo -o results/

# Use Bayesian method
python run_reconstructor.py alignment.fa --phylogeny tree.nwk --method bayesian -o results/
```

## Dependencies

- biopython>=1.81
- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0
- matplotlib>=3.7.0

## Status

Planned
