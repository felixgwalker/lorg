# Admixture Signal Scanner

Scans for admixture signals and estimates per-sample ancestry proportions by fitting parametric admixture models across a range of K values.

## Overview

Given a VCF of LD-pruned SNPs, this tool fits unsupervised (or supervised) admixture models for K = k_min…k_max using an EM algorithm, selects the best K by cross-validation error, and reports per-sample ancestry proportions in a structure-style output with optional bar chart visualisation.

## Approach

**Inputs:** VCF of LD-pruned bi-allelic SNPs (use `--ld-prune` upstream or provide pre-pruned data); optional reference population panel VCF for supervised mode.

**Core method:** An EM algorithm iterates between updating individual ancestry proportions (Q matrix) and population allele frequency estimates (P matrix) until convergence. Cross-validation is performed by masking 20 % of genotypes per fold and computing prediction error. The K with lowest CV error is selected as the best model. Unsupervised mode initialises randomly (seed-controlled); supervised mode fixes reference population frequencies. Admixed individuals are flagged when no single component exceeds 90 %.

**Outputs:** Q matrix TSV of ancestry proportions (`ancestry_proportions.tsv`); P matrix TSV of allele frequencies per component; CV error TSV; optional structure bar chart.

**How it ships:** `python run_scanner.py variants.vcf`; `main.py` delegates to `src.pipeline.main()` which loads `run_scanner.py` via `importlib`.

## Usage

```bash
# Scan admixture for K=2..8
python run_scanner.py variants.vcf --k-min 2 --k-max 8 -o results/

# Synthetic demo (no real input required)
python run_scanner.py --demo -o results/

# Supervised mode with reference panel
python run_scanner.py variants.vcf --model supervised --reference-panel ref.vcf -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0
- matplotlib>=3.7.0

## Status

Planned
