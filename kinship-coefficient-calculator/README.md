# Kinship Coefficient Calculator

Calculates pairwise kinship coefficients for all samples in a multi-sample VCF and classifies relationships by degree.

## Overview

Given a multi-sample VCF of LD-pruned SNPs, this tool computes pairwise kinship using the KING-robust estimator (or genomic relatedness matrix / IBD segments) and assigns each pair to a relationship class: identical twins, 1st-degree, 2nd-degree, 3rd-degree, or unrelated.

## Approach

**Inputs:** Multi-sample VCF of LD-pruned bi-allelic SNPs (recommended ≥ 10 k SNPs after pruning).

**Core method:** The KING-robust kinship estimator computes φ̂ = (IBS2 − 2·IBS0) / (4 · sqrt(Het_A · Het_B)) per pair, which is robust to population structure. IBD-mode estimates the probability of sharing 0, 1, or 2 alleles IBD (π0, π1, π2) from allele-frequency-adjusted IBS sharing. The genomic relatedness matrix (GRM) computes A = (1/L) · Σ (g − 2p)(g − 2p)ᵀ / 2p(1−p). Relationship thresholds follow KING conventions: φ ≥ 0.354 (identical), ≥ 0.177 (1st degree), ≥ 0.0884 (2nd degree), ≥ 0.0442 (3rd degree).

**Outputs:** TSV of pairwise kinship coefficients above the threshold (`kinship_pairs.tsv`); optional heatmap of the kinship matrix.

**How it ships:** `python run_calculator.py samples.vcf`; `main.py` delegates to `src.pipeline.main()` which loads `run_calculator.py` via `importlib`.

## Usage

```bash
# Calculate kinship for all sample pairs
python run_calculator.py samples.vcf -o results/

# Synthetic demo (no real input required)
python run_calculator.py --demo -o results/

# Use IBD method
python run_calculator.py samples.vcf --method IBD -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0
- matplotlib>=3.7.0

## Status

Planned
