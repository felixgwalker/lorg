# Genomic Diversity Index

Computes a comprehensive genomic diversity index (θW, θπ, Tajima's D, Ho, He, Fis) in sliding windows across the genome.

## Overview

Given a population VCF, this tool calculates classical nucleotide diversity statistics in sliding windows and produces a genome-wide summary index that quantifies the overall level and distribution of genetic diversity in the population.

## Approach

**Inputs:** VCF of bi-allelic SNPs (single- or multi-sample for Ho/He/Fis metrics).

**Core method:** For each window meeting the minimum SNP count threshold: Watterson's θW is computed from the number of segregating sites; θπ from average pairwise differences; Tajima's D from D = (θπ − θW) / sqrt(Var(θπ − θW)); observed heterozygosity (Ho) from the fraction of heterozygous genotypes; expected heterozygosity (He = 2pq summed across sites); and inbreeding coefficient Fis = 1 − Ho/He. Genome-wide values are computed as weighted means across windows. Windows with Tajima's D outside (−2, 2) are flagged as potentially under selection or demographic influence.

**Outputs:** Per-window TSV of all metrics (`diversity_windows.tsv`); genome-wide summary TSV (`diversity_index.tsv`); optional genome-wide diversity track plot.

**How it ships:** `python run_index.py variants.vcf`; `main.py` delegates to `src.pipeline.main()` which loads `run_index.py` via `importlib`.

## Usage

```bash
# Compute all diversity metrics
python run_index.py variants.vcf -o results/

# Synthetic demo (no real input required)
python run_index.py --demo -o results/

# Compute only θW and θπ with a 100 kb window
python run_index.py variants.vcf --metrics theta_w theta_pi --window-size 100000 -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0
- matplotlib>=3.7.0

## Status

Planned
