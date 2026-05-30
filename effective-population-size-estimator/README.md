# Effective Population Size Estimator

Estimates effective population size (Ne) from SNP data using Watterson's estimator, Tajima's estimator, LD-based Ne, or a PSMC-style temporal profile.

## Overview

Given a VCF of SNPs, this tool computes the site frequency spectrum (SFS), derives nucleotide diversity statistics (θW, θπ), and converts them to Ne using a user-supplied mutation rate and generation time. An LD-based method and a simplified PSMC-style temporal Ne trajectory are also supported.

## Approach

**Inputs:** VCF of bi-allelic SNPs (mono- or multi-sample).

**Core method:** The folded SFS is computed from allele count data. Watterson's θW uses the number of segregating sites; Tajima's θπ uses average pairwise differences. Both are converted to Ne via θ = 4Neμ. The LD-based estimator computes r² between SNP pairs within a sliding window; Ne is estimated from the expected r² under drift. The PSMC-style method uses a hidden Markov model on heterozygosity density to infer a temporal Ne trajectory. All methods return point estimates with bootstrap confidence intervals.

**Outputs:** TSV of Ne estimates per method (`ne_estimates.tsv`); optional SFS and temporal trajectory plots.

**How it ships:** `python run_estimator.py variants.vcf`; `main.py` delegates to `src.pipeline.main()` which loads `run_estimator.py` via `importlib`.

## Usage

```bash
# Estimate Ne from a real VCF
python run_estimator.py variants.vcf -o results/

# Synthetic demo (no real input required)
python run_estimator.py --demo -o results/

# Use LD-based method with custom generation time
python run_estimator.py variants.vcf --method ld_based --generation-time 25 -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0
- matplotlib>=3.7.0

## Status

Planned
