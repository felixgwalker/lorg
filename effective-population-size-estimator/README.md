# Effective Population Size Estimator

Estimates effective population size (Ne) from SNP data using Watterson's estimator, Tajima's estimator, LD-based Ne, or a PSMC-style temporal profile.

## Overview

Given a VCF of SNPs, this tool computes the site frequency spectrum (SFS), derives nucleotide diversity statistics (θW, θπ), and converts them to Ne using a user-supplied mutation rate and generation time. An LD-based method and a simplified PSMC-style temporal Ne trajectory are also supported.

## Approach

**Baseline:** NeEstimator (Do et al. 2014) and PSMC (Li & Durbin 2011) estimate
Ne from contemporary and diploid-heterozygosity data respectively. This tool does
not re-implement their algorithms — it calls NeEstimator and provides a PSMC-style
HMM that can call samtools/bcftools for input preparation.

**Novel layer:** Kept *only* for the temporal/aDNA Ne angle. The de-extinction-specific
additions: (1) ancient-tip sampling — Ne estimation from a time series of ancient
samples (pseudo-haploid calls at low coverage) using a temporal-samples extension
of the LD-Ne estimator; (2) low-coverage correction for heterozygosity bias; (3)
output framed as a Ne trajectory feeding `genetic-rescue-viability-estimator` and
`inbreeding-risk-forecaster` rather than as a standalone statistic. The parts of
this tool that duplicate NeEstimator on modern diploid data should not be published.

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
