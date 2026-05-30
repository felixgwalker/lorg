# Population Viability Genomics Estimator

Estimates population viability by integrating inbreeding, effective population size, genetic load, and adaptive diversity to project extinction probability and derive the genomic minimum viable population size.

## Overview

Given a population VCF (and optionally census size), this tool computes current genomic metrics, models inbreeding accumulation under current Ne, estimates genetic load, and projects the probability of population persistence at 50, 100, and 200-year horizons using an individual-based stochastic model.

## Approach

**Inputs:** Population VCF; optional census size.

**Core method:** Current metrics are computed: FROH (from ROH), Ne (from LD), adaptive diversity score (Fst outlier heterozygosity). Genetic load is estimated from the ratio of rare deleterious (based on CADD/PolyPhen) to synonymous variants if annotations are present, or from the overall heterozygosity of low-frequency variants as a proxy. A Wright-Fisher model is run forward for each time horizon: each generation, F increases by 1/(2Ne) and inbreeding depression is applied as a multiplicative fitness reduction (parameterised by 3.14 lethal equivalents, the mean mammalian value). Extinction probability is estimated from 1000 stochastic replicates. The genomic minimum viable population is the Ne at which 95 % persistence probability is maintained over 100 years.

**Outputs:** TSV of viability projections (`viability_projections.tsv`); genomic metrics summary; optional extinction probability trajectory plot.

**How it ships:** `python run_estimator.py variants.vcf`; `main.py` delegates to `src.pipeline.main()` which loads `run_estimator.py` via `importlib`.

## Usage

```bash
# Estimate population viability
python run_estimator.py variants.vcf --census-size 500 -o results/

# Synthetic demo (no real input required)
python run_estimator.py --demo -o results/

# Custom time horizons
python run_estimator.py variants.vcf --time-horizons 50 100 500 -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0
- matplotlib>=3.7.0

## Status

Planned
