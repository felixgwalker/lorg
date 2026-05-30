# Inbreeding Risk Forecaster

Forecasts inbreeding risk for a population by estimating current inbreeding from runs of homozygosity, inferring Ne, and projecting future inbreeding accumulation.

## Overview

Given a population VCF, this tool estimates FROH (from ROH), Fis (from heterozygosity), derives effective population size, and projects the expected increase in inbreeding over a user-defined number of generations to classify current and future risk.

## Approach

**Inputs:** VCF of SNPs for the target population.

**Core method:** Per-individual ROH are detected using a sliding-window approach (default ≥ 500 kb). FROH = Σ ROH_length / autosome_length. Fis is computed from observed vs. expected heterozygosity. Ne is estimated from the LD-based Ne estimator or from FROH time series if multiple generations are available. Future inbreeding is projected using the Wright formula: ΔF per generation = 1/(2Ne). Risk thresholds: critical (F > 0.25, Ne < 50); high (F > 0.125, Ne < 100); moderate (F > 0.0625, Ne < 200); low otherwise.

**Outputs:** TSV of per-individual ROH and FROH (`inbreeding_stats.tsv`); forecast TSV; optional Ne trajectory and inbreeding projection plot.

**How it ships:** `python run_forecaster.py variants.vcf`; `main.py` delegates to `src.pipeline.main()` which loads `run_forecaster.py` via `importlib`.

## Usage

```bash
# Forecast inbreeding risk
python run_forecaster.py variants.vcf -o results/

# Synthetic demo (no real input required)
python run_forecaster.py --demo -o results/

# Project over 100 generations with generation time 2 years
python run_forecaster.py variants.vcf --n-generations 100 --generation-time 2.0 -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0
- matplotlib>=3.7.0

## Status

Planned
