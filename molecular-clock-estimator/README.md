# Molecular Clock Estimator

Estimates the molecular clock rate and compares strict vs. relaxed clock models using Bayesian MCMC with fossil/biogeographic calibration points.

## Overview

Given a multiple sequence alignment, a starting phylogeny, and calibration constraints, this tool fits strict and relaxed lognormal/exponential clock models via Bayesian MCMC, selects the best model by Bayes factors, and reports per-branch substitution rates with credible intervals.

## Approach

**Inputs:** FASTA multiple sequence alignment; Newick starting tree; optional JSON of calibration constraints (node ID, minimum age, maximum age, prior distribution).

**Core method:** A Bayesian hierarchical model is set up with a prior on the overall substitution rate (gamma or log-normal) and a tree prior (Yule or birth-death). Under the strict clock, all branches share a single rate. Under the relaxed lognormal/exponential clock, each branch draws its rate independently from a distribution. MCMC samples from the joint posterior; ESS and convergence diagnostics are checked automatically. Bayes factors (from path sampling or stepping-stone) compare strict vs. relaxed models.

**Outputs:** MCMC log file; TSV of rate estimates per model (`clock_estimates.tsv`); optional trace plot and consensus tree.

**How it ships:** `python run_estimator.py alignment.fa --phylogeny tree.nwk --calibrations calibrations.json`; `main.py` delegates to `src.pipeline.main()` which loads `run_estimator.py` via `importlib`.

## Usage

```bash
# Estimate molecular clock with calibrations
python run_estimator.py alignment.fa --phylogeny tree.nwk --calibrations calibrations.json -o results/

# Synthetic demo (no real input required)
python run_estimator.py --demo -o results/

# Strict clock only
python run_estimator.py alignment.fa --phylogeny tree.nwk --clock-model strict -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0
- matplotlib>=3.7.0

## Status

Planned
