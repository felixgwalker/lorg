# Lineage Divergence Dater

Dates lineage divergence events using a molecular clock, Bayesian dated-tips approach for ancient DNA, or pairwise genetic distance.

## Overview

Given a sequence alignment, phylogeny, and optional calibration constraints or radiocarbon-dated tip ages, this tool estimates split times between lineages with confidence intervals, supporting both modern and ancient sample sets.

## Approach

**Inputs:** FASTA multiple sequence alignment; Newick phylogenetic tree; optional JSON of calibration points (node ID, age bounds) or dated-tip ages (sample ID, radiocarbon age BP, uncertainty).

**Core method:** (1) **Molecular clock** — calibrated node dating using a strict or relaxed clock, converting branch lengths to absolute time using the calibration points and a substitution rate. (2) **Bayesian dated tips** — for ancient DNA datasets, tip ages from radiocarbon dating are treated as calibrations, allowing the substitution rate and divergence dates to be co-estimated; useful when no fossil calibrations are available but multiple ancient samples are dated. (3) **Pairwise distance** — simple Jukes-Cantor corrected pairwise distances divided by 2× the substitution rate, providing quick uncalibrated estimates. All methods output mean divergence time and 95 % confidence/credible intervals.

**Outputs:** TSV of divergence dates per lineage pair (`divergence_dates.tsv`); optional time-calibrated phylogeny and timeline plot.

**How it ships:** `python run_dater.py alignment.fa --phylogeny tree.nwk --calibrations calibrations.json`; `main.py` delegates to `src.pipeline.main()` which loads `run_dater.py` via `importlib`.

## Usage

```bash
# Date lineage divergences with calibrations
python run_dater.py alignment.fa --phylogeny tree.nwk --calibrations calibrations.json -o results/

# Synthetic demo (no real input required)
python run_dater.py --demo -o results/

# Ancient DNA dated-tips approach
python run_dater.py alignment.fa --phylogeny tree.nwk --calibrations tip_ages.json --method bayesian_dated_tips -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0
- matplotlib>=3.7.0

## Status

Planned
