# Metabolic Pathway Balancer

Balances a heterologous metabolic pathway for optimal flux distribution using FBA, identifying bottlenecks, cofactor imbalances, and toxic intermediates.

## Overview

Given a pathway definition (reactions, stoichiometry, kinetics), this tool runs flux balance analysis, detects limiting steps, and recommends enzyme expression level adjustments to maximise product yield, productivity, or minimise by-products.

## Approach

**Inputs:** JSON of pathway reactions with substrates, products, stoichiometric coefficients, and optional Km/kcat parameters.

**Core method:** The stoichiometric matrix is assembled from all reactions. FBA is run with the specified objective by solving the linear programme: maximise c·v subject to S·v = 0 (steady state), lb ≤ v ≤ ub. Flux variability analysis (FVA) identifies reactions with constrained flux ranges — low-flux steps are flagged as bottlenecks. Cofactor balance (NADH, NADPH, ATP) is checked; imbalances that would limit product formation are reported. Intermediate accumulation is estimated from flux ratios; intermediates above the toxic threshold are flagged. Recommendations for upregulating/downregulating specific enzymes are generated.

**Outputs:** TSV of reaction fluxes (`pathway_fluxes.tsv`); bottleneck report (`bottlenecks.tsv`); optional pathway map with flux overlay.

**How it ships:** `python run_balancer.py --pathway pathway.json`; `main.py` delegates to `src.pipeline.main()` which loads `run_balancer.py` via `importlib`.

## Usage

```bash
# Balance a metabolic pathway
python run_balancer.py --pathway pathway.json -o results/

# Synthetic demo (no real input required)
python run_balancer.py --demo -o results/

# Maximise productivity
python run_balancer.py --pathway pathway.json --objective maximise_productivity -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0
- matplotlib>=3.7.0

## Status

Planned
