# Gene Circuit Stability Estimator

Estimates the dynamic stability and robustness of a synthetic gene circuit by simulating ODE dynamics and classifying behaviour under parameter perturbations.

## Overview

Given a gene circuit description (nodes, interactions, kinetic parameters), this tool simulates the ODE system from multiple initial conditions, identifies steady states, classifies the dynamic behaviour (stable, oscillating, bistable, or unstable), and computes a robustness score.

## Approach

**Inputs:** JSON of circuit nodes (protein, basal expression, degradation rate, Hill coefficient) and edges (source, target, interaction type, K_half, max_effect).

**Core method:** The circuit is represented as a system of ODEs: dX_i/dt = production_rate_i(X) − degradation_rate_i · X_i, where production is a Hill function of all regulators of node i. The system is integrated using RK45 from N random initial conditions. Steady states are identified as points where all derivatives ≈ 0; their number and values classify behaviour: 1 steady state → stable; 2 stable states → bistable; no stable state with sustained oscillation → oscillating; parameter-sensitive divergence → unstable. Lyapunov exponents are estimated from the Jacobian at each steady state. Robustness is quantified as the fraction of ±10 % parameter perturbations that preserve the original behaviour class.

**Outputs:** TSV of steady states and Lyapunov exponents; robustness summary; optional phase portrait plot.

**How it ships:** `python run_estimator.py --circuit circuit.json`; `main.py` delegates to `src.pipeline.main()` which loads `run_estimator.py` via `importlib`.

## Usage

```bash
# Estimate circuit stability
python run_estimator.py --circuit circuit.json -o results/

# Synthetic demo (no real input required)
python run_estimator.py --demo -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0
- matplotlib>=3.7.0

## Status

Planned
