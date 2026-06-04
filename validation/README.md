# validation/

Simulation-based test harness for the de-extinction toolkit. Turns "we claim
damage-awareness" into per-tool accuracy metrics backed by simulated ground truth.

## Structure

```
validation/
  post-mortem-damage-simulator/   # generates damaged reads with known parameters
  bench_adna_tools.py             # benchmark aDNA damage/authentication tools
  bench_popgen_tools.py           # benchmark popgen tools (msprime/SLiM)
  conftest.py                     # shared pytest fixtures
  requirements.txt
  README.md  (this file)
```

## Workflow

### aDNA tool benchmarks (`bench_adna_tools.py`)

1. `post-mortem-damage-simulator` generates synthetic reads at parameterised
   damage levels (Briggs model, known C→T rates, known fragment lengths).
2. Each PIVOT aDNA tool (ancient-dna-damage-classifier, contamination-estimator,
   palaeogenomic-coverage-assessor, …) is run on the simulated reads.
3. Tool outputs are compared to ground truth. Metrics: per-position C→T/G→A
   recovery error, authenticity classification accuracy, false-positive
   contamination rate.

### Population genetics benchmarks (`bench_popgen_tools.py`)

1. `msprime` simulates coalescent trees with known demographic history
   (Ne, bottleneck times, admixture proportions).
2. `SLiM` simulates forward-time selection on proxy-species models.
3. PIVOT popgen tools (roh-interpreter, effective-population-size-estimator,
   introgression-detector, …) are run on simulated VCFs.
4. Outputs compared to simulation truth. Metrics: Ne estimation RMSE, ROH
   recall/precision, D-stat power at known admixture proportions.

## Dependencies

```
msprime>=1.3
tskit>=0.5
SLiM  (external binary, ≥4.0)
pytest>=8.0
biopython>=1.81
numpy, pandas, scipy, matplotlib
```
(Install via `pip install -r validation/requirements.txt`; SLiM via separate
binary install.)

## Running

```bash
pytest validation/ -v                    # all benchmarks
pytest validation/bench_adna_tools.py    # aDNA tools only
pytest validation/bench_popgen_tools.py  # popgen tools only
```

## Non-negotiable for publication

Per-tool benchmark results from this harness are the evidence base for every
PIVOT tool's damage-aware claims. Without them, "damage-aware" is a marketing
word. With them, it is a falsifiable claim with a reported accuracy number.
