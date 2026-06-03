# Post-Mortem Damage Simulator

> **Internal benchmark utility — not a shipped tool.**
>
> gargammel and fragSim already provide production-grade ancient-DNA read simulators.
> This module lives in `validation/` and exists solely to generate ground-truth damaged
> reads for benchmarking every PIVOT tool's damage-aware claims. It is the mechanism
> that turns "we claim damage-awareness" into a falsifiable, per-tool assertion backed
> by simulated data with known parameters.  Do not publish or expose it as a standalone
> product.

Simulates post-mortem damage in ancient DNA reads by applying the Briggs model to generate C→T / G→A deamination patterns and realistic fragment length distributions.

## Role in the validation harness

Used by `validation/` benchmark scripts to:
1. Generate synthetic ancient reads at controlled damage levels from reference genomes.
2. Feed those reads into each PIVOT tool (damage classifier, contamination estimator, etc.).
3. Compare tool outputs against known ground truth to produce per-method accuracy metrics.

See `validation/README.md` for the full simulation pipeline.

## Approach

**Inputs:** Reference genome FASTA (for read sampling); damage model parameters.

**Core method:** Fragment positions are sampled uniformly or weighted by GC content from the reference. Fragment lengths are drawn from a geometric distribution parameterised by the mean fragment length. The Briggs model is applied: at each read end, a single-stranded overhang of geometric length is generated; cytosines in the overhang are deaminated (C→T at 5' end, G→A at 3' end) with probability `deamination_rate_ss`; double-stranded interior positions are deaminated with probability `deamination_rate_ds`. Nick frequencies introduce additional 5' overhangs stochastically. The output FASTQ includes both damaged and original sequences for validation.

**Outputs:** Damaged reads FASTQ (`simulated_reads.fastq.gz`); damage summary TSV; optional deamination profile plot.

**How it ships:** `python run_simulator.py reference.fa`; `main.py` delegates to `src.pipeline.main()` which loads `run_simulator.py` via `importlib`.

## Usage

```bash
# Simulate damaged reads from a reference genome
python run_simulator.py reference.fa --n-reads 100000 -o results/

# Synthetic demo (no real input required)
python run_simulator.py --demo -o results/

# Simulate with a uniform damage model
python run_simulator.py reference.fa --model uniform --deamination-rate-ss 0.3 -o results/
```

## Dependencies

- biopython>=1.81
- pandas>=2.0.0
- numpy>=1.24.0
- matplotlib>=3.7.0

## Status

Planned
