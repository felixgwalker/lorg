# Post-Mortem Damage Simulator

Simulates post-mortem damage in ancient DNA reads by applying the Briggs model to generate C→T / G→A deamination patterns and realistic fragment length distributions.

## Overview

Given a reference genome FASTA, this tool samples reads, introduces post-mortem damage according to a configurable damage model, and outputs a FASTQ file with authentic ancient DNA characteristics for use in tool benchmarking and contamination testing.

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
