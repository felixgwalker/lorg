# pegRNA Optimiser

Grid-searches PBS and RT template length combinations to find Pareto-optimal pegRNA designs.

## Overview

Given a target locus and desired edit, exhaustively enumerates all (PBS length, RT length) combinations
within configurable ranges and scores each by predicted efficiency features. Returns both the full
ranked candidate list and the Pareto-optimal front — designs where no other candidate is strictly better
on both efficiency and synthesis complexity simultaneously.

## Approach

**Inputs:** Target locus FASTA (≥200 bp); edit specification JSON; optional custom PBS/RT ranges.

**Core method:** Grid search over all PBS lengths in [8, 15] nt and RT lengths in [10, 16] nt (or
user-specified). For each combination: compute efficiency score using the same feature model as
prime-edit-efficiency-predictor (PBS GC, RT GC, MFE approximation, spacer score); compute synthesis
complexity score (function of total pegRNA length and internal repeats). Pareto-optimal front is
extracted by dominance comparison. Heatmap visualisation shows efficiency score across the PBS × RT
length grid.

**Outputs:** Ranked TSV of all candidates; TSV of Pareto-optimal front; optional heatmap PNG.

**Dependencies reused:** biopython for sequence handling; numpy for grid operations.

**How it ships:** `python run_optimiser.py target.fa --edit edit.json`; delegated from
`main.py → src.pipeline.main() → run_optimiser.py`.

## Usage

```bash
python run_optimiser.py target.fa --edit edit.json -o results/
python run_optimiser.py target.fa --edit edit.json --pbs-range 10-13 --rt-range 10-14 -o results/
python run_optimiser.py --demo -o results/
```

## Dependencies

- biopython>=1.81
- pandas>=2.0.0
- numpy>=1.24.0
- matplotlib>=3.7.0
- scipy>=1.10.0

## Status

Planned
