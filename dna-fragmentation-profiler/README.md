# DNA Fragmentation Profiler

Profiles DNA fragmentation length distribution and post-mortem deamination damage patterns in ancient DNA BAM files.

## Overview

Given a BAM file of aligned reads, this tool generates the fragment length histogram and 5'/3' terminal C→T / G→A substitution frequency plots characteristic of ancient DNA, classifying the sample's fragmentation pattern as ancient, modern-like, or degraded.

## Approach

**Inputs:** BAM file of aligned reads (indexed).

**Core method:** Fragment lengths are computed from properly paired reads or inferred from single-end read lengths. Terminal base mismatches are tabulated per position from each read end, computing C→T frequency at 5' ends and G→A frequency at 3' ends for the first `context_bases` positions. The ancient DNA signature is characterised by elevated C→T at position 1 of the 5' end (typically ≥ 15 % for authentically ancient samples). The fragmentation pattern is classified as: ancient (mean length < 80 bp and C→T rate ≥ 15 %), modern-like (mean length > 150 bp and low damage), or degraded (short but low damage — environmental rather than PMD).

**Outputs:** Fragment length histogram TSV (`fragment_lengths.tsv`); deamination profile TSV (`deamination_profile.tsv`); combined damage plot PNG.

**How it ships:** `python run_profiler.py sample.bam`; `main.py` delegates to `src.pipeline.main()` which loads `run_profiler.py` via `importlib`.

## Usage

```bash
# Profile an ancient DNA BAM
python run_profiler.py sample.bam -o results/

# Synthetic demo (no real input required)
python run_profiler.py --demo -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0
- matplotlib>=3.7.0

## Status

Planned
