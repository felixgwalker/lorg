# Genome Rearrangement Mapper

Maps chromosomal rearrangements (inversions, translocations, fusions, fissions) between two genomes from synteny block coordinates.

## Overview

Given a synteny block TSV (output from Conserved Synteny Detector or similar), this tool analyses block orientation, chromosomal assignment, and adjacency to classify each inter-block boundary as an inversion, inter-chromosomal translocation, chromosome fusion, chromosome fission, or transposition.

## Approach

**Inputs:** TSV of synteny blocks with coordinates in both species, orientation, and chromosome labels (output of conserved-synteny-detector or MCScan).

**Core method:** Blocks are sorted by chromosome and position in species A. Consecutive blocks are compared: if adjacent blocks in A map to opposite strands of the same chromosome in B, an inversion is inferred; if they map to different chromosomes in B, a translocation is inferred; if two blocks from different chromosomes in A map adjacently on one chromosome in B, a fusion or fission is inferred depending on direction. Breakpoint coordinates are reported in both species with flanking gene annotations from a provided gene BED. Rearrangements are assigned a unique ID and size estimate.

**Outputs:** TSV of rearrangement events (`rearrangements.tsv`); TSV of breakpoints (`breakpoints.tsv`); optional dot-plot with rearrangements highlighted.

**How it ships:** `python run_mapper.py --synteny-blocks blocks.tsv`; `main.py` delegates to `src.pipeline.main()` which loads `run_mapper.py` via `importlib`.

## Usage

```bash
# Map rearrangements from synteny blocks
python run_mapper.py --synteny-blocks blocks.tsv -o results/

# Synthetic demo (no real input required)
python run_mapper.py --demo -o results/

# Detect inversions only
python run_mapper.py --synteny-blocks blocks.tsv --rearrangement-types inversion -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- matplotlib>=3.7.0

## Status

Planned
