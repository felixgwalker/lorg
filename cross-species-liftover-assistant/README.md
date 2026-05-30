# Cross-Species Liftover Assistant

Maps genomic BED intervals from a source species to a target species using a UCSC chain file, with a synteny-based fallback for uncovered regions.

## Overview

Given a BED file and a UCSC-format pairwise chain file (e.g., hg38-to-mm39), this tool maps each interval to the target genome, reports status (success, partial, failed, multi-mapping), and optionally falls back to BLAST-based synteny inference for intervals not covered by the chain.

## Approach

**Inputs:** BED of intervals to lift over; UCSC chain file for the source→target species pair; optional source and target FASTAs for BLAST fallback.

**Core method:** Each BED interval is intersected with the chain alignment blocks. Intervals fully covered by a single chain block are mapped directly. Partially covered intervals are reported as partial. Intervals with no chain coverage trigger the BLAST fallback: the source interval sequence is extracted, BLASTed against the target genome, and the top synteny-consistent hit above the identity threshold is used. Multi-mapping intervals (multiple hits with similar scores) are flagged. The final output BED includes target coordinates, strand, and mapping confidence.

**Outputs:** BED of successfully mapped intervals (`lifted.bed`); BED of unmapped intervals (`unmapped.bed`); TSV summary of all intervals (`liftover_summary.tsv`).

**How it ships:** `python run_assistant.py intervals.bed --chain human_to_mouse.chain`; `main.py` delegates to `src.pipeline.main()` which loads `run_assistant.py` via `importlib`.

## Usage

```bash
# Lift BED intervals using a chain file
python run_assistant.py intervals.bed --chain hg38_to_mm39.chain -o results/

# Synthetic demo (no real input required)
python run_assistant.py --demo -o results/

# Use BLAST fallback only (no chain)
python run_assistant.py intervals.bed --source-fasta hg38.fa --target-fasta mm39.fa -o results/
```

## Dependencies

- biopython>=1.81
- pandas>=2.0.0
- numpy>=1.24.0

## Status

Planned
