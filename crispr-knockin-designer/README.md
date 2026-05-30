# CRISPR Knock-in Designer

Designs sgRNA and HDR donor template for CRISPR-mediated precise knock-in.

## Overview

Precision knock-in requires a guide that cuts close to the desired insertion point plus a homology-
directed repair (HDR) donor carrying the insert flanked by homology arms. This tool automates donor
template design including silent PAM mutation (to prevent re-cutting after successful integration)
and determines whether ssODN or dsDNA/AAV delivery is appropriate based on insert size.

## Approach

**Inputs:** Target locus FASTA (≥500 bp); insert sequence FASTA; desired insertion point (derived
from locus FASTA context or user-specified).

**Core method:** All NGG PAM sites within `--max-cut-distance` (default 30 bp) of the desired
insertion point are identified and scored by on-target features (Doench 2016 model). The
highest-scoring guide is selected. Left and right homology arms of `--arm-length` bp are extracted
from the target FASTA flanking the cut. The PAM sequence in both arms is silently mutated to prevent
re-cutting. Donor type recommendation: ssODN (≤200 bp total donor), dsDNA (201–400 bp), AAV
(>400 bp or in vivo delivery).

**Outputs:** Guide TSV; donor template FASTA with arm annotations; optional schematic.

**Dependencies reused:** biopython for FASTA I/O and reverse complement.

**How it ships:** `python run_designer.py locus.fa --insert insert.fa`; delegated from
`main.py → src.pipeline.main() → run_designer.py`.

## Usage

```bash
python run_designer.py locus.fa --insert insert.fa -o results/
python run_designer.py locus.fa --insert insert.fa --arm-length 200 -o results/
python run_designer.py --demo -o results/
```

## Dependencies

- biopython>=1.81
- pandas>=2.0.0
- numpy>=1.24.0
- matplotlib>=3.7.0

## Status

Planned
