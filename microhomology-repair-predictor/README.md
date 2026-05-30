# Microhomology Repair Predictor

Predicts MMEJ deletion products from microhomology sequences flanking a double-strand break.

## Overview

Microhomology-mediated end joining (MMEJ) produces predictable, templated deletions when short direct
repeats flank a DSB. By enumerating all microhomologies in the flanking sequence and scoring them, this
tool predicts the most likely deletion products and their relative frequencies — enabling CRISPR
experiments to be designed to either avoid or exploit MMEJ outcomes.

## Approach

**Inputs:** FASTA with ≥60 bp on each side of the cut site; cut position (default: sequence centre).

**Core method:** All pairs of matching subsequences (microhomologies) of length ≥ `--min-mh-length`
(default 2 bp) within `--search-window` bp on each side of the cut are enumerated by exact matching.
Each microhomology is scored: MH score = length² × (1 + GC_fraction), following the scoring scheme
of inDelphi (Shen et al. 2018). The implied deletion product (removal of the sequence between the two
matching ends) is computed. Scores are converted to predicted relative frequencies via softmax
normalisation. Products are ranked by predicted frequency.

**Outputs:** Ranked MMEJ product TSV (sequence, deletion size, MH length, predicted frequency);
optional lollipop frequency plot.

**How it ships:** `python run_predictor.py cutsite.fa`; delegated from
`main.py → src.pipeline.main() → run_predictor.py`.

## Usage

```bash
python run_predictor.py cutsite.fa -o results/
python run_predictor.py cutsite.fa --cut-position 60 --min-mh-length 3 -o results/
python run_predictor.py --demo -o results/
```

## Dependencies

- biopython>=1.81
- numpy>=1.24.0
- pandas>=2.0.0
- matplotlib>=3.7.0

## Status

Planned
