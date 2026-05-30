# Repeat Element Classifier

Classifies repeat elements in a genome assembly by parsing RepeatMasker output, producing a BED annotation and repeat landscape summary.

## Overview

Given a genome assembly FASTA or a pre-existing RepeatMasker .out file, this tool classifies all annotated repeats by class (LINE, SINE, LTR, DNA transposon, satellite, simple repeat) and family, computes the repeat landscape (element count by divergence from consensus), and summarises the total repeat content of the assembly.

## Approach

**Inputs:** Assembly FASTA (for de novo RepeatMasker run) or pre-computed RepeatMasker .out file.

**Core method:** RepeatMasker .out fields (class, family, divergence, position) are parsed. Elements below the minimum length threshold or above 50 % divergence are excluded. Each element is classified into one of seven top-level classes. The repeat landscape is computed by binning elements into 5 % divergence intervals per class — young elements (low divergence) represent recent insertions; old elements (high divergence) represent ancient transposable element activity. The total repeat fraction is computed as the cumulative element length divided by the assembly length.

**Outputs:** BED of annotated repeats (`repeat_elements.bed`); TSV of class-level summary; optional repeat landscape stacked bar chart.

**How it ships:** `python run_classifier.py assembly.fa`; `main.py` delegates to `src.pipeline.main()` which loads `run_classifier.py` via `importlib`.

## Usage

```bash
# Classify repeats from RepeatMasker output
python run_classifier.py --repeat-masker-out assembly.fa.out -o results/

# Synthetic demo (no real input required)
python run_classifier.py --demo -o results/
```

## Dependencies

- biopython>=1.81
- pandas>=2.0.0
- numpy>=1.24.0
- matplotlib>=3.7.0

## Status

Planned
