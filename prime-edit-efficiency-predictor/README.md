# Prime Edit Efficiency Predictor

Predicts the editing efficiency of pegRNA designs using DeepPrime-style sequence features.

## Overview

Takes a set of pegRNA designs (spacer, PBS, RT template) and returns a predicted efficiency score
for each. The score is based on a feature model derived from published prime editing efficiency datasets,
capturing the key determinants of successful prime editing without requiring a deep learning runtime.

## Approach

**Inputs:** TSV or JSON of pegRNA designs with columns `id`, `spacer`, `pbs`, `rt_template`;
optionally a target locus FASTA for positional context.

**Core method:** DeepPrime-style feature extraction: PBS GC content (optimal 40–60 %), RT template
length (optimal 10–14 nt), RT GC fraction, nick distance (optimal 40–70 bp), spacer melting
temperature (Tm), and an MFE approximation of the RT template secondary structure using a
nearest-neighbour thermodynamic model (Turner 2004 parameters approximated analytically). Features
are combined using a weighted linear model calibrated on published PE efficiency data from Anzalone
et al. 2019 and Chen et al. 2021. Outputs normalised efficiency scores 0–1.

**Outputs:** Scored TSV with per-feature breakdown; optional feature contribution bar chart.

**Dependencies reused:** biopython for sequence I/O; scipy for MFE approximation.

**How it ships:** `python run_predictor.py designs.tsv --target target.fa`; delegated from
`main.py → src.pipeline.main() → run_predictor.py`.

## Usage

```bash
python run_predictor.py designs.tsv --target target.fa -o results/
python run_predictor.py --demo -o results/
```

## Dependencies

- biopython>=1.81
- numpy>=1.24.0
- scipy>=1.10.0
- pandas>=2.0.0
- matplotlib>=3.7.0

## Status

Planned
