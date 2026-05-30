# Edit Outcome Simulator

Simulates the CRISPR indel outcome distribution at a target site using an inDelphi-style approach.

## Overview

After a CRISPR cut, the resulting indel spectrum depends on the local sequence context — particularly
the microhomologies flanking the cut and the +1 nucleotide that templates 1-bp insertions. This tool
predicts the full indel distribution and frameshift rate using an inDelphi-inspired sequence model.

## Approach

**Inputs:** FASTA with ≥60 bp flanking the cut site (or guide spacer + target for auto-location); Cas
variant selection; optional number of Monte Carlo simulation draws.

**Core method:** Based on the inDelphi model (Shen et al. 2018): (1) 1-bp insertions are predicted by
templating the nucleotide immediately 3′ of the cut (position +1 on the non-template strand), assigned
a base frequency of ~50 % of all insertions. (2) Deletions of 1–30 bp are enumerated; each deletion
is scored by the microhomology score (length² × GC factor) of any homology between the sequences
flanking the deletion endpoints. (3) Scores are converted to probabilities via softmax normalisation.
(4) Frameshift probability is computed as the fraction of outcomes with (indel size mod 3) ≠ 0.

**Outputs:** Indel distribution TSV; frameshift rate JSON; optional stacked bar chart.

**How it ships:** `python run_simulator.py target.fa --guide SPACER`; delegated from
`main.py → src.pipeline.main() → run_simulator.py`.

## Usage

```bash
python run_simulator.py target.fa --guide NNNNNNNNNNNNNNNNNNN -o results/
python run_simulator.py target.fa --cas Cas12a -o results/
python run_simulator.py --demo -o results/
```

## Dependencies

- biopython>=1.81
- numpy>=1.24.0
- pandas>=2.0.0
- matplotlib>=3.7.0
- scipy>=1.10.0

## Status

Planned
