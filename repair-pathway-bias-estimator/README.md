# Repair Pathway Bias Estimator

Estimates NHEJ, MMEJ, and HDR repair pathway probabilities at a CRISPR cut site.

## Overview

The dominant repair pathway at a DSB determines whether the editing outcome is a random indel (NHEJ),
a templated deletion (MMEJ), or a precise edit (HDR). This tool predicts the relative probabilities of
each pathway by combining microhomology analysis with cell-type-specific HDR bias estimates.

## Approach

**Inputs:** FASTA with ≥50 bp flanking sequence on each side of the cut site; cell type string.

**Core method:** MMEJ probability is estimated by enumerating all microhomologies of ≥2 bp within
±30 bp of the cut. Each microhomology is scored by MH score = length² × (1 + GC_fraction). The sum
of MH scores is normalised to an MMEJ probability using a sigmoid calibration. HDR probability is
looked up from a cell-type table (iPSC: 0.30; HEK293: 0.25; primary T cells: 0.10; neurons: 0.05)
reflecting the known S/G2-phase fraction of each cell type. NHEJ probability is the remainder
(1 − MMEJ − HDR). The most likely MMEJ deletion products are returned ranked by predicted frequency.

**Outputs:** Pathway probability JSON; MMEJ product TSV; optional pie chart.

**How it ships:** `python run_estimator.py cutsite.fa --cell-type iPSC`; delegated from
`main.py → src.pipeline.main() → run_estimator.py`.

## Usage

```bash
python run_estimator.py cutsite.fa --cell-type iPSC -o results/
python run_estimator.py cutsite.fa --cell-type HEK293 --cut-position 60 -o results/
python run_estimator.py --demo -o results/
```

## Dependencies

- biopython>=1.81
- numpy>=1.24.0
- pandas>=2.0.0
- matplotlib>=3.7.0

## Status

Planned
