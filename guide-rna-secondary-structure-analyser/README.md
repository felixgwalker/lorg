# Guide RNA Secondary Structure Analyser

Analyses secondary structure of sgRNA spacer+scaffold sequences to flag guides with poor accessibility.

## Overview

Guide RNA activity depends not only on target accessibility but also on the guide RNA's own secondary
structure. Misfolding within the spacer, or base-pairing between the spacer and the scaffold sequence,
reduces the fraction of active guide molecules. This tool predicts MFE secondary structure, reports
seed-region accessibility, and flags guide-scaffold duplex formation.

## Approach

**Inputs:** FASTA or TSV of 20 nt guide spacer sequences; sgRNA scaffold identifier.

**Core method:** Each spacer is concatenated with the full sgRNA scaffold sequence (76 nt for SpCas9)
forming the full guide RNA. RNA secondary structure is predicted using a Zuker nearest-neighbour MFE
approximation (Turner 2004 parameters implemented analytically without an external ViennaRNA call).
The seed region (positions 1–12 from PAM, i.e. the 3′ 12 nt of the spacer) accessibility is computed
as the fraction of seed bases predicted to be unpaired. A guide-scaffold duplex flag is set if any
spacer subsequence of ≥4 nt is complementary to the scaffold loop regions.

**Outputs:** Structure TSV with MFE, dot-bracket string, seed accessibility, duplex flag; optional
dot-bracket schematic per guide.

**How it ships:** `python run_analyser.py guides.fa`; delegated from
`main.py → src.pipeline.main() → run_analyser.py`.

## Usage

```bash
python run_analyser.py guides.fa -o results/
python run_analyser.py guides.fa --scaffold SaCas9 -o results/
python run_analyser.py --demo -o results/
```

## Dependencies

- biopython>=1.81
- numpy>=1.24.0
- pandas>=2.0.0
- matplotlib>=3.7.0
- scipy>=1.10.0

## Status

Planned
