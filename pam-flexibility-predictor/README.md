# PAM Flexibility Predictor

Scores PAM site availability for a panel of Cas variants at a target locus.

## Overview

Different Cas nucleases require different PAM sequences. When the canonical SpCas9 NGG PAM is absent
or poorly positioned at a target locus, this tool identifies which alternative Cas variants (SaCas9,
Cas9-NG, SpRY, AsCas12a, Cas12b, CasX) have usable PAM sites and ranks them by PAM density within
the target window.

## Approach

**Inputs:** Target locus FASTA (any length); optional list of Cas variants to score.

**Core method:** IUPAC-aware PAM site scanning using a position weight matrix for each variant's PAM
motif. Both strands are searched. PAM density is reported as sites per kb within the locus for each
variant, alongside the count of sites within 200 bp of the locus centre. SpRY and Cas9-NG near-PAMless
variants are included as fallback options with lower specificity penalties.

**Outputs:** Compatibility matrix TSV (variants × density metrics); optional bar chart.

**Dependencies reused:** biopython for FASTA I/O and IUPAC expansion; numpy for density calculation.

**How it ships:** `python run_predictor.py target.fa`; delegated from
`main.py → src.pipeline.main() → run_predictor.py`.

## Usage

```bash
python run_predictor.py target.fa -o results/
python run_predictor.py target.fa --cas-variants SpCas9 SaCas9 Cas9-NG -o results/
python run_predictor.py --demo -o results/
```

## Dependencies

- biopython>=1.81
- pandas>=2.0.0
- numpy>=1.24.0

## Status

Planned
