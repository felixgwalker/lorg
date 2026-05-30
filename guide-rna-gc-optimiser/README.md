# Guide RNA GC Optimiser

Scores and ranks guide RNAs by GC content features for efficient Cas9-mediated editing.

## Overview

GC content is one of the strongest predictors of guide RNA on-target activity. This tool scores guide
RNAs on total GC fraction, seed-region GC (the 12 nt proximal to the PAM), homopolymer run penalties,
and poly-T avoidance (which terminates U6 transcription), producing a composite GC optimality score.

## Approach

**Inputs:** FASTA or TSV of 20 nt guide spacer sequences.

**Core method:** Rule-based scoring:
- Total GC (0.30–0.70 optimal): scored as 1 − |GC − 0.50| / 0.20, capped at 1.
- Seed-region GC (positions 1–12 from PAM, 0.40–0.60 optimal): same normalisation.
- Homopolymer penalty: −0.2 per run of ≥4 identical bases.
- Poly-T penalty: −0.3 per run of ≥4 T bases (RNA Pol III terminator).
Composite score is a weighted sum (total GC 35 %, seed GC 35 %, homopolymer 15 %, poly-T 15 %).
Guides are ranked and flagged as passing/failing the GC filter.

**Outputs:** Scored/ranked TSV with per-feature subscores; optional bar chart.

**Dependencies reused:** biopython for sequence I/O; pandas/numpy for scoring.

**How it ships:** `python run_optimiser.py guides.fa`; delegated from
`main.py → src.pipeline.main() → run_optimiser.py`.

## Usage

```bash
python run_optimiser.py guides.fa -o results/
python run_optimiser.py guides.fa --gc-min 0.40 --gc-max 0.65 -o results/
python run_optimiser.py --demo -o results/
```

## Dependencies

- biopython>=1.81
- pandas>=2.0.0
- numpy>=1.24.0
- matplotlib>=3.7.0

## Status

Planned
