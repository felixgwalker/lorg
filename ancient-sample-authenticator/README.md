# Ancient Sample Authenticator

Authenticates ancient DNA samples by evaluating fragment length, terminal deamination, contamination, endogenous fraction, and coverage into a composite verdict.

## Overview

Given a BAM file of aligned ancient DNA reads, this tool integrates multiple authentication criteria into a weighted composite score and assigns a verdict: authentic, likely authentic, uncertain, or modern contamination.

## Approach

**Inputs:** BAM file of aligned ancient DNA reads; reference genome FASTA.

**Core method:** The following criteria are evaluated: (1) **Fragment length** — mean < 200 bp; ancient samples typically peak at 30–80 bp. (2) **Deamination damage** — C→T rate at 5'-most position ≥ 5 % (ancient) or ≥ 15 % (highly authentic). (3) **Contamination estimate** — from MT consensus and/or ANGSD; contamination < 5 % is required. (4) **Endogenous fraction** — fraction of reads mapping to the reference after quality filtering; ≥ 1 % is a weak pass. (5) **Coverage** — mean coverage of reference; ≥ 0.5× is required for reliable damage assessment. Each criterion is scored and weighted (deamination 30 %, contamination 25 %, fragment length 20 %, endogenous fraction 15 %, coverage 10 %). Composite score ≥ 0.75 → authentic; 0.5–0.75 → likely authentic; < 0.5 → uncertain or failed.

**Outputs:** Authentication report TSV (`authentication_report.tsv`); optional radar/spider chart.

**How it ships:** `python run_authenticator.py sample.bam --reference genome.fa`; `main.py` delegates to `src.pipeline.main()` which loads `run_authenticator.py` via `importlib`.

## Usage

```bash
# Authenticate an ancient DNA sample
python run_authenticator.py sample.bam --reference hg38.fa -o results/

# Synthetic demo (no real input required)
python run_authenticator.py --demo -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0
- matplotlib>=3.7.0

## Status

Planned
