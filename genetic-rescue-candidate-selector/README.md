# Genetic Rescue Candidate Selector

Ranks donor populations for genetic rescue of an inbred recipient by balancing expected heterozygosity gain, kinship distance, ecotype compatibility, and outbreeding depression risk.

## Overview

Given genomic data for a recipient (inbred) population and a set of candidate donor populations, this tool scores each donor on a composite rescue index and recommends the best candidates for managed translocations.

## Approach

**Inputs:** VCF of the recipient population; VCFs of candidate donor populations; optional population metadata (ecotype, location coordinates).

**Core method:** For each donor–recipient pair: (1) **Kinship** is computed using the KING estimator; donors with kinship ≥ 0.25 (close relatives) are penalised; donors with very low kinship (effectively unrelated) are favoured for heterozygosity gain but assessed for outbreeding depression risk. (2) **Expected heterozygosity gain** = He_recipient_post − He_recipient_pre under a simple mixing model. (3) **Ecotype compatibility** — donor and recipient are compared on available ecotype/subspecies labels; mismatches penalise the score. (4) **Geographic proximity** — a distance decay is applied if coordinates are available. A composite rescue score is computed; donors ranked first are recommended as transplant candidates.

**Outputs:** TSV of ranked donor candidates with sub-scores (`rescue_candidates.tsv`); optional bar chart of candidate scores.

**How it ships:** `python run_selector.py --recipient recipient.vcf --donors donor1.vcf donor2.vcf`; `main.py` delegates to `src.pipeline.main()` which loads `run_selector.py` via `importlib`.

## Usage

```bash
# Rank donors for genetic rescue
python run_selector.py --recipient recipient.vcf --donors donor1.vcf donor2.vcf -o results/

# Synthetic demo (no real input required)
python run_selector.py --demo -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0
- matplotlib>=3.7.0

## Status

Planned
