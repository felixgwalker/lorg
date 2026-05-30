# Conservation Priority Ranker

Ranks multiple populations by conservation urgency using a composite genomic priority score combining inbreeding, Ne, adaptive diversity, unique alleles, and threat status.

## Overview

Given VCFs for multiple populations and optional metadata, this tool computes per-population genomic metrics, combines them into a weighted composite score, assigns each population to a conservation priority tier (critical/high/medium/low), and recommends management actions.

## Approach

**Inputs:** One or more population VCFs; optional metadata TSV (population size, IUCN/national threat status, geographic data).

**Core method:** Per-population metrics are computed: (1) inbreeding level (FROH or Fis); (2) effective population size (LD-based); (3) adaptive diversity score (Fst outlier fraction × heterozygosity); (4) unique allele fraction (alleles found only in this population, contributing irreplaceable diversity); (5) threat status (numeric mapping of IUCN categories). Metrics are normalised [0, 1] and combined with weights (inbreeding 30 %, Ne 25 %, adaptive diversity 20 %, unique alleles 15 %, threat status 10 %). Populations with composite score > 0.75 are critical; 0.5–0.75 high; 0.25–0.5 medium; < 0.25 low.

**Outputs:** TSV of ranked populations with tier and score breakdown (`conservation_priorities.tsv`); optional heatmap of metrics.

**How it ships:** `python run_ranker.py --vcfs pop1.vcf pop2.vcf pop3.vcf`; `main.py` delegates to `src.pipeline.main()` which loads `run_ranker.py` via `importlib`.

## Usage

```bash
# Rank populations by conservation priority
python run_ranker.py --vcfs pop1.vcf pop2.vcf pop3.vcf -o results/

# Synthetic demo (no real input required)
python run_ranker.py --demo -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0
- matplotlib>=3.7.0

## Status

Planned
