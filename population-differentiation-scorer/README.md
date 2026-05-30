# Population Differentiation Scorer

Scores population differentiation using Weir-Cockerham Fst, Gst, and Jost's D across all population pairs genome-wide and in sliding windows.

## Overview

Given a multi-population VCF and a population map, this tool computes pairwise differentiation statistics genome-wide and in sliding windows, identifying Fst outlier regions that are candidates for divergent selection or local adaptation.

## Approach

**Inputs:** Multi-population VCF of bi-allelic SNPs; TSV population map (sample, population).

**Core method:** Weir-Cockerham Fst is computed for each SNP from within- and between-population variance components, then averaged across windows and genome-wide. Gst (Nei's standardised measure) accounts for within-population heterozygosity and is preferable for high-diversity loci. Jost's D is calculated as the effective number of alleles unique to each population. All pairwise combinations of populations are computed. Window-level Fst values above the outlier percentile threshold are flagged as differentiation outliers. Bootstrap confidence intervals are provided for genome-wide estimates.

**Outputs:** TSV of pairwise genome-wide statistics (`pairwise_differentiation.tsv`); TSV of per-window Fst (`window_fst.tsv`); optional Manhattan plot of window Fst.

**How it ships:** `python run_scorer.py variants.vcf --pop-map populations.tsv`; `main.py` delegates to `src.pipeline.main()` which loads `run_scorer.py` via `importlib`.

## Usage

```bash
# Score differentiation between all population pairs
python run_scorer.py variants.vcf --pop-map populations.tsv -o results/

# Synthetic demo (no real input required)
python run_scorer.py --demo -o results/

# Use 100 kb windows and flag top 1 % as outliers
python run_scorer.py variants.vcf --pop-map populations.tsv --window-size 100000 --outlier-percentile 99 -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0
- matplotlib>=3.7.0

## Status

Planned
