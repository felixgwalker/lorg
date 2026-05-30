# Selection Sweep Detector

Detects positive selection sweeps using iHS, XP-EHH, the composite likelihood ratio (CLR), and Tajima's D in sliding windows.

## Overview

Given a phased multi-population VCF, this tool computes haplotype homozygosity and diversity statistics genome-wide, identifies outlier windows consistent with recent positive selection, and merges adjacent outliers into candidate sweep regions annotated with nearby genes.

## Approach

**Inputs:** Phased multi-population VCF; TSV population map; optional gene annotation BED.

**Core method:** (1) **iHS** — standardised integrated haplotype homozygosity; large |iHS| indicates extended haplotype blocks at intermediate-frequency alleles typical of incomplete sweeps. (2) **XP-EHH** — cross-population EHH ratio; high values indicate a sweep in one population relative to another. (3) **CLR** — composite likelihood ratio test of the SFS against a selective sweep model. (4) **Tajima's D** — strongly negative D in a window indicates an excess of rare variants consistent with a recent sweep. A composite Z-score is computed across tests, and windows above the outlier percentile threshold are flagged. Adjacent outlier windows are merged into sweep regions; sweep type (hard/soft/incomplete) is inferred from haplotype homozygosity profiles.

**Outputs:** TSV of per-window sweep statistics (`sweep_windows.tsv`); TSV of merged sweep regions with candidate genes (`sweep_regions.tsv`); optional Manhattan plot.

**How it ships:** `python run_detector.py variants.vcf --pop-map populations.tsv`; `main.py` delegates to `src.pipeline.main()` which loads `run_detector.py` via `importlib`.

## Usage

```bash
# Detect selection sweeps
python run_detector.py variants.vcf --pop-map populations.tsv -o results/

# Synthetic demo (no real input required)
python run_detector.py --demo -o results/

# iHS and XP-EHH only, with gene annotation
python run_detector.py variants.vcf --pop-map populations.tsv --tests iHS XP-EHH --gene-annotation genes.bed -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0
- matplotlib>=3.7.0

## Status

Planned
