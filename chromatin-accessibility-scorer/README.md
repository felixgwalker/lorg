# Chromatin Accessibility Scorer

Scores chromatin accessibility at genomic peak regions from ATAC-seq BAM data, classifying regions as open, intermediate, or closed.

## Overview

Given an ATAC-seq BAM and a BED of peak regions, this tool counts Tn5 insertion events (fragment midpoints or cut sites) within each region, normalises by library size, and classifies each peak's accessibility by enrichment over a local background model.

## Approach

**Inputs:** ATAC-seq BAM file (coordinate-sorted, indexed); BED of peak regions.

**Core method:** For each peak, Tn5 cut sites (read start positions adjusted for +4/−5 bp Tn5 offset) are counted within the region. Raw counts are normalised to RPM, RPKM, or TMM. Local background is estimated from flanking 2 kb windows on each side; enrichment = peak_signal / background_signal. Peaks with enrichment ≥ 2 (adjustable) and FDR-corrected q-value < 0.05 are classified as open; 1–2 fold enrichment as intermediate; below 1 as closed. Nucleosome-free regions (NFRs) are identified from sub-nucleosomal fragments (< 200 bp).

**Outputs:** TSV of accessibility scores per peak (`accessibility_scores.tsv`); optional enrichment distribution plot.

**How it ships:** `python run_scorer.py sample.bam --peaks peaks.bed`; `main.py` delegates to `src.pipeline.main()` which loads `run_scorer.py` via `importlib`.

## Usage

```bash
# Score chromatin accessibility
python run_scorer.py sample.bam --peaks peaks.bed -o results/

# Synthetic demo (no real input required)
python run_scorer.py --demo -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0
- matplotlib>=3.7.0

## Status

Planned
