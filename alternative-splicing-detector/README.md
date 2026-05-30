# Alternative Splicing Detector

Detects differential alternative splicing events (exon skipping, intron retention, alt 5'/3' sites) between two RNA-seq conditions using PSI quantification.

## Overview

Given RNA-seq BAMs for two conditions and a gene annotation GTF, this tool quantifies percent spliced in (PSI) for each alternative splicing event, tests for significant differences using a likelihood ratio test, and reports significant events above a ΔPSI threshold.

## Approach

**Inputs:** RNA-seq BAM files for condition A and B; GTF gene annotation.

**Core method:** Splicing events are enumerated from the annotation: exon skipping (ES), intron retention (IR), alternative 5' splice site (A5SS), and alternative 3' splice site (A3SS). For each event, inclusion reads (spanning the exon/intron boundary) and exclusion reads (skipping) are counted from the BAM files. PSI = inclusion / (inclusion + exclusion) per condition. A Dirichlet-multinomial or beta-binomial test is applied to assess whether ΔPSI = PSI_A − PSI_B is significant. Events with |ΔPSI| ≥ 0.1 and FDR < 0.05 are flagged.

**Outputs:** TSV of differential splicing events (`splicing_events.tsv`); optional ψ scatter plot.

**How it ships:** `python run_detector.py --bam-a condA.bam --bam-b condB.bam --annotation genes.gtf`; `main.py` delegates to `src.pipeline.main()` which loads `run_detector.py` via `importlib`.

## Usage

```bash
# Detect differential splicing
python run_detector.py --bam-a treated.bam --bam-b control.bam --annotation genes.gtf -o results/

# Synthetic demo (no real input required)
python run_detector.py --demo -o results/

# Exon skipping and intron retention only
python run_detector.py --bam-a treated.bam --bam-b control.bam --annotation genes.gtf --event-types exon_skipping intron_retention -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0
- matplotlib>=3.7.0

## Status

Planned
