# Safe Harbour Integration Finder

Identifies genomic safe harbour sites for stable transgene integration.

## Overview

Not all genomic loci are equally safe for stable transgene integration. Insertion near oncogenes,
enhancers, or imprinted regions can cause silencing, activation of nearby genes, or tumourigenic
transformation. This tool screens for validated safe harbour loci (AAVS1, H11, Rosa26, CCR5) and
scores additional candidate intergenic regions by oncogene distance and regulatory context.

## Approach

**Inputs:** Genome FASTA (or chromosomal subset); optional BED files for regulatory elements and
oncogenes; species selection (human/mouse).

**Core method:** First, known validated safe harbour loci for the specified species are searched by
local FASTA alignment (AAVS1/PPP1R12C at chr19, H11/HPRT1, Rosa26 for mouse, CCR5). These are
returned as tier-1 (validated) candidates. For genome-wide candidate discovery, intergenic regions
are scored by: distance from the nearest oncogene in `--oncogene-bed` (optimal >1 Mb), distance from
nearest regulatory element in `--regulatory-bed` (optimal >50 kb), estimated repeat density as a
proxy for chromatin accessibility, and expression level of flanking genes.

**Outputs:** Ranked candidate BED; scored TSV; optional genome browser–style track plot.

**How it ships:** `python run_finder.py genome.fa --species human`; delegated from
`main.py → src.pipeline.main() → run_finder.py`.

## Usage

```bash
python run_finder.py genome.fa --species human -o results/
python run_finder.py genome.fa --regulatory-bed regulatory.bed --oncogene-bed oncogenes.bed -o results/
python run_finder.py --demo -o results/
```

## Dependencies

- biopython>=1.81
- pandas>=2.0.0
- numpy>=1.24.0
- matplotlib>=3.7.0

## Status

Planned
