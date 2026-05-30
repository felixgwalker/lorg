# Off-Target Cluster Detector

Detects genomic hotspots of CRISPR off-target activity using sliding-window density and clustering.

## Overview

While individual off-target sites are routinely reported, their genomic distribution is rarely
analysed. Clustering of off-target sites in gene-dense or regulatory regions represents a greater
safety risk than scattered sites. This tool identifies and annotates off-target hotspot clusters.

## Approach

**Inputs:** BED or TSV of off-target sites with columns: chromosome, position, strand, CFD score
(e.g. from Cas-OFFinder + CFD scoring, CRISPOR, or CIRCLE-seq).

**Core method:** (1) Sliding window density scan across each chromosome in `--window-size` bp steps,
counting off-target sites per window. Windows above a Z-score threshold (default 2.0) are candidate
hotspot regions. (2) Candidate sites within `--window-size` of each other are merged by single-linkage
clustering (DBSCAN-style with distance threshold = 50 kb). Clusters with ≥ `--min-cluster-size` sites
are reported. (3) Clusters are annotated if an annotation BED is provided (regulatory element overlap,
genic overlap, repeat-rich regions).

**Outputs:** Off-target cluster BED; cluster summary TSV; optional Manhattan-style plot of off-target
density across all chromosomes.

**How it ships:** `python run_detector.py offtargets.bed`; delegated from
`main.py → src.pipeline.main() → run_detector.py`.

## Usage

```bash
python run_detector.py offtargets.bed -o results/
python run_detector.py offtargets.bed --window-size 50000 --min-cluster-size 5 -o results/
python run_detector.py --demo -o results/
```

## Dependencies

- numpy>=1.24.0
- pandas>=2.0.0
- scipy>=1.10.0
- matplotlib>=3.7.0
- biopython>=1.81

## Status

Planned
