# Enhancer Target Linker

Links enhancer elements to candidate target genes using the activity-by-contact (ABC) model, expression correlation, or distance scoring.

## Overview

Given a BED of enhancer elements with activity scores and a TSS annotation, this tool scores enhancer-gene pairs within a genomic window and reports the most likely target gene per enhancer, optionally incorporating Hi-C contact data (ABC model) or expression correlation across samples.

## Approach

**Inputs:** BED of enhancer elements with activity scores (e.g., from ATAC-seq or H3K27ac); BED of gene TSS positions; optional enhancer activity matrix (samples × enhancers) for correlation; optional Hi-C cool/hic matrix.

**Core method:** For each enhancer, all genes with TSS within `max_distance_bp` are candidate targets. Under the **ABC model**: score = (activity × contact_frequency) / sum(activity × contact_frequency over all nearby elements); contact frequency is from the Hi-C matrix or approximated by a power-law distance decay. Under **correlation**: Pearson r between the enhancer activity vector and each gene's expression vector across samples; genes with r ≥ threshold are linked. Under **distance**: score decays as 1/d². All methods normalise scores per enhancer; the top-scoring gene is the primary predicted target.

**Outputs:** TSV of enhancer-gene links with scores and distances (`enhancer_target_links.tsv`); optional arc plot.

**How it ships:** `python run_linker.py --enhancers enhancers.bed --gene-tss tss.bed`; `main.py` delegates to `src.pipeline.main()` which loads `run_linker.py` via `importlib`.

## Usage

```bash
# Link enhancers to target genes
python run_linker.py --enhancers enhancers.bed --gene-tss tss.bed -o results/

# Synthetic demo (no real input required)
python run_linker.py --demo -o results/

# ABC model with Hi-C
python run_linker.py --enhancers enhancers.bed --gene-tss tss.bed --hic-matrix contacts.hic --method activity_by_contact -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0
- matplotlib>=3.7.0

## Status

Planned
