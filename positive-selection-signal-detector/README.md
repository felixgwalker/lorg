# Positive Selection Signal Detector

Detects genomic regions and individual genes showing evidence of positive (adaptive) selection across a multi-species alignment, identifying candidate loci likely responsible for lineage-specific adaptive traits. In de-extinction and conservation genomics, pinpointing positively selected genes in the target extinct species highlights the functionally important changes that must be faithfully recapitulated — not merely the most divergent sequences — when engineering a proxy genome.

## Inputs

- Multi-species codon-aware alignment in FASTA or MAF format covering coding sequences of interest
- A species tree topology in Newick format
- Optional: GFF3 gene annotation for the reference species
- Parameters: foreground branch label (the lineage of interest for branch-site tests), significance p-value threshold, minimum alignment length filter

## Outputs

- A per-gene dN/dS ratio table (ω values) in CSV format with likelihood ratio test statistics and adjusted p-values
- A list of genes with significant evidence of positive selection, annotated with gene name and functional description
- A visualisation of dN/dS ratios across the genome or across a gene set (PNG/SVG)
- A branch-site test results table for foreground-lineage-specific selection if a foreground branch is specified

## Method

Extracts codon-aligned sequences for each gene from the input alignment. Computes maximum-likelihood estimates of dN/dS under the M1a (neutral) and M2a (positive selection) models using a codon substitution framework. Applies a likelihood ratio test (LRT) comparing nested models to identify genes with ω > 1. Optionally applies the branch-site model to test for episodic positive selection on the specified foreground lineage. Corrects for multiple testing using the Benjamini-Hochberg FDR procedure.

## Dependencies

- `biopython` — codon alignment I/O and sequence utilities
- `scipy` — likelihood ratio test statistics and p-value computation
- `pandas` — results table construction and FDR correction
- `numpy` — numerical likelihood optimisation
- `matplotlib` — dN/dS genome-wide and gene-set plots
- `ete3` — phylogenetic tree manipulation and branch labelling
