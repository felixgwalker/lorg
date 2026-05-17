# Proxy Species Edit Burden Calculator

Calculates the total number and nature of genome edits required to transform a living proxy species genome toward the reconstructed genomic sequence of a target extinct species. This feasibility metric — the edit burden — is a key input for de-extinction project planning, helping teams assess whether a given proxy-target pair is tractable for CRISPR-based genome engineering within practical limits of current technology.

## Inputs

- Proxy species genome assembly in FASTA format
- Target species genome or consensus sequence in FASTA format (e.g. reconstructed from ancient DNA)
- Gene annotation files in GFF3 format for both species (optional but improves prioritisation)
- A BED file or gene list defining the regions of interest (e.g. protein-coding exons only)
- Parameters: minimum variant quality threshold, indel size classification boundaries

## Outputs

- A VCF file of all variants distinguishing proxy from target within the specified regions
- An edit burden summary table in CSV format broken down by variant class (SNV, small indel, structural variant)
- A prioritised edit list ranked by predicted phenotypic impact (using gene essentiality and conservation scores)
- A visualisation of edit density across chromosomes or target regions (PNG/SVG)

## Method

Performs pairwise genome alignment between proxy and target assemblies using an anchor-based approach. Calls variants in the aligned regions. Classifies variants by type and size. Annotates each variant with gene context, conservation score (phastCons or phyloP), and predicted functional impact (synonymous, nonsynonymous, splice, intergenic). Aggregates counts by category to produce the edit burden score. Optionally filters to a minimum edit set targeting coding regions only.

## Dependencies

- `biopython` — FASTA I/O and sequence comparison
- `pysam` — BAM/VCF handling
- `pandas` — variant table and summary statistics
- `numpy` — numerical aggregation
- `pyvcf` — VCF parsing and filtering
- `matplotlib` — edit density chromosome plot
