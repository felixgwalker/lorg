# Lineage-Specific Gene Finder

Identifies orphan (taxonomically restricted) genes in a query species by screening for proteins with no detectable homolog in outgroup proteomes.

## Overview

Given a focal species protein FASTA and a set of outgroup proteomes, this tool runs BLAST against each outgroup, retains genes with no significant hit in a minimum number of outgroups, and classifies likely origins (de novo, rapid divergence, horizontal transfer, or unknown).

## Approach

**Inputs:** Focal species protein FASTA; one or more outgroup species protein FASTAs.

**Core method:** Each focal gene is BLASTed against each outgroup proteome. A gene is classified as lineage-specific if it has no hit with e-value < 1e-3 and query coverage ≥ 30 % in at least `min_outgroup_species` outgroups. The remaining candidates are further screened: genes with no hits at all are provisional orphans; genes with hits only in distant taxa are possible horizontal gene transfers (flagged by codon usage bias); genes with weak hits (rapid divergence candidates) are flagged by dN/dS analysis. Expression evidence (if a transcriptome is provided) raises confidence in functional genes.

**Outputs:** TSV of lineage-specific genes with origin classification (`lineage_specific_genes.tsv`); optional Venn diagram of hit coverage across outgroups.

**How it ships:** `python run_finder.py query.fa --outgroups out1.fa out2.fa`; `main.py` delegates to `src.pipeline.main()` which loads `run_finder.py` via `importlib`.

## Usage

```bash
# Find lineage-specific genes
python run_finder.py focal.fa --outgroups out1.fa out2.fa out3.fa -o results/

# Synthetic demo (no real input required)
python run_finder.py --demo -o results/
```

## Dependencies

- biopython>=1.81
- pandas>=2.0.0
- numpy>=1.24.0

## Status

Planned
