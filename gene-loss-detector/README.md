# Gene Loss Detector

Detects gene losses on specific lineages by identifying genes present in most species but absent or pseudogenised in others, mapped onto a species phylogeny.

## Overview

Given an ortholog presence/absence table and a Newick species tree, this tool applies parsimony to infer the branch on which each gene loss occurred, classifies the loss type (deletion, pseudogenisation, truncation, exon loss), and grades confidence based on the number of species supporting the ancestral presence.

## Approach

**Inputs:** TSV of ortholog presence/absence per species (rows = genes, columns = species, values = gene ID or absent); optional Newick species tree; optional protein FASTAs for pseudogenisation scanning.

**Core method:** For each gene, the species coverage vector is mapped onto the phylogeny. Dollo parsimony is applied: a gene is assumed to have been present in the ancestor if ≥ 3 species in independent clades have it. Loss is inferred on the most parsimonious branch. If protein FASTAs are provided, the query genome is searched for remnants (TBLASTN); frameshifts and premature stop codons in the remnant confirm pseudogenisation. Complete absence of any BLAST hit (< 30 % coverage of any exon) is classified as deletion.

**Outputs:** TSV of gene losses with branch, type, and confidence (`gene_losses.tsv`); optional phylogeny with losses mapped.

**How it ships:** `python run_detector.py --ortholog-table orthologs.tsv --phylogeny tree.nwk`; `main.py` delegates to `src.pipeline.main()` which loads `run_detector.py` via `importlib`.

## Usage

```bash
# Detect gene losses across a clade
python run_detector.py --ortholog-table orthologs.tsv --phylogeny tree.nwk -o results/

# Synthetic demo (no real input required)
python run_detector.py --demo -o results/
```

## Dependencies

- biopython>=1.81
- pandas>=2.0.0
- numpy>=1.24.0

## Status

Planned
