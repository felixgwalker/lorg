# Ortholog Mapper

Maps orthologous genes between a query and one or more target species using reciprocal best BLAST hits or OMA, with optional synteny validation.

## Overview

Given protein FASTAs for a query species and one or more target species, this tool identifies orthologous gene pairs, classifies the relationship type (1:1, 1:N, N:1, N:N), and optionally filters to those supported by conserved synteny context.

## Approach

**Inputs:** Query species protein FASTA; one or more target species protein FASTAs.

**Core method:** All-vs-all BLAST is run between the query and each target. Reciprocal best hits (RBH) are identified: gene A in species X is an RBH of gene B in species Y if A is the best hit of B and B is the best hit of A. Hits are filtered by e-value (< 1e-10) and identity (≥ 30 %). Relationship type is inferred from the count of RBHs per query gene per target. Synteny validation (when enabled) checks that neighbouring genes within a ±10-gene window also show homology, providing confidence that the RBH reflects orthology rather than in-paralogy. OMA and Inparanoid methods use graph-based clustering instead.

**Outputs:** TSV of ortholog pairs per species pair (`orthologs.tsv`); ortholog group table (`ortholog_groups.tsv`); optional upset plot of species coverage.

**How it ships:** `python run_mapper.py query.fa --targets target1.fa target2.fa`; `main.py` delegates to `src.pipeline.main()` which loads `run_mapper.py` via `importlib`.

## Usage

```bash
# Map orthologs between query and two target species
python run_mapper.py human.fa --targets mouse.fa zebrafish.fa -o results/

# Synthetic demo (no real input required)
python run_mapper.py --demo -o results/

# Use OMA with synteny support
python run_mapper.py human.fa --targets mouse.fa --method OMA --synteny-support -o results/
```

## Dependencies

- biopython>=1.81
- pandas>=2.0.0
- numpy>=1.24.0

## Status

Planned
