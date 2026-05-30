# Paralog Cluster Builder

Builds clusters of paralogous genes from a single-species proteome using all-vs-self BLAST and Markov Cluster Algorithm (MCL).

## Overview

Given a protein FASTA for a single species, this tool runs all-vs-self BLAST, filters paralog pairs by identity and coverage, clusters them with MCL, and annotates each cluster with its likely duplication mode (tandem, segmental, dispersed, or retroposed).

## Approach

**Inputs:** Single-species protein FASTA.

**Core method:** BLAST all-vs-self generates pairwise bit scores for all protein pairs. Self-hits are removed; hits below the identity and coverage thresholds are filtered. Remaining pairs are used to build a similarity graph, which MCL clusters using random walks with controllable inflation. Cluster duplication type is inferred from genomic context: tandem if members lie within 100 kb on the same chromosome, segmental if within a known segmental duplication block (≥ 10 kb, ≥ 90 % identity), retroposed if intron-less, dispersed otherwise.

**Outputs:** TSV of paralog clusters (`paralog_clusters.tsv`); TSV of pairwise paralog pairs (`paralog_pairs.tsv`); optional cluster size distribution plot.

**How it ships:** `python run_builder.py proteome.fa`; `main.py` delegates to `src.pipeline.main()` which loads `run_builder.py` via `importlib`.

## Usage

```bash
# Build paralog clusters from a proteome
python run_builder.py proteome.fa -o results/

# Synthetic demo (no real input required)
python run_builder.py --demo -o results/

# Higher MCL inflation for finer clusters
python run_builder.py proteome.fa --inflation 3.0 -o results/
```

## Dependencies

- biopython>=1.81
- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0

## Status

Planned
