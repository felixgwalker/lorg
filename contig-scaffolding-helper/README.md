# Contig Scaffolding Helper

Orders and orients assembly contigs into scaffolds using paired-end read links, Hi-C contacts, or reference-guided scaffolding.

## Overview

Given assembly contigs and linking evidence (read pairs, Hi-C, or reference alignment), this tool builds a contig graph, resolves order and orientation of each contig group, and outputs a scaffold FASTA with N-filled gap sequences between joined contigs.

## Approach

**Inputs:** Assembly contigs FASTA; BAM of paired-end or Hi-C reads (or reference alignment for guided scaffolding).

**Core method:** Read links are extracted from the BAM: for paired-end reads, pairs mapping to different contigs provide ordering evidence; for Hi-C, contact counts between contig pairs form a proximity matrix. A contig link graph is built; edges with fewer than `min_link_support` links are pruned. The graph is traversed to identify linear scaffold paths (resolving branchpoints by link weight). Gap size between adjacent contigs is estimated from the insert size distribution (paired-end) or from Hi-C contact decay. Final scaffolds are emitted as FASTA with Ns for gaps. Unplaced contigs are reported separately.

**Outputs:** Scaffold FASTA (`scaffolds.fa`); AGP describing contig-scaffold mapping; TSV of contig link evidence.

**How it ships:** `python run_helper.py contigs.fa --links reads.bam`; `main.py` delegates to `src.pipeline.main()` which loads `run_helper.py` via `importlib`.

## Usage

```bash
# Scaffold contigs with paired-end reads
python run_helper.py contigs.fa --links mapped_reads.bam -o results/

# Synthetic demo (no real input required)
python run_helper.py --demo -o results/

# Hi-C scaffolding
python run_helper.py contigs.fa --links hic_contacts.cool --evidence-type hi_c -o results/
```

## Dependencies

- biopython>=1.81
- pandas>=2.0.0
- numpy>=1.24.0

## Status

Planned
