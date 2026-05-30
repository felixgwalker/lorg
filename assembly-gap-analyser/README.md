# Assembly Gap Analyser

Analyses gaps (N-runs) in a genome assembly, classifying their genomic context and potential functional impact.

## Overview

Given a genome assembly FASTA, this tool locates all N-runs, classifies each gap by type (contig, scaffold, centromere, telomere), and optionally determines whether gaps interrupt annotated gene models or lie between synteny-supported gene positions.

## Approach

**Inputs:** Assembly FASTA; optional gene annotation GTF; optional reference annotation for synteny-based gap prioritisation.

**Core method:** Each sequence is scanned for runs of N or n of length ≥ `min_gap_length`. Gap context is inferred: gaps ≥ 50 kbp near sequence ends may be centromeric; sequences starting/ending with (TTAGGG)ₙ are telomeric. If a gene annotation is provided, gaps that overlap gene models are flagged as gene-disrupting. If a reference annotation is provided, synteny analysis identifies gap positions where orthologous genes are present in the reference but absent in the query — these are prioritised as gaps likely to contain missing gene content.

**Outputs:** BED of gap coordinates (`assembly_gaps.bed`); TSV of gap classification and gene impact (`gap_summary.tsv`); optional gap length distribution plot.

**How it ships:** `python run_analyser.py assembly.fa --annotation genes.gtf`; `main.py` delegates to `src.pipeline.main()` which loads `run_analyser.py` via `importlib`.

## Usage

```bash
# Analyse gaps in an assembly
python run_analyser.py assembly.fa --annotation genes.gtf -o results/

# Synthetic demo (no real input required)
python run_analyser.py --demo -o results/
```

## Dependencies

- biopython>=1.81
- pandas>=2.0.0
- numpy>=1.24.0
- matplotlib>=3.7.0

## Status

Planned
