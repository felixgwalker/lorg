# Assembly Quality Assessor

Computes N50, N90, L50, L90, GC content, gap statistics, and ambiguous base count from a genome assembly FASTA.

## Overview

Given a genome assembly FASTA, this tool calculates the full set of standard assembly quality metrics and classifies the assembly as reference-quality, chromosome-level, scaffold-level, or contig-level based on N50 thresholds.

## Approach

**Inputs:** Genome assembly FASTA (any number of sequences).

**Core method:** Sequence lengths are extracted and sorted in descending order. N50/N90 are computed as the sequence length at which 50 %/90 % of the total assembly is covered. L50/L90 are the number of sequences needed to reach that cumulative length. Gap content (runs of N/n) and ambiguous base counts are computed per sequence and summed. GC content is computed excluding Ns. Quality class is assigned: N50 ≥ 25 Mbp → reference-quality; ≥ 1 Mbp → chromosome-level; ≥ 10 kbp → scaffold-level; below → contig-level.

**Outputs:** TSV of assembly statistics (`assembly_stats.tsv`); optional sequence length distribution plot.

**How it ships:** `python run_assessor.py assembly.fa`; `main.py` delegates to `src.pipeline.main()` which loads `run_assessor.py` via `importlib`.

## Usage

```bash
# Assess assembly quality
python run_assessor.py assembly.fa -o results/

# Synthetic demo (no real input required)
python run_assessor.py --demo -o results/
```

## Dependencies

- biopython>=1.81
- pandas>=2.0.0
- numpy>=1.24.0
- matplotlib>=3.7.0

## Status

Planned
