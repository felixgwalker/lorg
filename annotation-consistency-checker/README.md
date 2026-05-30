# Annotation Consistency Checker

Validates a GTF gene annotation for internal consistency: coordinate hierarchy, duplicate IDs, overlapping features, strand mismatches, and chromosome name agreement.

## Overview

Given a GTF annotation and optionally a reference FASTA, this tool performs comprehensive structural and semantic validation, reporting all detected issues and a pass/fail verdict.

## Approach

**Inputs:** GTF annotation file; optional reference genome FASTA (for chromosome name validation).

**Core method:** The GTF is parsed line by line. Checks performed: (1) **Duplicate IDs** — gene_id and transcript_id uniqueness; (2) **Coordinate hierarchy** — exon coordinates fall within their transcript, transcript within gene; (3) **Strand consistency** — all features of a gene share the same strand; (4) **Overlapping features** — exons within the same transcript do not overlap; CDSs of different genes on the same strand do not overlap (configurable); (5) **Missing CDS** — protein-coding transcripts contain at least one CDS feature; (6) **Coordinate validity** — start ≤ end and start ≥ 1; (7) **Chromosome names** — all chromosomes in the GTF are present in the FASTA (if provided).

**Outputs:** TSV of detected issues with coordinates and issue type (`annotation_issues.tsv`); summary statistics; optional issue-type bar chart.

**How it ships:** `python run_checker.py --annotation genes.gtf`; `main.py` delegates to `src.pipeline.main()` which loads `run_checker.py` via `importlib`.

## Usage

```bash
# Check a GTF annotation
python run_checker.py --annotation genes.gtf --fasta genome.fa -o results/

# Synthetic demo (no real input required)
python run_checker.py --demo -o results/

# Strict mode — non-zero exit on any issue
python run_checker.py --annotation genes.gtf --strict -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- matplotlib>=3.7.0

## Status

Planned
