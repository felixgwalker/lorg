# Promoter Strength Estimator

Estimates promoter strength from sequence features (TATA box, Inr, CpG island, TFBS density) and optional H3K4me3 ChIP-seq signal.

## Overview

Given a gene annotation and reference genome, this tool scores each gene's promoter region across multiple sequence and epigenomic features and produces a composite strength estimate classified as strong, moderate, weak, or silent.

## Approach

**Inputs:** GTF or BED of gene TSS positions; reference genome FASTA; optional H3K4me3 bigWig.

**Core method:** For each gene, the promoter sequence (default 2 kb upstream to 200 bp downstream of TSS) is extracted. Features scored: (1) TATA box PWM score at the −30 to −25 window; (2) Initiator element (Inr) PWM score at the TSS; (3) CpG island overlap (CpG O/E ratio ≥ 0.6 within the promoter); (4) GC content (≥ 55 % correlates with constitutive expression); (5) TFBS density (total PWM hits per kb above IC threshold); (6) H3K4me3 mean signal in the −500 to +500 window if provided. Features are normalised and combined into a composite score (0–1). Score ≥ 0.7 → strong; 0.4–0.7 → moderate; 0.2–0.4 → weak; < 0.2 → silent.

**Outputs:** TSV of per-gene promoter strength scores (`promoter_strengths.tsv`); optional score distribution plot.

**How it ships:** `python run_estimator.py --annotation genes.gtf --fasta genome.fa`; `main.py` delegates to `src.pipeline.main()` which loads `run_estimator.py` via `importlib`.

## Usage

```bash
# Estimate promoter strength
python run_estimator.py --annotation genes.gtf --fasta genome.fa -o results/

# Synthetic demo (no real input required)
python run_estimator.py --demo -o results/

# Include H3K4me3 signal
python run_estimator.py --annotation genes.gtf --fasta genome.fa --chip-h3k4me3 h3k4me3.bw -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0
- matplotlib>=3.7.0

## Status

Planned
