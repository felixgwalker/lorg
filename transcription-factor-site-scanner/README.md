# Transcription Factor Site Scanner

Scans FASTA sequences against a JASPAR/MEME PWM database for transcription factor binding sites, reporting hits with p-values against a background model.

## Overview

Given FASTA sequences and a PWM database, this tool scores every position against each PWM and reports all hits exceeding the p-value threshold, with positions in BED format and optional sequence context.

## Approach

**Inputs:** FASTA of sequences to scan; JASPAR or MEME-format PWM database.

**Core method:** Log-odds scoring matrices are precomputed from each PWM with pseudocount = 0.1 against a uniform background. The score distribution is computed analytically (for small PWMs) or by simulation (for large ones) to derive p-values. For each sequence and each PWM, all positions where the log-odds score exceeds the threshold corresponding to p < cutoff are reported. Both strands are scanned by default. FDR correction (Benjamini-Hochberg) is applied across all sequence × PWM combinations.

**Outputs:** BED of TFBS hits (`tfbs_hits.bed`); TSV with scores and p-values (`tfbs_hits.tsv`); optional sequence logo summary.

**How it ships:** `python run_scanner.py sequences.fa --pwm-db jaspar.meme`; `main.py` delegates to `src.pipeline.main()` which loads `run_scanner.py` via `importlib`.

## Usage

```bash
# Scan sequences for TFBS
python run_scanner.py sequences.fa --pwm-db jaspar.meme -o results/

# Synthetic demo (no real input required)
python run_scanner.py --demo -o results/

# Stricter p-value cutoff, forward strand only
python run_scanner.py sequences.fa --pwm-db jaspar.meme --p-value-threshold 1e-6 --single-strand -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0
- matplotlib>=3.7.0

## Status

Planned
