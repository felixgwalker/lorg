# Codon Optimisation Engine

Optimises protein-coding DNA sequences for expression in a target host organism by replacing each codon with the most suitable synonym from the host codon usage table.

## Overview

Given a protein FASTA and a target host, this tool generates a synonymous DNA sequence with maximised codon adaptation index (CAI), GC content within target bounds, and no specified restriction enzyme recognition sites.

## Approach

**Inputs:** Protein FASTA; host organism name (for codon usage table lookup).

**Core method:** A codon usage table is retrieved for the target host from an embedded database (based on Kazusa/GenBank coding sequences). For each amino acid in the input protein: under **CAI_maximised**, the codon with the highest weight (w = freq_i / max_freq) is selected; under **most_frequent**, the globally most common codon is used; under **harmonised**, the codon frequency pattern of the input is matched to the host; under **random_weighted**, codons are drawn by frequency. After the initial sequence is built, restriction site avoidance is applied by locally swapping synonymous codons. GC content is checked; if out of range, codons are iteratively adjusted. CAI is computed as the geometric mean of all codon weights.

**Outputs:** Optimised DNA FASTA (`optimised_sequences.fa`); TSV of CAI before/after and GC content; optional CAI distribution plot.

**How it ships:** `python run_engine.py proteins.fa --host "Homo sapiens"`; `main.py` delegates to `src.pipeline.main()` which loads `run_engine.py` via `importlib`.

## Usage

```bash
# Optimise sequences for human expression
python run_engine.py proteins.fa --host "Homo sapiens" -o results/

# Synthetic demo (no real input required)
python run_engine.py --demo -o results/

# E. coli expression, avoid EcoRI and BamHI sites
python run_engine.py proteins.fa --host "Escherichia coli" --avoid-sites EcoRI BamHI -o results/
```

## Dependencies

- biopython>=1.81
- pandas>=2.0.0
- numpy>=1.24.0
- matplotlib>=3.7.0

## Status

Planned
