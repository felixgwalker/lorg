# Biosafety Risk Assessor

Screens synthetic biology sequences against select agent, virulence factor, antibiotic resistance, and toxin databases to recommend biosafety containment levels.

## Overview

Given a FASTA of DNA or protein sequences, this tool performs BLAST-based homology screening against curated biosafety databases and assigns each sequence a recommended biosafety level (BSL-1 to BSL-4) with supporting flags and containment recommendations.

## Approach

**Inputs:** DNA or protein FASTA of synthetic sequences to assess.

**Core method:** Each sequence is BLASTed against: (1) **Select agent database** — CDC/USDA select agent and toxin sequences; any hit triggers a high-risk flag. (2) **VFDB** — virulence factor database; hits are classified by factor type (adhesin, toxin, immune evasion, etc.). (3) **CARD** — Comprehensive Antibiotic Resistance Database; hits flag potential antibiotic resistance transfer risk. (4) **ToxProt** — animal toxin database. Hits passing identity ≥ 50 % and coverage ≥ 50 % at e-value < 1e-5 are included. Flags are combined into a composite risk score; recommended BSL is assigned based on the highest-severity flag. Containment recommendations (physical, operational, institutional) are generated per flag. Note: this tool provides a computational screening assessment only and does not replace regulatory compliance review.

**Outputs:** TSV of per-sequence biosafety assessments (`biosafety_assessments.tsv`); flag detail TSV; optional risk summary table.

**How it ships:** `python run_assessor.py sequences.fa`; `main.py` delegates to `src.pipeline.main()` which loads `run_assessor.py` via `importlib`.

## Usage

```bash
# Screen sequences for biosafety risks
python run_assessor.py sequences.fa -o results/

# Synthetic demo (no real input required)
python run_assessor.py --demo -o results/

# Screen without antibiotic resistance check
python run_assessor.py sequences.fa --no-resistance -o results/
```

## Dependencies

- biopython>=1.81
- pandas>=2.0.0
- numpy>=1.24.0

## Status

Planned
