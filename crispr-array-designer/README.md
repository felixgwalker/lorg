# CRISPR Array Designer

Assembles multi-spacer CRISPR arrays for Cas12a, Cas9, or Cas12b from a list of target sequences.

## Overview

Designed for applications requiring simultaneous targeting of multiple loci from a single array
construct — transcriptional regulation panels, pooled screens, or multiplex knockouts. The tool
identifies PAM sites for each target, extracts spacers, checks uniqueness by k-mer matching within
the target set, and assembles the final array sequence by interspersing spacers with the appropriate
system-specific direct repeat.

## Approach

**Inputs:** FASTA or TSV of target sequences (one sequence per target); Cas system selection.

**Core method:** For each target, PAM sites are scanned using the system PAM motif (TTTV for Cas12a,
NGG for Cas9). The adjacent spacer sequence is extracted (23 nt for Cas12a, 20 nt for Cas9) and scored
on GC content (25–75 %). Spacer uniqueness is verified by checking for shared 12-mers with all other
spacers in the set to flag potential cross-reactivity. Valid spacers are assembled into the final array
as: `[DR][spacer1][DR][spacer2]...[spacerN][DR]` using the system direct repeat.

**Outputs:** Array sequence FASTA; per-spacer summary TSV including PAM position, GC, and uniqueness
score; optional array schematic.

**Dependencies reused:** biopython for sequence I/O; pandas for spacer table.

**How it ships:** `python run_designer.py targets.fa --cas Cas12a`; delegated from
`main.py → src.pipeline.main() → run_designer.py`.

## Usage

```bash
python run_designer.py targets.fa --cas Cas12a -o results/
python run_designer.py targets.fa --cas Cas9 -o results/
python run_designer.py --demo -o results/
```

## Dependencies

- biopython>=1.81
- pandas>=2.0.0
- numpy>=1.24.0
- matplotlib>=3.7.0

## Status

Planned
