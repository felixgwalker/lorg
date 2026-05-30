# Prime Edit Design Assistant

Designs pegRNAs for prime editing from a target locus sequence and desired edit specification.

## Overview

Given a target locus FASTA and a desired nucleotide change, this tool identifies all usable PAM sites,
enumerates reverse transcriptase (RT) template and primer binding site (PBS) length combinations, and
ranks resulting pegRNA designs by predicted activity features. Optionally designs PE3 nicking guides.

## Approach

**Inputs:** Target locus FASTA (≥200 bp centred on the edit site); edit specification JSON
`{"position": int, "ref": "A", "alt": "G"}`.

**Core method:** PAM scanning (default NGG for SpCas9 PE2/PE3) on both strands of the target. For each
PAM site within range of the desired edit, a 20 nt spacer is selected and the RT template is built
encoding the edit plus downstream flanking sequence (10–16 nt). The PBS is the reverse complement of
the 3′ end of the protospacer (8–15 nt). Each (PBS length, RT length) combination is scored on:
PBS GC content (optimal 30–70 %), RT GC, RT MFE approximation, spacer on-target score, and
homopolymer penalties. PE3 nicking guides are searched in the 40–90 bp window on the non-edited strand.

**Outputs:** TSV of ranked pegRNA designs (`designs.tsv`); optional schematic PNG.

**Dependencies reused:** biopython for FASTA I/O and reverse complement; numpy/scipy for thermodynamic
approximations.

**How it ships:** `python run_assistant.py target.fa --edit edit.json`; `main.py` delegates to
`src.pipeline.main()` which loads `run_assistant.py` via `importlib`.

## Usage

```bash
# Design pegRNAs for a real target
python run_assistant.py target.fa --edit edit.json -o results/

# Synthetic demo (no real input required)
python run_assistant.py --demo -o results/

# Use SaCas9 PAM
python run_assistant.py target.fa --edit edit.json --pam NNGRRT -o results/
```

## Dependencies

- biopython>=1.81
- pandas>=2.0.0
- numpy>=1.24.0
- matplotlib>=3.7.0
- scipy>=1.10.0

## Status

Planned
