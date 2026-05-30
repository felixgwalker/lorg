# Base Edit Outcome Predictor

Predicts per-base editing probabilities for CBE and ABE base editors using position- and context-dependent efficiency models.

## Overview

Base editors introduce targeted C→T (CBE) or A→G (ABE) transitions without DSBs. Predicting which
bases within the editing window will be converted — and at what frequency — is critical for avoiding
bystander edits. This tool implements a BE-Hive-inspired model scoring each editable base by its
window position and trinucleotide context.

## Approach

**Inputs:** TSV with columns `id`, `spacer` (20 nt), `target_sequence` (30 nt context); base editor
type selection (CBE3, BE4max, ABE8e, ABEmax).

**Core method:** Editing window positions (default 4–8 counting from spacer 5′ end) are identified.
For CBE editors, all C bases in the window are candidates; for ABE, all A bases. Each candidate site
is scored by: (1) position efficiency lookup (position 5–6 most active, 4 and 8 least), (2)
trinucleotide context weight (e.g. TC motif strongly preferred for CBE), (3) bystander probability if
multiple editable bases are present. Indel frequency is looked up from editor-specific tables.

**Outputs:** Per-base editing probability TSV; top product prediction; optional per-target heatmap.

**Dependencies reused:** biopython for sequence manipulation; numpy for position scoring.

**How it ships:** `python run_predictor.py targets.tsv --editor ABE8e`; delegated from
`main.py → src.pipeline.main() → run_predictor.py`.

## Usage

```bash
python run_predictor.py targets.tsv --editor CBE3 -o results/
python run_predictor.py targets.tsv --editor ABE8e -o results/
python run_predictor.py --demo -o results/
```

## Dependencies

- biopython>=1.81
- pandas>=2.0.0
- numpy>=1.24.0
- matplotlib>=3.7.0

## Status

Planned
