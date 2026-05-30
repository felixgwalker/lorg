# Genome Edit Feasibility Scorer

Scores the overall feasibility of a genome editing project from a composite of locus and context factors.

## Overview

Before committing to a genome editing experiment, it is useful to have a rapid feasibility assessment
covering PAM availability, target accessibility, and practical constraints. This tool reads a project
specification and returns a composite feasibility score (0–1) with a per-component breakdown and
band classification (unfeasible / challenging / feasible / highly feasible).

## Approach

**Inputs:** Project specification JSON with keys: `locus_fasta` (path), `edit_type` (knockout /
base-edit / prime-edit), `cell_type` (string), `delivery_method` (string); optional
`chromatin_score` (0–1, from ATAC-seq) and `essentiality_score` (0–1, from DepMap).

**Core method:** Five weighted components:
1. **PAM density** (20 %): number of NGG PAM sites within 200 bp of the target region / 20 (capped at 1).
2. **GC content** (15 %): protospacer GC score — Gaussian centred at 50 %, σ = 15 %.
3. **Chromatin accessibility** (25 %): user-supplied ATAC score or repeat-density proxy.
4. **Essentiality risk** (20 %): inverted DepMap essentiality score (high essentiality → low feasibility for loss-of-function edits).
5. **Delivery suitability** (20 %): lookup table of cell type × delivery method compatibility.

**Outputs:** Feasibility report JSON; component breakdown TSV; optional radar chart.

**How it ships:** `python run_scorer.py project_spec.json`; delegated from
`main.py → src.pipeline.main() → run_scorer.py`.

## Usage

```bash
python run_scorer.py project_spec.json -o results/
python run_scorer.py --demo -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- biopython>=1.81
- matplotlib>=3.7.0

## Status

Planned
