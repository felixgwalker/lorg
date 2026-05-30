# Multiplex Edit Planner

Plans and validates a multiplex CRISPR editing strategy across multiple simultaneous targets.

## Overview

Introducing multiple edits simultaneously raises risks of guide cross-reactivity, overlapping cut
windows, and chromosomal translocations. This tool reads an edit manifest, checks pairwise guide
compatibility by CFD score, flags guides with overlapping cut windows, identifies translocation-risk
pairs (same chromosome, <50 Mb apart), and batches guides into compatible delivery groups.

## Approach

**Inputs:** Edit manifest JSON or TSV with columns: `target_id`, `guide_spacer`, `chromosome`,
`position`, `edit_type`.

**Core method:** For all guide pairs, CFD cross-reactivity score is computed (threshold 0.10 for
caution, 0.25 for incompatible). Pairs with cut windows within 100 bp on the same chromosome are
flagged for window overlap. Pairs on the same chromosome within 50 Mb are flagged for translocation
risk (DSB proximity drives translocation events). A greedy batching algorithm assigns guides to
delivery batches of ≤ `--max-guides-per-batch` such that no batch contains incompatible pairs.

**Outputs:** Compatibility matrix TSV; batched delivery plan JSON; compatibility heatmap.

**Dependencies reused:** pandas for matrix operations; numpy for CFD computation.

**How it ships:** `python run_planner.py manifest.json`; delegated from
`main.py → src.pipeline.main() → run_planner.py`.

## Usage

```bash
python run_planner.py manifest.json -o results/
python run_planner.py manifest.json --max-guides-per-batch 6 -o results/
python run_planner.py --demo -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- biopython>=1.81
- matplotlib>=3.7.0

## Status

Planned
