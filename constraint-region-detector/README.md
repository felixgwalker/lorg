# Constraint Region Detector

Detects whether variants fall in genomically constrained regions using gnomAD constraint metrics.

## Overview

Given a VCF of variants and a gnomAD-format constraint annotation table, this tool intersects each variant with gene-level constraint scores (LOEUF, pLI, missense Z-score) and flags those residing in highly constrained genes or sub-genic regions.

## Approach

**Inputs:** VCF of variants; gnomAD-format constraint TSV (with columns: gene, transcript, LOEUF, pLI, mis_z, oe_lof).

**Core method:** Each variant is mapped to its overlapping gene(s) via a coordinate interval tree. The gene's constraint metrics are retrieved and evaluated against configurable thresholds (default LOEUF < 0.35, pLI > 0.9, missense Z > 3.09). The primary metric (default LOEUF) determines the constraint call; all available metrics are reported for transparency. Variants in intergenic regions are reported as unconstrained unless a constrained non-coding element BED is provided.

**Outputs:** TSV of variant constraint annotations (`constraint_overlaps.tsv`); optional scatter plot of LOEUF vs pLI for evaluated genes.

**How it ships:** `python run_detector.py variants.vcf --constraint-file gnomad_constraint.tsv`; `main.py` delegates to `src.pipeline.main()` which loads `run_detector.py` via `importlib`.

## Usage

```bash
# Detect constraint for a real VCF
python run_detector.py variants.vcf --constraint-file gnomad_constraint.tsv -o results/

# Synthetic demo (no real input required)
python run_detector.py --demo -o results/

# Use pLI as the primary metric
python run_detector.py variants.vcf --constraint-file gnomad_constraint.tsv --primary-metric pLI -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0

## Status

Planned
