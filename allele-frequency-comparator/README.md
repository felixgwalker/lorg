# Allele Frequency Comparator

Compares allele frequencies across gnomAD populations to identify population-specific variants and quantify differentiation.

## Overview

Given a VCF and a gnomAD-format population allele frequency table, this tool retrieves per-population AC/AN/AF for each variant, computes pairwise fold changes between all population pairs, tests for significant frequency differences using Fisher's exact test (with Bonferroni correction), and estimates Fst as a summary of differentiation.

## Approach

**Inputs:** VCF or tab-separated variant list (chrom, pos, ref, alt); gnomAD-format population AF TSV with columns per population (e.g. `afr_ac`, `afr_an`).

**Core method:** Variants are intersected with the AF table. For each variant, a fold change matrix across all requested population pairs is computed (max/min AF ratio). Variants with fold change ≥ threshold and Fisher test p < 0.05 (after correction) are flagged as population-specific or population-enriched. Weir–Cockerham Fst is estimated from AC/AN values across populations as a genome-wide differentiation metric. All results are reported with confidence intervals where feasible.

**Outputs:** TSV of population-stratified AFs and comparison statistics (`af_comparisons.tsv`); optional barplot of AFs per variant.

**How it ships:** `python run_comparator.py variants.vcf --af-table gnomad_af.tsv`; `main.py` delegates to `src.pipeline.main()` which loads `run_comparator.py` via `importlib`.

## Usage

```bash
# Compare allele frequencies across all gnomAD populations
python run_comparator.py variants.vcf --af-table gnomad_af.tsv -o results/

# Synthetic demo (no real input required)
python run_comparator.py --demo -o results/

# Compare only AFR and EAS populations
python run_comparator.py variants.vcf --af-table gnomad_af.tsv --populations afr eas -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0
- matplotlib>=3.7.0

## Status

Planned
