# Compound Heterozygosity Detector

Detects compound heterozygous variant pairs in the same gene from phased genotypes or trio data.

## Overview

Given a VCF (phased or unphased) and an optional family PED file, this tool identifies all pairs of heterozygous variants in the same gene where the two alleles are on opposite haplotypes, consistent with a recessive inheritance pattern that would not be apparent from homozygosity alone.

## Approach

**Inputs:** VCF of variants (ideally annotated with gnomAD AF and consequence); optional PLINK PED file for trio-based phasing.

**Core method:** Heterozygous variants per gene per sample are collected after AF filtering (default < 1 %). If the VCF is phased (PS/HP tags), haplotype blocks are used directly. With a trio PED, phase-by-transmission is applied: a variant inherited from each parent is inferred to reside on opposite haplotypes. For unphased singletons, all pairs of heterozygous variants in the same gene are reported as "possible" comp-hets. Confidence is graded: confirmed (phased), likely (trio-consistent), possible (unphased). Only high-impact consequences (stop-gained, frameshift, missense, canonical splice) are included by default.

**Outputs:** TSV of comp-het pairs with gene, variant coordinates, consequence, AF, phase status, and confidence (`comp_het_pairs.tsv`); optional gene-level summary.

**How it ships:** `python run_detector.py variants.vcf --ped family.ped`; `main.py` delegates to `src.pipeline.main()` which loads `run_detector.py` via `importlib`.

## Usage

```bash
# Detect comp-hets in a trio
python run_detector.py variants.vcf --ped family.ped -o results/

# Synthetic demo (no real input required)
python run_detector.py --demo -o results/

# Singleton (unphased) — reports all possible pairs
python run_detector.py variants.vcf --sample-id PROBAND_01 -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0

## Status

Planned
