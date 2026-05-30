# Variant Pathogenicity Aggregator

Aggregates evidence from ClinVar and in silico predictors into ACMG/AMP criteria and produces a composite five-tier pathogenicity classification.

## Overview

Given an annotated VCF (with ClinVar, CADD, REVEL, SpliceAI INFO fields), this tool maps each available evidence signal to the appropriate ACMG/AMP criterion, applies the weighted point-based scoring scheme from Richards et al. 2015, and outputs a five-tier classification (pathogenic / likely pathogenic / VUS / likely benign / benign) alongside the full evidence breakdown.

## Approach

**Inputs:** Annotated VCF with INFO fields for ClinVar significance, CADD Phred, REVEL, SpliceAI delta scores, gnomAD AF; optional supplementary TSV of additional tool scores.

**Core method:** Each criterion is evaluated against pre-defined rules: PVS1 for predicted null variants in haploinsufficient genes; PS1/PM5 for amino acid changes at known pathogenic positions; PM2 for absent/extremely rare in gnomAD; PP3 for concordant in silico predictions (CADD ≥ 25, REVEL ≥ 0.7, SpliceAI ≥ 0.5); BA1 for common variants (AF ≥ 5 %). Points are summed per the ACMG point framework; ties and conflicts are reported. ClinVar evidence at star level ≥ 1 contributes PP5/BP6.

**Outputs:** TSV of classified variants with evidence codes and points (`pathogenicity_classifications.tsv`); optional bar chart of classification distribution.

**How it ships:** `python run_aggregator.py variants.vcf`; `main.py` delegates to `src.pipeline.main()` which loads `run_aggregator.py` via `importlib`.

## Usage

```bash
# Classify variants from a real annotated VCF
python run_aggregator.py variants.vcf -o results/

# Synthetic demo (no real input required)
python run_aggregator.py --demo -o results/

# Require 2-star ClinVar evidence
python run_aggregator.py variants.vcf --clinvar-stars 2 -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0

## Status

Planned
