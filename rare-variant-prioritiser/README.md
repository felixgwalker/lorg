# Rare Variant Prioritiser

Prioritises rare variants by combining allele frequency, CADD score, gene constraint, and phenotype matching into a three-tier ranking.

## Overview

Given an annotated VCF, this tool filters to rare variants (default gnomAD AF < 1 %) and scores each on a weighted composite of evidence: ultra-rarity, CADD Phred score, gene LOEUF constraint, HPO phenotype-to-gene matching, ClinVar pathogenicity, and optional gene panel membership. Variants are assigned to Tier 1 (most compelling), Tier 2, or Tier 3.

## Approach

**Inputs:** Annotated VCF (VEP or SnpEff output with gnomAD AF, CADD, ClinVar INFO fields); optional HPO term list; optional gene panel text file.

**Core method:** After AF filtering, each variant accumulates evidence codes (ULTRA_RARE, HIGH_CADD, CONSTRAINED_GENE, HPO_MATCH, CLINVAR_PATHOGENIC, IN_PANEL). Each code carries a weight; the sum is normalised to [0, 1]. High-impact consequences (stop-gained, frameshift, canonical splice) receive a +0.15 bonus. Variants above 0.75 are Tier 1; 0.50–0.75 are Tier 2; below are Tier 3.

**Outputs:** TSV of prioritised variants with tier assignments and evidence codes (`prioritised_variants.tsv`); optional rank plot.

**How it ships:** `python run_prioritiser.py variants.vcf`; `main.py` delegates to `src.pipeline.main()` which loads `run_prioritiser.py` via `importlib`.

## Usage

```bash
# Prioritise variants from a real annotated VCF
python run_prioritiser.py variants.vcf -o results/

# Synthetic demo (no real input required)
python run_prioritiser.py --demo -o results/

# Add HPO phenotype matching and gene panel
python run_prioritiser.py variants.vcf --hpo-terms HP:0001250 HP:0000256 --gene-panel panel.txt -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0

## Status

Planned
