# Missense Impact Scorer

Scores the functional impact of missense variants using conservation, substitution cost, and physicochemical property changes.

## Overview

Given a VCF of annotated missense variants, this tool combines three evidence streams — evolutionary conservation (PhyloP/GERP/phastCons), BLOSUM62 amino acid substitution cost, and physicochemical property group changes — into a weighted composite score. Variants are classified on a five-tier scale aligned to ACMG/AMP terminology.

## Approach

**Inputs:** VCF with HGVS protein-level annotations (HGVS_P INFO field or CSQ/ANN from VEP/SnpEff); optional transcript-level conservation BED.

**Core method:** For each missense variant, the pipeline (1) retrieves the per-base conservation score at the affected codon position, (2) looks up the BLOSUM62 entry for the ref→alt amino acid pair, (3) determines whether the substitution crosses a physicochemical property boundary (nonpolar / polar / aromatic / charged), and (4) computes a weighted composite score (conservation 40 %, BLOSUM 30 %, physicochemical 30 %). Variants are binned: benign (< 0.3), likely benign (0.3–0.45), uncertain (0.45–0.55), likely pathogenic (0.55–0.7), pathogenic (≥ 0.7).

**Outputs:** TSV of scored variants (`scored_variants.tsv`); optional scatter plot of composite scores.

**How it ships:** `python run_scorer.py variants.vcf`; `main.py` delegates to `src.pipeline.main()` which loads `run_scorer.py` via `importlib`.

## Usage

```bash
# Score a real VCF
python run_scorer.py variants.vcf -o results/

# Synthetic demo (no real input required)
python run_scorer.py --demo -o results/

# Use GERP conservation
python run_scorer.py variants.vcf --conservation-tool GERP -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0

## Status

Planned
