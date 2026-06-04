# Adaptive Diversity Scorer

Scores adaptive genetic diversity in a conservation population by identifying Fst-outlier loci, testing environment–genotype associations, and computing adaptive heterozygosity relative to a neutral baseline.

## Overview

Given a population VCF, this tool identifies putatively adaptive loci and quantifies the population's adaptive diversity, classifying its adaptive capacity as high, moderate, low, or critically low.

## Approach

**Baseline:** scikit-allel (Miles et al.) computes Fst, Hudson's estimator, and
genotype–environment association statistics. This tool does not re-implement those
statistics — it calls scikit-allel for all Fst and GEA computations.

**Novel layer:** The novel contribution is the *adaptive-to-neutral classification
layer and composite score*. scikit-allel returns raw Fst values; this tool
classifies loci as adaptive or neutral using the Fst outlier distribution, combines
three independent evidence streams (Fst outlier, GEA, prior-knowledge immune/MHC
regions), and outputs a single composite adaptive capacity verdict calibrated for
conservation decision-making. The 0–1 score and the tiered verdict (critically low /
low / moderate / high) are the publishable output, not the Fst values themselves.

**Inputs:** VCF of the population; optional environmental variable TSV (for genotype–environment association); optional pre-identified adaptive loci BED.

**Core method:** Adaptive loci are identified by (1) Fst outlier detection (top 1 % of Fst distribution relative to reference populations if available, or Z-score-based); (2) genotype–environment associations (Spearman r ≥ 0.3 between allele frequency and each environmental variable); (3) prior knowledge (immune gene regions, MHC). Adaptive heterozygosity is computed at the identified loci; neutral heterozygosity is computed genome-wide excluding adaptive loci. The adaptive-to-neutral ratio scores the population's evolutionary potential. Score ≥ 0.8 → high; 0.5–0.8 → moderate; 0.3–0.5 → low; < 0.3 → critically low.

**Outputs:** TSV of adaptive loci (`adaptive_loci.tsv`); adaptive diversity summary (`adaptive_diversity_score.tsv`); optional scatter plot.

**How it ships:** `python run_scorer.py variants.vcf`; `main.py` delegates to `src.pipeline.main()` which loads `run_scorer.py` via `importlib`.

## Usage

```bash
# Score adaptive diversity
python run_scorer.py variants.vcf --environment env.tsv -o results/

# Synthetic demo (no real input required)
python run_scorer.py --demo -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0
- matplotlib>=3.7.0

## Status

Planned
