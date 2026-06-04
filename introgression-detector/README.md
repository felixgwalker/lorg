# Introgression Detector

Detects gene flow between populations using Patterson's D-statistic and f4-ratio in a four-population ABBA-BABA framework.

## Overview

Given a multi-population VCF and a population map, this tool computes genome-wide and sliding-window ABBA-BABA statistics to identify evidence of introgression between a putative donor population (P3) and one of two target populations (P1/P2), using an outgroup to polarise alleles.

## Approach

**Baseline:** Dsuite (Malinsky et al. 2021) and ADMIXTOOLS2 (Patterson et al. 2012;
Maier et al. 2023) implement D-statistic and f4-ratio tests. This tool does not
re-implement those statistics — it calls Dsuite as a dependency.

**Novel layer:** Standard D-statistics assume genotype-quality equivalent to modern
WGS. aDNA sequences have elevated error rates and systematic deamination; the
standard ABBA-BABA test mistakes C→T damage artefacts for derived alleles, inflating
D in the direction of the damaged sequence. The novel contribution is a
*damage-aware null model*: expected ABBA/BABA counts are modelled under the observed
damage profile (from `DamageProfile`) to adjust the null distribution, improving
power for genuinely ancient introgression while reducing false positives from damage.
Without this correction, D-stat results on low-quality aDNA are not publishable.

**Inputs:** Multi-population VCF of bi-allelic SNPs; TSV mapping sample IDs to population labels (columns: sample, population).

**Core method:** For each SNP, ABBA (ancestral in P1, derived in P2 and P3) and BABA (derived in P1 and P3, ancestral in P2) patterns are counted. D = (nABBA − nBABA) / (nABBA + nBABA); significant positive D indicates P3→P2 gene flow. Block jackknife (default 100 blocks) provides standard errors and Z-scores. The f4-ratio estimates the introgressed fraction. Sliding-window D values localise introgressed segments. Dfoil is supported for five-taxon topologies.

**Outputs:** TSV of window-level D statistics (`introgression_windows.tsv`); genome-wide summary (`introgression_summary.tsv`); optional Manhattan plot.

**How it ships:** `python run_detector.py variants.vcf --pop-map populations.tsv`; `main.py` delegates to `src.pipeline.main()` which loads `run_detector.py` via `importlib`.

## Usage

```bash
# Detect introgression across populations
python run_detector.py variants.vcf --pop-map populations.tsv --p1 EUR --p2 EAS --p3 NEA --outgroup CHI -o results/

# Synthetic demo (no real input required)
python run_detector.py --demo -o results/

# Use f4-ratio instead of D-statistic
python run_detector.py variants.vcf --pop-map populations.tsv --test f4_ratio -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0
- matplotlib>=3.7.0

## Status

Planned
