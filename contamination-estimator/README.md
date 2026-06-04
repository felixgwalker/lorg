# Contamination Estimator

Estimates modern human contamination in ancient DNA BAM files using mitochondrial consensus deviation, X-chromosome heterozygosity, and ANGSD-based likelihood approaches.

## Overview

Given a BAM file of ancient DNA reads and a reference genome, this tool applies multiple complementary methods to estimate the fraction of reads derived from modern human contamination, combining their estimates into a single QC verdict.

## Approach

**Baseline:** ANGSD (Korneliussen et al. 2014), schmutzi (Renaud et al. 2015), and
ContamMix (Fu et al. 2013) each provide likelihood-based contamination estimates from
different evidence streams. This tool does not re-derive those likelihoods — it calls
each as a dependency (subprocess or library).

**Novel layer:** The novel contribution is the *composite verdict layer*: combining
MT-consensus deviation, X-heterozygosity, ANGSD GLF, and schmutzi deamination-aware
estimates into a single weighted estimate with a QC decision threshold. A plain call
to any single tool leaves the user to integrate the results manually; this tool
surfaces a pass/fail verdict calibrated for de-extinction reference-sample QC
pipelines.

**Inputs:** BAM file of ancient DNA reads (mapped to reference genome); reference genome FASTA.

**Core method:** (1) **MT consensus** — reads mapping to the mitochondrial genome are piled up; non-consensus alleles inconsistent with post-mortem damage patterns (C→T, G→A) are flagged as potential contamination. (2) **X-chromosome heterozygosity** — for male samples, excess heterozygosity on the X chromosome (which should be hemizygous) indicates contamination. (3) **ANGSD GLF** — genotype likelihood-based estimation of contamination fraction using a panel of high-frequency polymorphisms. (4) **schmutzi** — deamination-aware contamination estimation. Estimates from all methods are combined (weighted mean); the combined estimate is compared to the threshold (default 3 %).

**Outputs:** TSV of per-method estimates with confidence intervals (`contamination_estimates.tsv`); QC summary; optional plot of deamination vs. contamination.

**How it ships:** `python run_estimator.py sample.bam --reference genome.fa`; `main.py` delegates to `src.pipeline.main()` which loads `run_estimator.py` via `importlib`.

## Usage

```bash
# Estimate contamination in an ancient DNA BAM
python run_estimator.py sample.bam --reference hg38.fa -o results/

# Synthetic demo (no real input required)
python run_estimator.py --demo -o results/

# MT consensus and X-chromosome methods only
python run_estimator.py sample.bam --reference hg38.fa --methods mt_consensus nuclear_X -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0
- matplotlib>=3.7.0

## Status

Planned
