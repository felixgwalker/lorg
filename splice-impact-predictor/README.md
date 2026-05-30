# Splice Impact Predictor

Predicts the impact of variants on splicing by scoring donor and acceptor sites against position weight matrices.

## Overview

Given a VCF of variants, a reference genome FASTA, and a gene annotation GTF, this tool evaluates each variant's proximity to known splice sites, scores both the reference and alternate allele sequences with MaxEntScan-style PWMs, and reports the delta score alongside an effect classification (disruption, creation, weakening, or neutral).

## Approach

**Inputs:** VCF of variants; reference genome FASTA (indexed); gene annotation GTF.

**Core method:** For each variant, the pipeline identifies all donor and acceptor sites within the search window. It extracts the reference sequence context (donor: 3 nt exonic + 6 nt intronic; acceptor: 20 nt intronic + 3 nt exonic) and scores both alleles against pre-trained position weight matrices. The delta score (alt − ref) is used to classify the effect: large negative delta = disruption; large positive delta at a novel position = creation; intermediate = weakening. Branch point scanning uses a separate AG-distance heuristic.

**Outputs:** TSV of variant-level splice scores (`splice_scores.tsv`); optional delta score histogram.

**How it ships:** `python run_predictor.py variants.vcf --fasta genome.fa --gtf genes.gtf`; `main.py` delegates to `src.pipeline.main()` which loads `run_predictor.py` via `importlib`.

## Usage

```bash
# Predict splicing impact for a real VCF
python run_predictor.py variants.vcf --fasta genome.fa --gtf genes.gtf -o results/

# Synthetic demo (no real input required)
python run_predictor.py --demo -o results/

# Restrict to canonical GT-AG sites only
python run_predictor.py variants.vcf --fasta genome.fa --gtf genes.gtf --canonical-only -o results/
```

## Dependencies

- biopython>=1.81
- pandas>=2.0.0
- numpy>=1.24.0

## Status

Planned
