# Promoter Variant Scorer

Scores variants in promoter regions for transcription factor binding site disruption using position weight matrices.

## Overview

Given a VCF of promoter variants, a reference genome FASTA, and a JASPAR-format PWM database, this tool scans the sequence context around each variant against all PWMs above an information content threshold, computes reference and alternate allele scores, and reports TFBS disruptions and de novo TFBS creations.

## Approach

**Inputs:** VCF of variants within the promoter window (default 2 kb upstream of TSS); reference genome FASTA; PWM database in MEME or TRANSFAC format.

**Core method:** For each variant, the surrounding sequence (matrix width centred on the variant position) is extracted. Both the reference and alternate alleles are scored against each PWM using log-odds matrices calibrated to a background nucleotide frequency. A delta score (alt − ref) is computed; variants exceeding the disruption or creation threshold are flagged. PWMs with information content < 8 bits are excluded to avoid low-specificity hits. The TATA box window (−35 to −25 from TSS) is always scanned regardless of IC.

**Outputs:** TSV of TFBS events (`tfbs_events.tsv`); optional heatmap of delta scores by TF family.

**How it ships:** `python run_scorer.py variants.vcf --fasta genome.fa --pwm-db jaspar.meme`; `main.py` delegates to `src.pipeline.main()` which loads `run_scorer.py` via `importlib`.

## Usage

```bash
# Score promoter variants against JASPAR PWMs
python run_scorer.py variants.vcf --fasta genome.fa --pwm-db jaspar.meme -o results/

# Synthetic demo (no real input required)
python run_scorer.py --demo -o results/

# Extend promoter window to 5 kb
python run_scorer.py variants.vcf --fasta genome.fa --pwm-db jaspar.meme --promoter-window 5000 -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0

## Status

Planned
