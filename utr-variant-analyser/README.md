# UTR Variant Analyser

Analyses variants in 5' and 3' UTR regions for uORF creation/disruption, Kozak context changes, and polyadenylation signal disruptions.

## Overview

Given a VCF of UTR variants, a reference genome FASTA, and a UTR annotation, this tool classifies each variant by its regulatory consequence: upstream ORF (uORF) creation or disruption in 5' UTRs, Kozak consensus strength changes, and hexamer polyadenylation signal disruption in 3' UTRs.

## Approach

**Inputs:** VCF of variants in UTR regions; reference genome FASTA; UTR annotation BED or GTF.

**Core method:** For **5' UTR** variants: the reference and alternate sequences within the UTR are scanned for ATG codons in all three frames to count uORFs (minimum 3 aa). Kozak context strength is scored by comparing the −3 and +4 positions to the GCCACC consensus. For **3' UTR** variants: canonical polyadenylation signals (AATAAA and variants) within 40 nt of the variant are checked for disruption. Each variant receives an effect classification and an impact summary string.

**Outputs:** TSV of UTR variant classifications (`utr_variants.tsv`); optional lollipop plot.

**How it ships:** `python run_analyser.py variants.vcf --fasta genome.fa --annotation utrs.bed`; `main.py` delegates to `src.pipeline.main()` which loads `run_analyser.py` via `importlib`.

## Usage

```bash
# Analyse a real VCF
python run_analyser.py variants.vcf --fasta genome.fa --annotation utrs.bed -o results/

# Synthetic demo (no real input required)
python run_analyser.py --demo -o results/

# Analyse 5' UTRs only
python run_analyser.py variants.vcf --fasta genome.fa --annotation utrs.bed --utr-type 5prime -o results/
```

## Dependencies

- biopython>=1.81
- pandas>=2.0.0
- numpy>=1.24.0

## Status

Planned
