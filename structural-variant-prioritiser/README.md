# Structural Variant Prioritiser

Prioritises structural variants by breakpoint gene impact, ClinGen dosage sensitivity, and overlap with known pathogenic SVs.

## Overview

Given a VCF or BED of structural variants, a gene annotation, and ClinGen dosage sensitivity scores, this tool intersects each SV with gene exons, retrieves haploinsufficiency (HI) and triplosensitivity (TS) scores for overlapping genes, queries for DECIPHER and ClinVar SV matches, and computes a composite priority score for three-tier ranking.

## Approach

**Inputs:** SV VCF or BED (SVTYPE, SVLEN, AF INFO fields); gene annotation BED or GTF; ClinGen dosage sensitivity TSV; optional DECIPHER and ClinVar SV databases.

**Core method:** SVs smaller than the minimum size or more common than the AF threshold are excluded. Each remaining SV is intersected with exon coordinates; SVs disrupting coding exons score higher. For DEL/DUP calls, the overlapping gene's ClinGen HI/TS score contributes the main weight. Database lookup against DECIPHER pathogenic CNVs and ClinVar pathogenic SVs (≥ 50 % reciprocal overlap) adds direct evidence. A weighted composite score determines tier assignment.

**Outputs:** TSV of prioritised SVs with tier, gene overlaps, dosage scores, and database matches (`prioritised_svs.tsv`); optional size vs. score scatter plot.

**How it ships:** `python run_prioritiser.py svs.vcf --gene-annotation genes.bed --dosage-scores clingen_hi.tsv`; `main.py` delegates to `src.pipeline.main()` which loads `run_prioritiser.py` via `importlib`.

## Usage

```bash
# Prioritise SVs with gene annotation and dosage scores
python run_prioritiser.py svs.vcf --gene-annotation genes.bed --dosage-scores clingen_hi.tsv -o results/

# Synthetic demo (no real input required)
python run_prioritiser.py --demo -o results/

# Only include deletions and duplications ≥ 1 kb
python run_prioritiser.py svs.vcf --gene-annotation genes.bed --sv-types DEL DUP --min-size 1000 -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0
- matplotlib>=3.7.0

## Status

Planned
