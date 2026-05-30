# Gene Model Validator

Validates gene models in a GTF annotation against a reference genome FASTA, checking start/stop codons, splice sites, and structural integrity.

## Overview

Given a GTF annotation and reference genome, this tool extracts each transcript's CDS, checks for structural validity, and reports all error types alongside the translated protein sequence for valid models.

## Approach

**Inputs:** GTF annotation; reference genome FASTA.

**Core method:** For each coding transcript, the CDS is assembled from exon coordinates and extracted from the FASTA. The pipeline checks: (1) presence of a start codon (ATG) at the CDS 5' end; (2) canonical stop codon (TAA/TAG/TGA) at the CDS 3' end; (3) no internal stop codons in the ORF; (4) canonical GT-AG splice sites at each intron boundary (GC-AG and AT-AC flagged as non-canonical unless `--allow-non-canonical`); (5) minimum intron length ≥ 60 bp; (6) no overlapping exons within the same transcript. Valid transcripts have their protein sequence reported.

**Outputs:** TSV of per-transcript validation results (`gene_model_validation.tsv`); FASTA of valid protein sequences (`valid_proteins.fa`); optional error-type bar chart.

**How it ships:** `python run_validator.py --annotation genes.gtf --fasta genome.fa`; `main.py` delegates to `src.pipeline.main()` which loads `run_validator.py` via `importlib`.

## Usage

```bash
# Validate gene models
python run_validator.py --annotation genes.gtf --fasta genome.fa -o results/

# Synthetic demo (no real input required)
python run_validator.py --demo -o results/
```

## Dependencies

- biopython>=1.81
- pandas>=2.0.0
- numpy>=1.24.0

## Status

Planned
