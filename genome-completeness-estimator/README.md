# Genome Completeness Estimator

Estimates genome assembly completeness by searching BUSCO conserved gene benchmarks against a lineage-specific database.

## Overview

Given an assembly FASTA and a BUSCO lineage dataset, this tool runs HMMER searches for each BUSCO conserved orthogroup, classifies each as complete (single-copy or duplicated), fragmented, or missing, and produces the standard BUSCO summary statistics.

## Approach

**Inputs:** Assembly FASTA (genome, proteome, or transcriptome); BUSCO lineage dataset (e.g. vertebrata_odb10).

**Core method:** The assembly is scanned with HMMER profile HMMs for each BUSCO gene family in the selected lineage. Hits are clustered by locus and classified: complete single-copy (one locus, full-length hit); complete duplicated (≥ 2 loci); fragmented (partial hit covering < 95 % of the HMM length); missing (no hit above e-value threshold). The completeness fraction is reported as (complete_single + complete_duplicated) / n_buscos. Standard BUSCO summary line is reported (C:X%[S:Y%,D:Z%],F:A%,M:B%,n:N).

**Outputs:** BUSCO full table TSV (`busco_results.tsv`); summary statistics TSV (`busco_summary.tsv`); optional bar chart.

**How it ships:** `python run_estimator.py assembly.fa --lineage vertebrata_odb10`; `main.py` delegates to `src.pipeline.main()` which loads `run_estimator.py` via `importlib`.

## Usage

```bash
# Estimate completeness against vertebrate lineage
python run_estimator.py assembly.fa --lineage vertebrata_odb10 -o results/

# Synthetic demo (no real input required)
python run_estimator.py --demo -o results/

# Protein mode
python run_estimator.py proteins.fa --mode proteins --lineage vertebrata_odb10 -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- matplotlib>=3.7.0

## Status

Planned
