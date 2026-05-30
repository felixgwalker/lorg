# Palaeogenomic Coverage Assessor

Assesses genome-wide coverage statistics for ancient DNA BAMs: mapping rate, depth, breadth at multiple thresholds, endogenous fraction, and duplication rate.

## Overview

Given a BAM file of aligned ancient DNA reads, this tool computes a comprehensive set of coverage statistics per chromosome and genome-wide, classifying the dataset's utility for downstream analyses (high/medium/low/insufficient coverage).

## Approach

**Inputs:** Sorted and indexed BAM file of aligned ancient DNA reads.

**Core method:** Read-level statistics (total reads, mapped reads, duplicate reads, mapping quality distribution) are extracted from BAM flags and MAPQ fields. Per-base depth is computed with configurable MAPQ (default 25) and base quality (default 20) thresholds, with optional duplicate exclusion. Breadth metrics report the fraction of the reference covered at ≥ 1×, 5×, 10×, and 20× depth. Endogenous fraction is the mapping rate after quality filters. Per-chromosome depth and breadth are reported. Coverage class is assigned by mean depth: ≥ 5× (high), 1–5× (medium), 0.1–1× (low), < 0.1× (insufficient).

**Outputs:** TSV of genome-wide coverage statistics (`coverage_summary.tsv`); per-chromosome TSV (`per_chromosome_coverage.tsv`); optional depth histogram plot.

**How it ships:** `python run_assessor.py sample.bam`; `main.py` delegates to `src.pipeline.main()` which loads `run_assessor.py` via `importlib`.

## Usage

```bash
# Assess coverage for an ancient DNA BAM
python run_assessor.py sample.bam -o results/

# Synthetic demo (no real input required)
python run_assessor.py --demo -o results/

# Keep duplicates and lower MAPQ threshold
python run_assessor.py sample.bam --keep-duplicates --min-mapq 20 -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- matplotlib>=3.7.0

## Status

Planned
