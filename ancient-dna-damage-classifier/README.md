# Ancient DNA Damage Classifier

Analyses sequencing reads from ancient or historically derived specimens to characterise post-mortem DNA damage patterns and classify reads as authentically ancient, potentially modern-contaminated, or ambiguous. Accurate damage profiling is essential for validating ancient genome assemblies used as reference targets in de-extinction projects, and for distinguishing true ancestral variants from artefactual substitutions introduced by chemical degradation.

**Scope note:** This tool absorbs the planned `dna-fragmentation-profiler` module. Fragment length distribution profiling and terminal-base deamination pattern visualisation (C→T at 5' / G→A at 3') are included here alongside the Bayesian per-read classifier, so both concerns live in a single aDNA characterisation step.

## Inputs

- Aligned reads in BAM format mapped to a reference genome, or raw reads in FASTQ format
- Reference genome in FASTA format (required if FASTQ is provided)
- Optional mapDamage2 output tables for comparison or validation
- Parameters: minimum mapping quality, read length filter, number of terminal bases to profile

## Outputs

- Damage frequency tables (C→T and G→A rates per read position) in CSV format
- A damage profile plot showing substitution frequency at 5' and 3' read termini (PNG/SVG)
- A per-read classification output in TSV format (authentic / contaminated / ambiguous) with posterior probability
- A summary statistics report including mean fragment length, deamination rates, and overall authenticity estimate

## Method

Computes position-specific C→T substitution frequencies at the 5' end and G→A at the 3' end, which are the hallmarks of cytosine deamination in ancient DNA. Fits a geometric decay model to observed damage rates. Applies a Bayesian classifier to assign each read a posterior probability of being authentically ancient, using a prior derived from the library-level damage profile. Fragment length distribution is also computed as a secondary authenticity indicator.

## Dependencies

- `pysam` — BAM/SAM file parsing and read iteration
- `biopython` — FASTQ I/O and sequence utilities
- `numpy` — positional frequency array construction
- `pandas` — tabular output handling
- `scipy` — curve fitting for geometric decay model and Bayesian inference
- `matplotlib` — damage profile visualisation
