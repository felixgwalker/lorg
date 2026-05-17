# Guide RNA Off-Target Scorer

Predicts and scores potential off-target cleavage sites across a genome for a given set of CRISPR guide RNA sequences, enabling researchers to select guides with the highest on-target specificity before proceeding to experimental validation. This is particularly important in genome engineering projects for endangered or extinct-proxy species, where off-target edits in germ-line cells are irreversible and may have cascading phenotypic effects.

## Inputs

- Guide RNA sequences in FASTA or CSV format (one or more guides)
- Target genome assembly in FASTA format
- PAM sequence string (e.g. `NGG` for SpCas9, `TTTN` for AsCas12a)
- Parameters: maximum number of mismatches to enumerate (typically 3–5), seed region length, bulge tolerance (RNA or DNA)

## Outputs

- A genome-wide off-target site table in CSV format with columns: guide ID, chromosome, position, strand, mismatches, mismatch positions, CFD score
- A per-guide specificity summary (on-target rank, number of predicted off-targets by mismatch tier)
- A Manhattan-style plot of off-target scores across the genome (PNG/SVG)
- A ranked guide comparison table for guide selection

## Method

Performs exhaustive seed-region search of the genome for sequences matching each guide within the specified mismatch budget, using a hash-index of genome k-mers. Scores each candidate off-target site using the Cutting Frequency Determination (CFD) model, which weights mismatches by position and nucleotide identity based on empirical cleavage data. Optionally computes MIT specificity scores for cross-reference. Sites are ranked and filtered by score threshold.

## Dependencies

- `biopython` — FASTA I/O and sequence complement operations
- `numpy` — mismatch matrix construction and scoring
- `pandas` — off-target table construction and filtering
- `matplotlib` — Manhattan plot and guide comparison visualisation
- `pysam` — optional genome index access for large assemblies
