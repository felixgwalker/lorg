# ERV Risk Mapper

Identifies and risk-scores endogenous retrovirus (ERV) elements within a target genome assembly, assessing their potential for reactivation or mobilisation during genome engineering. In de-extinction and gene therapy contexts, inadvertent activation of ERV loci can pose biosafety risks or disrupt gene expression; this tool provides a systematic landscape of ERV locations, families, estimated age, and activity risk before editing campaigns begin.

## Inputs

- Target genome assembly in FASTA format
- Optional genome annotation file in GFF3 or BED format
- An ERV/repeat reference database (e.g. Dfam HMM library or a custom RepeatMasker library in FASTA)
- A minimum LTR identity threshold parameter for classifying recently active elements

## Outputs

- A BED file of all annotated ERV loci with family and subfamily labels
- A per-element activity risk score table in CSV format (columns: locus, family, LTR integrity, estimated insertion age, risk tier)
- A genome-wide risk map plot showing ERV density and risk tiers per chromosome (PNG/SVG)
- A summary report listing high-risk loci recommended for manual review

## Method

Aligns the query genome against the ERV reference database using profile HMM searches. Clusters hits into full-length or solo-LTR elements. Estimates insertion age from the LTR-LTR divergence under a neutral substitution rate. Scores activity risk based on LTR completeness, open reading frame integrity, and proximity to active chromatin regions if annotation is supplied. Elements are tiered into low, moderate, and high risk categories.

## Dependencies

- `biopython` — FASTA I/O and sequence utilities
- `pandas` — locus table construction and filtering
- `numpy` — numerical scoring
- `matplotlib` — risk map and chromosome density plots
- `pybedtools` — BED interval operations and overlap analysis
- `hmmer` (external binary, called via subprocess) — profile HMM searches against ERV database
