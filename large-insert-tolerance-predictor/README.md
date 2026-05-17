# Large Insert Tolerance Predictor

Predicts whether a specified genomic locus can tolerate the insertion of a large exogenous DNA segment — such as a multi-gene expression cassette or a reconstructed genomic interval — without disrupting local gene expression, chromatin architecture, or regulatory networks. This tool addresses the practical challenge that large inserts used in de-extinction workflows (e.g. inserting entire mammoth-specific gene clusters) carry a higher risk of position-effect variegation and regulatory interference than single-nucleotide edits.

## Inputs

- Target locus coordinates in BED format (chromosome, start, end)
- Reference genome assembly in FASTA format
- Insert size in base pairs
- Optional: ATAC-seq or FAIRE-seq chromatin accessibility signal in BigWig or BED format
- Optional: existing gene annotation in GFF3 format
- Parameters: window radius around insertion site (default 500 kb), gene density weight, regulatory proximity penalty radius

## Outputs

- A tolerance score (0–100) for each candidate insertion site within the specified window
- A feature annotation table listing nearby genes, enhancers, CTCF sites, and repeat elements per candidate site (CSV)
- A locus context visualisation showing the scored insertion window with annotated features (PNG/SVG)
- A ranked list of recommended insertion sites with rationale

## Method

Analyses a window around the target locus for features that predict insertion intolerance: gene density, proximity to known regulatory elements (promoters, enhancers, insulator sequences), repetitive element content, and, if provided, chromatin accessibility signal. Each feature contributes a weighted penalty or bonus to the tolerance score. Sites within DNase-I hypersensitive or ATAC-seq accessible regions are penalised as likely active regulatory space. Sites in gene deserts with low regulatory element density are scored as high tolerance. Scores are aggregated into a composite tolerance index.

## Dependencies

- `biopython` — FASTA I/O and sequence feature extraction
- `pybedtools` — interval intersection and annotation overlap
- `pyBigWig` — BigWig chromatin accessibility signal reading
- `pandas` — feature annotation table construction
- `numpy` — weighted scoring calculations
- `matplotlib` — locus context and score visualisation
