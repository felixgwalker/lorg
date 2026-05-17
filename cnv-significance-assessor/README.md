# CNV Significance Assessor

Assesses the functional and population-level significance of copy number variants (CNVs) detected in a genome of interest, annotating each variant with gene content, dosage sensitivity, and cross-population frequency to distinguish likely benign polymorphisms from potentially deleterious changes. In conservation genomics and de-extinction feasibility studies, CNV assessment clarifies whether structural differences between proxy and target species — or among individuals of an endangered population — represent adaptive variation, neutral drift, or fitness-relevant losses and gains.

## Inputs

- CNV calls in BED or VCF format (deletions, duplications, or complex events)
- Reference genome gene annotation in GFF3 format
- Optional: population-level CNV frequency database in BED or VCF format (e.g. DGV, gnomAD-SV, or a custom panel)
- Optional: haploinsufficiency and triplosensitivity scores for genes of interest (CSV)
- Parameters: minimum CNV size filter, gene overlap fraction threshold, population frequency cutoff for filtering common variants

## Outputs

- An annotated CNV table in CSV format with columns: CNV ID, type, size, overlapping genes, dosage sensitivity score, population frequency, significance tier
- A significance classification summary (likely benign, variant of uncertain significance, likely pathogenic)
- A visualisation of CNV distribution across chromosomes with significance colour coding (PNG/SVG)
- A gene impact report listing all genes wholly or partially spanned by high-significance CNVs

## Method

Intersects CNV intervals with gene bodies and regulatory elements from the annotation. For each overlapping gene, retrieves haploinsufficiency (pLI, pHaplo) and triplosensitivity (pTriplo) scores if available. Cross-references CNV coordinates against the population frequency database to identify common variants. Classifies each CNV using a rule-based scoring system combining size, gene content, dosage sensitivity, and population frequency. Generates visualisations using chromosome ideogram plots.

## Dependencies

- `pandas` — CNV annotation table and classification logic
- `pybedtools` — interval intersection with gene and regulatory element annotations
- `pyvcf` — VCF format CNV input parsing
- `numpy` — score aggregation and filtering
- `matplotlib` — chromosome ideogram and significance distribution plots
