# Regulatory Element Conservation Scorer

Scores the conservation of non-coding regulatory elements — including promoters, enhancers, silencers, and insulator sequences — across a set of species, distinguishing elements that are functionally conserved (motif intact, accessible chromatin) from those that have diverged. In de-extinction projects, regulatory conservation scoring identifies which control regions of a proxy genome can be retained without modification and which must be edited to match target-species regulatory logic, directly informing the non-coding component of the edit burden.

## Inputs

- Regulatory element coordinates in BED format for the reference species
- Genome assemblies in FASTA format for each species to compare
- A transcription factor motif database in MEME or JASPAR format
- Optional: chromatin accessibility data (ATAC-seq or DNase-seq peaks) in BED or BigWig format for one or more species
- Parameters: motif score threshold, minimum alignment identity for ortholog detection, species subset selection

## Outputs

- A per-element conservation score table in CSV format with columns: element ID, coordinates, motif presence per species, sequence identity, optional epigenomic signal correlation
- A heatmap of conservation scores across elements and species (PNG/SVG)
- A ranked list of poorly conserved elements recommended for manual review or targeted editing
- A motif presence/absence matrix per element per species (CSV)

## Method

For each regulatory element, extracts the orthologous sequence in each target species via pairwise alignment. Scans all orthologous sequences for the presence of transcription factor binding motifs using position weight matrix (PWM) scoring. Computes a conservation score as a weighted combination of: pairwise sequence identity, motif retention rate across species, and if provided, cross-species ATAC-seq signal correlation. Elements are classified as conserved, partially conserved, or diverged based on configurable thresholds.

## Dependencies

- `biopython` — FASTA I/O, pairwise alignment, and motif scanning
- `pandas` — score table and motif matrix construction
- `numpy` — PWM scoring and weighted aggregation
- `pyBigWig` — ATAC-seq BigWig signal extraction
- `pybedtools` — BED interval operations and cross-species overlap
- `matplotlib` — heatmap and ranking visualisation
