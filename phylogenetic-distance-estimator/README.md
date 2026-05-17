# Phylogenetic Distance Estimator

Estimates evolutionary distance between a target extinct or endangered species and a set of candidate donor or proxy species. Useful for prioritising which extant species to use as a genomic reference or cell donor in de-extinction and conservation genomics workflows. Outputs a pairwise distance matrix and an optional phylogenetic tree, enabling researchers to rapidly rank species by genomic proximity without requiring a full-scale phylogenomic study.

## Inputs

- One or more multiple-sequence alignment files in FASTA or CLUSTAL format, covering shared orthologous gene regions or whole-genome alignments
- A plain-text list of species names corresponding to the sequences
- A substitution model selection parameter (e.g. `JC69`, `K2P`, `GTR+G`)
- Optional: a calibration file specifying node ages in Mya for time-scaled divergence estimation

## Outputs

- A pairwise distance matrix in CSV format
- A Newick-format phylogenetic tree file
- A matplotlib heatmap and dendrogram visualisation (PNG/SVG)
- A ranked species comparison table listing candidate donors ordered by proximity to the target

## Method

Performs multiple sequence alignment if unaligned input is provided, using a progressive alignment strategy. Computes pairwise substitution distances under the selected nucleotide or amino acid model. Constructs a distance-based tree using neighbour-joining. Optionally applies a molecular clock calibration to convert branch lengths to divergence times. Distance scores are normalised to a 0–1 scale for downstream ranking.

## Dependencies

- `biopython` — sequence I/O, alignment, and distance calculation
- `ete3` — tree construction, manipulation, and visualisation
- `scipy` — hierarchical clustering for dendrogram generation
- `pandas` — distance matrix and ranking table handling
- `matplotlib` — heatmap and dendrogram plotting
- `numpy` — numerical operations
