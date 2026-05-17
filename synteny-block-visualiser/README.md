# Synteny Block Visualiser

Identifies and visualises conserved synteny blocks between two or more genome assemblies, revealing large-scale chromosomal rearrangements, inversions, and translocations that distinguish species. In de-extinction genomics, synteny analysis informs which genomic regions of a proxy species are architecturally equivalent to the target, guiding the placement of engineered loci and the interpretation of genome-wide edit burden calculations.

## Inputs

- Two or more genome assemblies in FASTA format
- Species name labels corresponding to each assembly
- Optional: pre-computed pairwise alignment in PAF or MAF format (skips internal alignment step)
- Parameters: minimum synteny block length (bp), minimum number of anchors per block, dot-plot or ribbon diagram mode selection

## Outputs

- A synteny block table in BED-pair format (CSV) listing matched intervals between each genome pair with orientation and block identity
- A dot-plot visualisation for pairwise genome comparisons (PNG/SVG)
- A ribbon/Circos-style diagram for multi-genome synteny overview (PNG/SVG)
- A rearrangement summary listing inversions, translocations, and duplications detected between each pair

## Method

Identifies orthologous anchors between assemblies using k-mer based seeding or by parsing a provided alignment file. Chains collinear anchors into synteny blocks using a dynamic-programming chaining algorithm, filtering by minimum length and anchor density. Detects inversions as antiparallel chains and translocations as inter-chromosomal chains. Renders dot-plots with coloured orientation indicators and ribbon diagrams using proportional chromosome representations. Block coordinates are reported in both assembly coordinate systems.

## Dependencies

- `biopython` — FASTA I/O and k-mer indexing
- `matplotlib` — dot-plot and ribbon diagram rendering
- `pandas` — synteny block table and rearrangement summary
- `numpy` — chaining algorithm and coordinate arithmetic
- `pysam` — PAF/MAF alignment file parsing
