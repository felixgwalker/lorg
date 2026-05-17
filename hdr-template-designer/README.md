# HDR Template Designer

Designs homology-directed repair (HDR) donor templates for precise CRISPR-mediated genome edits, generating ready-to-synthesise DNA sequences with correctly configured homology arms flanking the desired edit. Precise HDR is the mechanism of choice for installing point mutations, short insertions, or multi-nucleotide changes in de-extinction and conservation genomics editing campaigns where exact sequence outcome is required.

## Inputs

- Reference genome region in FASTA format (a window around the intended cut site)
- Desired edit specification: either a VCF-format variant record or a short replacement sequence string
- Cut site coordinates in BED format or as a genomic coordinate string (chr:start-end)
- Parameters: homology arm length (default 100 bp), minimum GC content threshold, repeat-masking flag, silent mutation flag for PAM disruption

## Outputs

- HDR donor template sequence in FASTA format, ready for synthesis as ssDNA oligo or dsDNA fragment
- Annotated template diagram showing homology arm boundaries, edit position, and PAM-disrupting silent mutations (PNG/SVG)
- A QC report flagging potential issues: low GC, homopolymer runs, repetitive elements within arms, secondary structure risk
- Alternative arm-length variants table for optimisation

## Method

Extracts flanking sequences of the specified arm length on each side of the cut site from the reference genome. Inserts the desired edit at the appropriate position. Optionally introduces a silent PAM-disrupting mutation to prevent Cas9 from re-cutting after successful HDR. Screens homology arms for repetitive sequences using repeat content analysis. Calculates GC content and flags arms outside the 40–60% optimal range. Outputs the full template as a single FASTA record with embedded annotations in the header.

## Dependencies

- `biopython` — FASTA I/O, sequence manipulation, and complement operations
- `pysam` — reference genome access and region extraction
- `pandas` — QC report and arm variant table
- `numpy` — GC content and secondary structure calculations
- `matplotlib` — annotated template diagram rendering
