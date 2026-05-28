import os
import numpy as np


def read_codon_alignment(fasta_path):
    """Parse a codon-aligned nucleotide FASTA file.

    Returns a list of (name, seq) tuples.  Validates that:
    - Every sequence length is a multiple of 3.
    - All sequences have the same length (aligned).
    - No gap characters break codon boundaries (gaps only at codon-triplet
      positions, i.e. gap runs whose length is a multiple of 3 *and* aligned
      to codon boundaries).

    Raises ValueError on validation failure.
    """
    sequences = []
    current_id = None
    current_seq = []
    with open(fasta_path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    sequences.append((current_id, "".join(current_seq).upper()))
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())
    if current_id is not None:
        sequences.append((current_id, "".join(current_seq).upper()))

    if not sequences:
        raise ValueError(f"No sequences found in {fasta_path}")

    # Validate lengths
    lengths = {len(seq) for _, seq in sequences}
    if len(lengths) > 1:
        raise ValueError(
            f"Sequences in {fasta_path} have different lengths: {sorted(lengths)}"
        )
    aln_len = lengths.pop()
    if aln_len % 3 != 0:
        raise ValueError(
            f"Alignment length {aln_len} is not a multiple of 3 in {fasta_path}"
        )

    # Validate that gaps do not break codon boundaries
    for name, seq in sequences:
        for codon_start in range(0, aln_len, 3):
            codon = seq[codon_start:codon_start + 3]
            gap_count = codon.count("-")
            if 0 < gap_count < 3:
                raise ValueError(
                    f"Gap breaks codon boundary at position {codon_start} "
                    f"in sequence '{name}' (codon: '{codon}')"
                )

    return sequences


def make_demo_alignment():
    """Return a short demo codon alignment as a list of (name, seq) tuples.

    The alignment contains 5 sequences over 15 codons (45 nucleotides),
    with a mix of synonymous and non-synonymous differences to exercise the
    dN/dS machinery.
    """
    seqs = [
        ("species_1", "ATGGCCAAAGTTCTGCAGCACGACTTCAACGGTTCGTAA"),
        ("species_2", "ATGGCCAAAGTTCTGCAGCACGACTTCAACGGTTCGTAA"),
        ("species_3", "ATGGCTAAAGTTTTGCAACACGACTTCAATGGCTCGTAA"),
        ("species_4", "ATGGCCAAAGTTCTGCAGCACGATTTCAACGGTTCGTAA"),
        ("species_5", "ATGGCCAAGGTTCTGCAGCATGACTTCAACGGTTCCTAA"),
    ]
    # Pad/trim so each is exactly 39 nt (13 codons, multiple of 3)
    result = []
    for name, seq in seqs:
        # Ensure length is multiple of 3
        trimmed = seq[: (len(seq) // 3) * 3]
        result.append((name, trimmed))
    return result


def parse_fasta(path):
    sequences = {}
    current_id = None
    current_seq = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())
    if current_id is not None:
        sequences[current_id] = "".join(current_seq)
    return sequences


def load_gene_alignments(input_path):
    genes = {}
    if os.path.isdir(input_path):
        for fname in sorted(os.listdir(input_path)):
            if fname.endswith(".fa") or fname.endswith(".fasta") or fname.endswith(".fas"):
                gene_name = os.path.splitext(fname)[0]
                fpath = os.path.join(input_path, fname)
                seqs = parse_fasta(fpath)
                if seqs:
                    genes[gene_name] = seqs
    else:
        seqs = parse_fasta(input_path)
        genes["gene_1"] = seqs
    return genes


def generate_synthetic_alignments(n_genes=20, n_species=5, seq_len=300, n_positive=3, rng=None):
    if rng is None:
        rng = np.random.default_rng(42)

    bases = list("ACGT")
    codon_table = _get_codon_table()
    genes = {}

    for gene_idx in range(n_genes):
        gene_name = f"gene_{gene_idx+1:02d}"
        is_positive = gene_idx < n_positive

        seq_len_codons = seq_len // 3
        ref_codons = []
        for _ in range(seq_len_codons):
            codon = "".join(rng.choice(list("ACGT"), size=3))
            ref_codons.append(codon)
        ref_seq = "".join(ref_codons)

        seqs = {"species_1": ref_seq}
        for sp_idx in range(1, n_species):
            if is_positive:
                mutation_rate = 0.08
                syn_bias = 0.2
            else:
                mutation_rate = 0.04
                syn_bias = 0.8

            mutant_codons = list(ref_codons)
            for ci, codon in enumerate(ref_codons):
                if rng.random() < mutation_rate:
                    if rng.random() < syn_bias:
                        mutant_codons[ci] = _synonymous_mutant(codon, rng, codon_table)
                    else:
                        mutant_codons[ci] = _nonsynonymous_mutant(codon, rng, codon_table)
            seqs[f"species_{sp_idx+1}"] = "".join(mutant_codons)

        genes[gene_name] = seqs

    return genes


def _get_codon_table():
    from .genetic_code import CODON_TABLE
    return CODON_TABLE


def _synonymous_mutant(codon, rng, codon_table):
    from .genetic_code import BASES
    ref_aa = codon_table.get(codon, "X")
    candidates = []
    for pos in range(3):
        for base in BASES:
            if base == codon[pos]:
                continue
            mutant = codon[:pos] + base + codon[pos+1:]
            if codon_table.get(mutant) == ref_aa:
                candidates.append(mutant)
    if candidates:
        return candidates[rng.integers(0, len(candidates))]
    return codon


def _nonsynonymous_mutant(codon, rng, codon_table):
    from .genetic_code import BASES
    ref_aa = codon_table.get(codon, "X")
    candidates = []
    for pos in range(3):
        for base in BASES:
            if base == codon[pos]:
                continue
            mutant = codon[:pos] + base + codon[pos+1:]
            mut_aa = codon_table.get(mutant, "*")
            if mut_aa != ref_aa and mut_aa != "*":
                candidates.append(mutant)
    if candidates:
        return candidates[rng.integers(0, len(candidates))]
    return codon
