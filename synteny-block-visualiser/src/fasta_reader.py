import random


def parse_fasta(fasta_path):
    """Parse a FASTA file and return {seq_name: sequence_str}."""
    sequences = {}
    current = None
    parts = []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if current is not None:
                    sequences[current] = "".join(parts)
                current = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.upper())
    if current is not None:
        sequences[current] = "".join(parts)
    return sequences


def make_demo_genomes():
    """Return two synthetic genome dicts demonstrating synteny rearrangements.

    Genome 1: 3 chromosomes x 100 kb each.
    Genome 2: genome1 with:
      - chromosome-scale inversion of chr2 (entire 100 kb reversed)
      - translocation of 20 kb from chr1 (positions 40000-60000) inserted at
        start of chr3 (chr3 original sequence follows)
      - ~2% random SNP noise applied across all chromosomes
    """
    random.seed(42)
    bases = list("ACGT")
    chrom_len = 100_000

    def rand_seq(n):
        return "".join(random.choices(bases, k=n))

    def apply_snp_noise(seq, rate=0.02):
        seq = list(seq)
        for i in range(len(seq)):
            if random.random() < rate:
                seq[i] = random.choice([b for b in bases if b != seq[i]])
        return "".join(seq)

    chr1 = rand_seq(chrom_len)
    chr2 = rand_seq(chrom_len)
    chr3 = rand_seq(chrom_len)

    genome1 = {
        "g1_chr1": chr1,
        "g1_chr2": chr2,
        "g1_chr3": chr3,
    }

    # Chromosome-scale inversion on chr2: reverse-complement the whole thing
    chr2_inverted = _reverse_complement(chr2)

    # Translocation: move 20 kb from chr1 (positions 40000-60000) to start of chr3
    trans_block = chr1[40_000:60_000]
    chr1_after_trans = chr1[:40_000] + chr1[60_000:]   # chr1 with the 20 kb removed
    chr3_after_trans = trans_block + chr3               # chr3 with the 20 kb prepended

    # Apply ~2% SNP noise
    genome2 = {
        "g2_chr1": apply_snp_noise(chr1_after_trans),
        "g2_chr2": apply_snp_noise(chr2_inverted),
        "g2_chr3": apply_snp_noise(chr3_after_trans),
    }

    return genome1, genome2


def _reverse_complement(seq):
    comp = str.maketrans("ACGT", "TGCA")
    return seq.translate(comp)[::-1]
