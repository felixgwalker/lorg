import random


def parse_fasta(fasta_path):
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
    random.seed(7)
    bases = list("ACGT")

    def rand_seq(n):
        return "".join(random.choices(bases, k=n))

    block_a = rand_seq(3000)
    block_b = rand_seq(2000)
    inversion_block = rand_seq(1500)
    translocation_block = rand_seq(1200)
    filler = [rand_seq(1000), rand_seq(800), rand_seq(600), rand_seq(900)]

    genome1 = {
        "g1_chr1": filler[0] + block_a + filler[1] + block_b,
        "g1_chr2": filler[2] + inversion_block + filler[3],
        "g1_chr3": rand_seq(2000) + translocation_block + rand_seq(800),
    }

    inv_rc = _reverse_complement(inversion_block)

    genome2 = {
        "g2_chr1": filler[0] + block_a + filler[1] + block_b,
        "g2_chr2": filler[2] + inv_rc + filler[3],
        "g2_chr3": translocation_block + rand_seq(2000) + rand_seq(800),
    }

    return genome1, genome2


def _reverse_complement(seq):
    comp = str.maketrans("ACGT", "TGCA")
    return seq.translate(comp)[::-1]
