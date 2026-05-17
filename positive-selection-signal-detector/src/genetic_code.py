CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

BASES = ["T", "C", "A", "G"]


def synonymous_sites_per_codon(codon):
    codon = codon.upper()
    if codon not in CODON_TABLE or CODON_TABLE[codon] == "*":
        return 0.0, 1.0
    ref_aa = CODON_TABLE[codon]
    s_sites = 0.0
    n_sites = 0.0
    for pos in range(3):
        syn = 0
        total = 0
        for base in BASES:
            if base == codon[pos]:
                continue
            mutant = codon[:pos] + base + codon[pos+1:]
            if mutant not in CODON_TABLE:
                continue
            mut_aa = CODON_TABLE[mutant]
            total += 1
            if mut_aa == ref_aa:
                syn += 1
        if total > 0:
            frac_syn = syn / total
            s_sites += frac_syn / 3.0
            n_sites += (1 - frac_syn) / 3.0
        else:
            n_sites += 1.0 / 3.0
    return s_sites, n_sites


def count_sites(sequence):
    seq = sequence.upper().replace("-", "")
    S_total = 0.0
    N_total = 0.0
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if len(codon) < 3 or "N" in codon:
            continue
        s, n = synonymous_sites_per_codon(codon)
        S_total += s
        N_total += n
    return S_total, N_total
