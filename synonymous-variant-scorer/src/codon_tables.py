"""Hardcoded human codon usage frequencies (fraction of synonymous codons used per amino acid)."""

HUMAN_CODON_FREQ: dict[str, float] = {
    "TTT": 0.45, "TTC": 0.55,
    "TTA": 0.07, "TTG": 0.13, "CTT": 0.13, "CTC": 0.20, "CTA": 0.07, "CTG": 0.40,
    "ATT": 0.36, "ATC": 0.48, "ATA": 0.16,
    "ATG": 1.00,
    "GTT": 0.18, "GTC": 0.24, "GTA": 0.11, "GTG": 0.47,
    "TCT": 0.15, "TCC": 0.22, "TCA": 0.15, "TCG": 0.06, "AGT": 0.15, "AGC": 0.24,
    "CCT": 0.28, "CCC": 0.33, "CCA": 0.27, "CCG": 0.11,
    "ACT": 0.24, "ACC": 0.36, "ACA": 0.28, "ACG": 0.12,
    "GCT": 0.26, "GCC": 0.40, "GCA": 0.23, "GCG": 0.11,
    "TAT": 0.43, "TAC": 0.57,
    "TAA": 0.28, "TAG": 0.20, "TGA": 0.52,
    "CAT": 0.41, "CAC": 0.59,
    "CAA": 0.25, "CAG": 0.75,
    "AAT": 0.46, "AAC": 0.54,
    "AAA": 0.42, "AAG": 0.58,
    "GAT": 0.46, "GAC": 0.54,
    "GAA": 0.42, "GAG": 0.58,
    "TGT": 0.45, "TGC": 0.55,
    "TGG": 1.00,
    "CGT": 0.08, "CGC": 0.19, "CGA": 0.11, "CGG": 0.21, "AGA": 0.20, "AGG": 0.20,
    "GGT": 0.16, "GGC": 0.34, "GGA": 0.25, "GGG": 0.25,
}

RARE_CODON_THRESHOLD = 0.15

RARE_CODONS: set[str] = {c for c, f in HUMAN_CODON_FREQ.items() if f < RARE_CODON_THRESHOLD}
