"""Codon tables: standard genetic code, human and E. coli codon usage frequencies."""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Standard genetic code (NCBI transl_table=1)
# Maps DNA codon (uppercase) -> single-letter amino acid, or "*" for stop.
# ---------------------------------------------------------------------------
STANDARD_GENETIC_CODE: dict[str, str] = {
    # Phenylalanine (F)
    "TTT": "F", "TTC": "F",
    # Leucine (L)
    "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    # Isoleucine (I)
    "ATT": "I", "ATC": "I", "ATA": "I",
    # Methionine / Start (M)
    "ATG": "M",
    # Valine (V)
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    # Serine (S)
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "AGT": "S", "AGC": "S",
    # Proline (P)
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    # Threonine (T)
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    # Alanine (A)
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    # Tyrosine (Y)
    "TAT": "Y", "TAC": "Y",
    # Stop codons
    "TAA": "*", "TAG": "*", "TGA": "*",
    # Histidine (H)
    "CAT": "H", "CAC": "H",
    # Glutamine (Q)
    "CAA": "Q", "CAG": "Q",
    # Asparagine (N)
    "AAT": "N", "AAC": "N",
    # Lysine (K)
    "AAA": "K", "AAG": "K",
    # Aspartate (D)
    "GAT": "D", "GAC": "D",
    # Glutamate (E)
    "GAA": "E", "GAG": "E",
    # Cysteine (C)
    "TGT": "C", "TGC": "C",
    # Tryptophan (W)
    "TGG": "W",
    # Arginine (R)
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGA": "R", "AGG": "R",
    # Glycine (G)
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

assert len(STANDARD_GENETIC_CODE) == 64, (
    f"Genetic code table has {len(STANDARD_GENETIC_CODE)} entries, expected 64"
)


# ---------------------------------------------------------------------------
# Human codon usage frequencies
# Relative Synonymous Codon Usage (RSCU) values normalised to [0, 1] per
# synonymous family, derived from Kazusa Homo sapiens (hg, GenBank release).
# Values represent the fraction of usage within each synonymous group (tAI proxy).
# ---------------------------------------------------------------------------
HUMAN_CODON_FREQ: dict[str, float] = {
    # Phe (F)
    "TTT": 0.45, "TTC": 0.55,
    # Leu (L)
    "TTA": 0.07, "TTG": 0.13, "CTT": 0.13, "CTC": 0.20, "CTA": 0.07, "CTG": 0.40,
    # Ile (I)
    "ATT": 0.36, "ATC": 0.48, "ATA": 0.16,
    # Met (M) — single codon
    "ATG": 1.00,
    # Val (V)
    "GTT": 0.18, "GTC": 0.24, "GTA": 0.11, "GTG": 0.47,
    # Ser (S)
    "TCT": 0.15, "TCC": 0.22, "TCA": 0.15, "TCG": 0.06, "AGT": 0.15, "AGC": 0.24,
    # Pro (P)
    "CCT": 0.28, "CCC": 0.33, "CCA": 0.27, "CCG": 0.11,
    # Thr (T)
    "ACT": 0.24, "ACC": 0.36, "ACA": 0.28, "ACG": 0.12,
    # Ala (A)
    "GCT": 0.26, "GCC": 0.40, "GCA": 0.23, "GCG": 0.11,
    # Tyr (Y)
    "TAT": 0.43, "TAC": 0.57,
    # Stop codons
    "TAA": 0.28, "TAG": 0.20, "TGA": 0.52,
    # His (H)
    "CAT": 0.41, "CAC": 0.59,
    # Gln (Q)
    "CAA": 0.25, "CAG": 0.75,
    # Asn (N)
    "AAT": 0.46, "AAC": 0.54,
    # Lys (K)
    "AAA": 0.42, "AAG": 0.58,
    # Asp (D)
    "GAT": 0.46, "GAC": 0.54,
    # Glu (E)
    "GAA": 0.42, "GAG": 0.58,
    # Cys (C)
    "TGT": 0.45, "TGC": 0.55,
    # Trp (W) — single codon
    "TGG": 1.00,
    # Arg (R)
    "CGT": 0.08, "CGC": 0.19, "CGA": 0.11, "CGG": 0.21, "AGA": 0.20, "AGG": 0.20,
    # Gly (G)
    "GGT": 0.16, "GGC": 0.34, "GGA": 0.25, "GGG": 0.25,
}

assert len(HUMAN_CODON_FREQ) == 64, (
    f"Human codon frequency table has {len(HUMAN_CODON_FREQ)} entries, expected 64"
)

# Rare codon threshold and set (human)
RARE_CODON_THRESHOLD = 0.15
RARE_CODONS: set[str] = {c for c, f in HUMAN_CODON_FREQ.items() if f < RARE_CODON_THRESHOLD}


# ---------------------------------------------------------------------------
# E. coli K-12 codon usage frequencies
# Relative Synonymous Codon Usage (RSCU) normalised to [0, 1] per synonymous
# family, derived from Kazusa Escherichia coli K-12 (GenBank release).
# ---------------------------------------------------------------------------
ECOLI_CODON_FREQ: dict[str, float] = {
    # Phe (F)
    "TTT": 0.58, "TTC": 0.42,
    # Leu (L)
    "TTA": 0.14, "TTG": 0.13, "CTT": 0.12, "CTC": 0.10, "CTA": 0.04, "CTG": 0.47,
    # Ile (I)
    "ATT": 0.51, "ATC": 0.42, "ATA": 0.07,
    # Met (M)
    "ATG": 1.00,
    # Val (V)
    "GTT": 0.28, "GTC": 0.20, "GTA": 0.17, "GTG": 0.35,
    # Ser (S)
    "TCT": 0.17, "TCC": 0.15, "TCA": 0.14, "TCG": 0.14, "AGT": 0.16, "AGC": 0.25,
    # Pro (P)
    "CCT": 0.18, "CCC": 0.13, "CCA": 0.20, "CCG": 0.49,
    # Thr (T)
    "ACT": 0.19, "ACC": 0.40, "ACA": 0.17, "ACG": 0.25,
    # Ala (A)
    "GCT": 0.18, "GCC": 0.26, "GCA": 0.23, "GCG": 0.33,
    # Tyr (Y)
    "TAT": 0.59, "TAC": 0.41,
    # Stop codons
    "TAA": 0.61, "TAG": 0.09, "TGA": 0.30,
    # His (H)
    "CAT": 0.57, "CAC": 0.43,
    # Gln (Q)
    "CAA": 0.34, "CAG": 0.66,
    # Asn (N)
    "AAT": 0.49, "AAC": 0.51,
    # Lys (K)
    "AAA": 0.74, "AAG": 0.26,
    # Asp (D)
    "GAT": 0.63, "GAC": 0.37,
    # Glu (E)
    "GAA": 0.68, "GAG": 0.32,
    # Cys (C)
    "TGT": 0.46, "TGC": 0.54,
    # Trp (W)
    "TGG": 1.00,
    # Arg (R)
    "CGT": 0.36, "CGC": 0.36, "CGA": 0.07, "CGG": 0.11, "AGA": 0.07, "AGG": 0.04,
    # Gly (G)
    "GGT": 0.35, "GGC": 0.37, "GGA": 0.13, "GGG": 0.15,
}

assert len(ECOLI_CODON_FREQ) == 64, (
    f"E. coli codon frequency table has {len(ECOLI_CODON_FREQ)} entries, expected 64"
)


def get_codon_usage_human() -> dict[str, float]:
    """Return human codon usage frequencies (tAI proxy) for all 64 codons.

    Values are Relative Synonymous Codon Usage (RSCU) normalised within each
    synonymous family to [0, 1], derived from Kazusa Homo sapiens database.
    A value of 1.0 indicates the sole or most preferred codon; rare codons
    approach 0.

    Returns
    -------
    dict[str, float]
        Mapping of uppercase DNA codon triplet to relative usage frequency.
    """
    return dict(HUMAN_CODON_FREQ)


def get_codon_usage_ecoli() -> dict[str, float]:
    """Return E. coli K-12 codon usage frequencies (tAI proxy) for all 64 codons.

    Values are Relative Synonymous Codon Usage (RSCU) normalised within each
    synonymous family to [0, 1], derived from Kazusa Escherichia coli K-12
    database. High values indicate preferred codons matching abundant tRNA
    pools; low values indicate rare, potentially pausing codons.

    Returns
    -------
    dict[str, float]
        Mapping of uppercase DNA codon triplet to relative usage frequency.
    """
    return dict(ECOLI_CODON_FREQ)
