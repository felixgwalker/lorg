import numpy as np
from .genetic_code import CODON_TABLE, count_sites


def jukes_cantor(p):
    if p <= 0:
        return 0.0
    if p >= 0.75:
        return 10.0
    val = 1 - (4 * p / 3)
    if val <= 0:
        return 10.0
    return -0.75 * np.log(val)


def count_differences(seq1, seq2):
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    min_len = min(len(seq1), len(seq2))
    min_len = (min_len // 3) * 3
    Sd = 0.0
    Nd = 0.0

    for i in range(0, min_len, 3):
        c1 = seq1[i:i+3]
        c2 = seq2[i:i+3]
        if len(c1) < 3 or len(c2) < 3:
            continue
        if "N" in c1 or "N" in c2 or "-" in c1 or "-" in c2:
            continue
        aa1 = CODON_TABLE.get(c1, "X")
        aa2 = CODON_TABLE.get(c2, "X")
        if aa1 == "*" or aa2 == "*":
            continue
        if c1 == c2:
            continue
        if aa1 == aa2:
            Sd += 1
        else:
            Nd += 1

    return Sd, Nd


def pairwise_dnds(seq1, seq2):
    S1, N1 = count_sites(seq1)
    S2, N2 = count_sites(seq2)
    S = (S1 + S2) / 2
    N = (N1 + N2) / 2

    Sd, Nd = count_differences(seq1, seq2)

    pS = Sd / S if S > 0 else 0.0
    pN = Nd / N if N > 0 else 0.0

    dS = jukes_cantor(pS)
    dN = jukes_cantor(pN)

    if dS < 1e-9:
        omega = 10.0 if dN > 0 else 1.0
    else:
        omega = min(dN / dS, 10.0)

    return dN, dS, omega, S, N, Sd, Nd


def compute_gene_dnds(sequences):
    species = list(sequences.keys())
    dN_vals = []
    dS_vals = []
    omega_vals = []

    for i in range(len(species)):
        for j in range(i + 1, len(species)):
            s1 = sequences[species[i]]
            s2 = sequences[species[j]]
            dN, dS, omega, S, N, Sd, Nd = pairwise_dnds(s1, s2)
            dN_vals.append(dN)
            dS_vals.append(dS)
            omega_vals.append(omega)

    if not omega_vals:
        return {"mean_dN": 0.0, "mean_dS": 0.0, "omega": 1.0, "n_pairs": 0}

    return {
        "mean_dN": float(np.mean(dN_vals)),
        "mean_dS": float(np.mean(dS_vals)),
        "omega": float(np.mean(omega_vals)),
        "n_pairs": len(omega_vals),
    }
