from dataclasses import dataclass
import numpy as np
from .genetic_code import CODON_TABLE, count_sites


@dataclass
class DNDSResult:
    """Result of a pairwise dN/dS calculation."""
    dN: float
    dS: float
    omega: float
    n_sites_S: float  # average synonymous sites (S̄)
    n_sites_N: float  # average non-synonymous sites (N̄)


def calculate_dnds(seq1, seq2):
    """Calculate dN/dS between two codon-aligned sequences using Nei-Gojobori (1986).

    Steps:
    1. Count synonymous (S) and non-synonymous (N) sites in each sequence;
       average to give S̄ and N̄.
    2. Count synonymous (Sd) and non-synonymous (Nd) differences between the
       two sequences.
    3. Apply Jukes-Cantor correction:
           pS = Sd / S̄,  pN = Nd / N̄
           dS = -3/4 * ln(1 - 4*pS/3)
           dN = -3/4 * ln(1 - 4*pN/3)
    4. omega = dN / dS.  If dS ≈ 0 and dN > 0 return np.inf; if dN ≈ 0
       return 0.

    Returns a DNDSResult dataclass.
    """
    S1, N1 = count_sites(seq1)
    S2, N2 = count_sites(seq2)
    S_bar = (S1 + S2) / 2.0
    N_bar = (N1 + N2) / 2.0

    Sd, Nd = count_differences(seq1, seq2)

    pS = Sd / S_bar if S_bar > 0 else 0.0
    pN = Nd / N_bar if N_bar > 0 else 0.0

    dS = jukes_cantor(pS)
    dN = jukes_cantor(pN)

    if dS < 1e-9:
        omega = np.inf if dN > 1e-9 else 0.0
    elif dN < 1e-9:
        omega = 0.0
    else:
        omega = dN / dS

    return DNDSResult(dN=dN, dS=dS, omega=omega, n_sites_S=S_bar, n_sites_N=N_bar)


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
