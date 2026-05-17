"""Jukes-Cantor (JC69) and Kimura 2-Parameter (K2P) pairwise distance models."""

from __future__ import annotations

import math

import numpy as np


PURINES = {"A", "G"}
PYRIMIDINES = {"C", "T"}


def compute_distance_matrix(
    sequences: list[tuple[str, str]],
    model: str = "K2P",
) -> tuple[list[str], np.ndarray]:
    """Compute pairwise distance matrix using the specified substitution model.

    Args:
        sequences: List of (name, sequence) pairs from alignment_reader.
        model: "JC69" or "K2P" (default).

    Returns:
        Tuple of (names list, symmetric distance matrix as ndarray).
    """
    names = [n for n, _ in sequences]
    seqs = [s for _, s in sequences]
    n = len(names)
    matrix = np.zeros((n, n), dtype=np.float64)

    dist_fn = _jc69 if model.upper() == "JC69" else _k2p

    for i in range(n):
        for j in range(i + 1, n):
            d = dist_fn(seqs[i], seqs[j])
            matrix[i, j] = d
            matrix[j, i] = d

    return names, matrix


def _jc69(seq1: str, seq2: str) -> float:
    """Jukes-Cantor distance: d = -3/4 * ln(1 - 4p/3)."""
    comparable = [(a, b) for a, b in zip(seq1, seq2) if a != "-" and b != "-" and a in "ACGT" and b in "ACGT"]
    if not comparable:
        return 0.0
    n = len(comparable)
    diffs = sum(1 for a, b in comparable if a != b)
    p = diffs / n
    if p >= 0.75:
        return 5.0
    arg = 1.0 - (4.0 * p / 3.0)
    if arg <= 0:
        return 5.0
    return -0.75 * math.log(arg)


def _k2p(seq1: str, seq2: str) -> float:
    """Kimura 2-Parameter distance: d = -1/2*ln(1-2P-Q) - 1/4*ln(1-2Q)."""
    comparable = [(a, b) for a, b in zip(seq1, seq2) if a != "-" and b != "-" and a in "ACGT" and b in "ACGT"]
    if not comparable:
        return 0.0
    n = len(comparable)
    transitions = 0
    transversions = 0
    for a, b in comparable:
        if a == b:
            continue
        if (a in PURINES and b in PURINES) or (a in PYRIMIDINES and b in PYRIMIDINES):
            transitions += 1
        else:
            transversions += 1

    P = transitions / n
    Q = transversions / n

    term1 = 1.0 - 2.0 * P - Q
    term2 = 1.0 - 2.0 * Q

    if term1 <= 0 or term2 <= 0:
        return _jc69(seq1, seq2)

    return -0.5 * math.log(term1) - 0.25 * math.log(term2)


def normalize_matrix(matrix: np.ndarray) -> np.ndarray:
    """Normalize distance matrix to [0, 1] by dividing by the max off-diagonal value."""
    mask = ~np.eye(matrix.shape[0], dtype=bool)
    max_val = matrix[mask].max() if mask.any() and matrix[mask].max() > 0 else 1.0
    return matrix / max_val
