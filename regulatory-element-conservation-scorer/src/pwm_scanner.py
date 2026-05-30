import math
import numpy as np

# ---------------------------------------------------------------------------
# Numpy PWM matrices (4 × L, rows: A C G T) for the 5 canonical motifs
# required by the public API.  Values are log-odds scores relative to a
# uniform 0.25 background — positive means enriched, negative depleted.
# ---------------------------------------------------------------------------
_BASES = "ACGT"
_BASE_IDX = {b: i for i, b in enumerate(_BASES)}


def _prob_to_logodds(prob_matrix):
    """Convert a 4×L probability array to log2 log-odds vs uniform 0.25."""
    return np.log2(np.maximum(prob_matrix, 1e-10) / 0.25)


def _build_numpy_pwm(prob_rows):
    """Build a 4×L numpy log-odds array from a list of {A,C,G,T} dicts."""
    L = len(prob_rows)
    mat = np.zeros((4, L), dtype=float)
    for col_idx, col in enumerate(prob_rows):
        for base, row_idx in _BASE_IDX.items():
            mat[row_idx, col_idx] = math.log2(max(col.get(base, 1e-10), 1e-10) / 0.25)
    return mat


# Five hardcoded numpy PWMs (SP1, CTCF, TATA-box, E-box, AP1).
NUMPY_PWMS = {
    "SP1": _build_numpy_pwm([
        {"A": 0.03, "C": 0.05, "G": 0.89, "T": 0.03},
        {"A": 0.03, "C": 0.05, "G": 0.89, "T": 0.03},
        {"A": 0.03, "C": 0.05, "G": 0.89, "T": 0.03},
        {"A": 0.03, "C": 0.88, "G": 0.06, "T": 0.03},
        {"A": 0.03, "C": 0.05, "G": 0.89, "T": 0.03},
        {"A": 0.03, "C": 0.05, "G": 0.89, "T": 0.03},
    ]),
    "CTCF": _build_numpy_pwm([
        {"A": 0.02, "C": 0.88, "G": 0.05, "T": 0.05},
        {"A": 0.02, "C": 0.89, "G": 0.05, "T": 0.04},
        {"A": 0.02, "C": 0.05, "G": 0.88, "T": 0.05},
        {"A": 0.02, "C": 0.87, "G": 0.06, "T": 0.05},
        {"A": 0.02, "C": 0.05, "G": 0.86, "T": 0.07},
        {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
        {"A": 0.02, "C": 0.05, "G": 0.88, "T": 0.05},
        {"A": 0.02, "C": 0.05, "G": 0.88, "T": 0.05},
    ]),
    "TATA": _build_numpy_pwm([
        {"A": 0.04, "C": 0.04, "G": 0.04, "T": 0.88},
        {"A": 0.88, "C": 0.04, "G": 0.04, "T": 0.04},
        {"A": 0.04, "C": 0.04, "G": 0.04, "T": 0.88},
        {"A": 0.88, "C": 0.04, "G": 0.04, "T": 0.04},
        {"A": 0.88, "C": 0.04, "G": 0.04, "T": 0.04},
        {"A": 0.88, "C": 0.04, "G": 0.04, "T": 0.04},
    ]),
    "EBOX": _build_numpy_pwm([
        {"A": 0.04, "C": 0.88, "G": 0.04, "T": 0.04},
        {"A": 0.88, "C": 0.04, "G": 0.04, "T": 0.04},
        {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
        {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
        {"A": 0.04, "C": 0.04, "G": 0.04, "T": 0.88},
        {"A": 0.04, "C": 0.04, "G": 0.88, "T": 0.04},
    ]),
    "AP1": _build_numpy_pwm([
        {"A": 0.04, "C": 0.04, "G": 0.04, "T": 0.88},
        {"A": 0.04, "C": 0.04, "G": 0.88, "T": 0.04},
        {"A": 0.88, "C": 0.04, "G": 0.04, "T": 0.04},
        {"A": 0.04, "C": 0.88, "G": 0.04, "T": 0.04},
        {"A": 0.04, "C": 0.04, "G": 0.04, "T": 0.88},
        {"A": 0.04, "C": 0.88, "G": 0.04, "T": 0.04},
        {"A": 0.88, "C": 0.04, "G": 0.04, "T": 0.04},
    ]),
}

_RC_MAP = str.maketrans("ACGTacgt", "TGCAtgca")


def _revcomp(seq):
    return seq.translate(_RC_MAP)[::-1]


def _score_sequence_against_numpy_pwm(seq, pwm_matrix):
    """Score a plain DNA string against a 4×L numpy log-odds PWM.

    Returns the raw log-odds sum at every valid position as a 1-D array.
    Positions containing non-ACGT characters contribute 0 (neutral).
    """
    L = pwm_matrix.shape[1]
    n = len(seq)
    if n < L:
        return np.array([])

    scores = np.zeros(n - L + 1, dtype=float)
    for i in range(n - L + 1):
        s = 0.0
        for k in range(L):
            base = seq[i + k].upper()
            row = _BASE_IDX.get(base, -1)
            if row >= 0:
                s += pwm_matrix[row, k]
        scores[i] = s
    return scores


def scan_pwm(sequence, pwm_matrix):
    """Score every L-mer in *sequence* against a 4×L numpy log-odds PWM.

    The PWM is scanned on both strands (forward and reverse complement).
    Each position score is the sum of position-specific log-odds values.

    Parameters
    ----------
    sequence   : str — DNA sequence (any case; non-ACGT treated as neutral)
    pwm_matrix : numpy.ndarray of shape (4, L) — rows are A, C, G, T

    Returns
    -------
    (best_score, best_position, best_strand) where
        best_score    : float — highest log-odds score found
        best_position : int  — 0-based start of the best-scoring L-mer
        best_strand   : str  — '+' or '-'
    """
    if not sequence or pwm_matrix is None or pwm_matrix.size == 0:
        return (float("-inf"), -1, "+")

    fwd_scores = _score_sequence_against_numpy_pwm(sequence, pwm_matrix)
    rc_seq = _revcomp(sequence)
    rev_scores = _score_sequence_against_numpy_pwm(rc_seq, pwm_matrix)

    best_score = float("-inf")
    best_pos = -1
    best_strand = "+"

    if fwd_scores.size > 0:
        idx = int(np.argmax(fwd_scores))
        if fwd_scores[idx] > best_score:
            best_score = float(fwd_scores[idx])
            best_pos = idx
            best_strand = "+"

    if rev_scores.size > 0:
        idx = int(np.argmax(rev_scores))
        # Convert reverse-complement position back to forward-strand coordinates.
        fwd_equiv = len(sequence) - pwm_matrix.shape[1] - idx
        if rev_scores[idx] > best_score:
            best_score = float(rev_scores[idx])
            best_pos = fwd_equiv
            best_strand = "-"

    return (best_score, best_pos, best_strand)


PWM_LIBRARY = {
    "CTCF": {
        "consensus": "CCGCGNGGNGGCAG",
        "matrix": [
            {"A": 0.02, "C": 0.88, "G": 0.05, "T": 0.05},
            {"A": 0.02, "C": 0.89, "G": 0.05, "T": 0.04},
            {"A": 0.02, "C": 0.05, "G": 0.88, "T": 0.05},
            {"A": 0.02, "C": 0.87, "G": 0.06, "T": 0.05},
            {"A": 0.02, "C": 0.05, "G": 0.86, "T": 0.07},
            {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
            {"A": 0.02, "C": 0.05, "G": 0.88, "T": 0.05},
            {"A": 0.02, "C": 0.05, "G": 0.88, "T": 0.05},
        ],
    },
    "SP1": {
        "consensus": "GGGCGG",
        "matrix": [
            {"A": 0.03, "C": 0.05, "G": 0.89, "T": 0.03},
            {"A": 0.03, "C": 0.05, "G": 0.89, "T": 0.03},
            {"A": 0.03, "C": 0.05, "G": 0.89, "T": 0.03},
            {"A": 0.03, "C": 0.88, "G": 0.06, "T": 0.03},
            {"A": 0.03, "C": 0.05, "G": 0.89, "T": 0.03},
            {"A": 0.03, "C": 0.05, "G": 0.89, "T": 0.03},
        ],
    },
    "NF1": {
        "consensus": "TTGGCN",
        "matrix": [
            {"A": 0.04, "C": 0.04, "G": 0.04, "T": 0.88},
            {"A": 0.04, "C": 0.04, "G": 0.04, "T": 0.88},
            {"A": 0.04, "C": 0.05, "G": 0.87, "T": 0.04},
            {"A": 0.04, "C": 0.05, "G": 0.87, "T": 0.04},
            {"A": 0.04, "C": 0.88, "G": 0.04, "T": 0.04},
            {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
        ],
    },
    "TATA": {
        "consensus": "TATAAA",
        "matrix": [
            {"A": 0.04, "C": 0.04, "G": 0.04, "T": 0.88},
            {"A": 0.88, "C": 0.04, "G": 0.04, "T": 0.04},
            {"A": 0.04, "C": 0.04, "G": 0.04, "T": 0.88},
            {"A": 0.88, "C": 0.04, "G": 0.04, "T": 0.04},
            {"A": 0.88, "C": 0.04, "G": 0.04, "T": 0.04},
            {"A": 0.88, "C": 0.04, "G": 0.04, "T": 0.04},
        ],
    },
    "EBOX": {
        "consensus": "CANNTG",
        "matrix": [
            {"A": 0.04, "C": 0.88, "G": 0.04, "T": 0.04},
            {"A": 0.88, "C": 0.04, "G": 0.04, "T": 0.04},
            {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
            {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
            {"A": 0.04, "C": 0.04, "G": 0.04, "T": 0.88},
            {"A": 0.04, "C": 0.04, "G": 0.88, "T": 0.04},
        ],
    },
    "AP1": {
        "consensus": "TGACTCA",
        "matrix": [
            {"A": 0.04, "C": 0.04, "G": 0.04, "T": 0.88},
            {"A": 0.04, "C": 0.04, "G": 0.88, "T": 0.04},
            {"A": 0.88, "C": 0.04, "G": 0.04, "T": 0.04},
            {"A": 0.04, "C": 0.88, "G": 0.04, "T": 0.04},
            {"A": 0.04, "C": 0.04, "G": 0.04, "T": 0.88},
            {"A": 0.04, "C": 0.88, "G": 0.04, "T": 0.04},
            {"A": 0.88, "C": 0.04, "G": 0.04, "T": 0.04},
        ],
    },
    "GATA": {
        "consensus": "WGATAR",
        "matrix": [
            {"A": 0.44, "C": 0.06, "G": 0.06, "T": 0.44},
            {"A": 0.04, "C": 0.04, "G": 0.88, "T": 0.04},
            {"A": 0.88, "C": 0.04, "G": 0.04, "T": 0.04},
            {"A": 0.04, "C": 0.04, "G": 0.04, "T": 0.88},
            {"A": 0.88, "C": 0.04, "G": 0.04, "T": 0.04},
            {"A": 0.44, "C": 0.06, "G": 0.44, "T": 0.06},
        ],
    },
    "ETS": {
        "consensus": "CCGGAA",
        "matrix": [
            {"A": 0.04, "C": 0.88, "G": 0.04, "T": 0.04},
            {"A": 0.04, "C": 0.88, "G": 0.04, "T": 0.04},
            {"A": 0.04, "C": 0.04, "G": 0.88, "T": 0.04},
            {"A": 0.04, "C": 0.04, "G": 0.88, "T": 0.04},
            {"A": 0.88, "C": 0.04, "G": 0.04, "T": 0.04},
            {"A": 0.88, "C": 0.04, "G": 0.04, "T": 0.04},
        ],
    },
}


def _score_position(seq, pos, matrix):
    score = 0.0
    for i, col in enumerate(matrix):
        if pos + i >= len(seq):
            return -999.0
        base = seq[pos + i]
        prob = col.get(base, 0.01)
        score += math.log2(max(prob, 1e-10) / 0.25)
    return score


def _max_min_score(matrix):
    max_s = sum(math.log2(max(col.values()) / 0.25) for col in matrix)
    min_s = sum(math.log2(min(max(v, 1e-10) for v in col.values()) / 0.25) for col in matrix)
    return max_s, min_s


def scan_pwms(seq, threshold=0.7):
    if not seq:
        return {}
    seq_upper = seq.upper()
    results = {}
    for motif_name, motif_data in PWM_LIBRARY.items():
        matrix = motif_data["matrix"]
        motif_len = len(matrix)
        max_s, min_s = _max_min_score(matrix)
        best_norm = 0.0
        best_pos = -1
        for pos in range(len(seq_upper) - motif_len + 1):
            raw = _score_position(seq_upper, pos, matrix)
            span = max_s - min_s
            if span == 0:
                norm = 0.0
            else:
                norm = (raw - min_s) / span
            norm = max(0.0, min(1.0, norm))
            if norm > best_norm:
                best_norm = norm
                best_pos = pos
        results[motif_name] = {"score": best_norm, "pos": best_pos, "found": best_norm >= threshold}
    return results
