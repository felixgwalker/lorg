import math

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
