import difflib
from collections import namedtuple
import numpy as np


# ---------------------------------------------------------------------------
# Public namedtuple returned by align_sequences()
# ---------------------------------------------------------------------------
Alignment = namedtuple("Alignment", ["aligned_seq1", "aligned_seq2", "cigar"])


# ---------------------------------------------------------------------------
# CIGAR helpers
# ---------------------------------------------------------------------------

def _build_cigar(aligned1: str, aligned2: str) -> str:
    """Build a CIGAR string from two gapped alignment strings.

    Operations used:
        M — match or mismatch (both columns are bases)
        I — insertion in seq2 relative to seq1 (seq1 has a gap '-')
        D — deletion in seq2 relative to seq1 (seq2 has a gap '-')
    """
    if not aligned1:
        return ""
    ops = []
    for a, b in zip(aligned1, aligned2):
        if a == "-":
            ops.append("I")
        elif b == "-":
            ops.append("D")
        else:
            ops.append("M")

    # Run-length encode.
    cigar_parts = []
    current_op = ops[0]
    count = 1
    for op in ops[1:]:
        if op == current_op:
            count += 1
        else:
            cigar_parts.append(f"{count}{current_op}")
            current_op = op
            count = 1
    cigar_parts.append(f"{count}{current_op}")
    return "".join(cigar_parts)


# ---------------------------------------------------------------------------
# difflib-based pairwise aligner (fast, pure Python)
# ---------------------------------------------------------------------------

def _difflib_align(seq1: str, seq2: str):
    """Align two sequences using difflib.SequenceMatcher.

    Returns (aligned_seq1, aligned_seq2) with '-' gap characters.
    """
    sm = difflib.SequenceMatcher(None, seq1, seq2, autojunk=False)
    al1: list[str] = []
    al2: list[str] = []
    for tag, i1, i2, j1, j2 in sm.get_opcodes():
        s1 = seq1[i1:i2]
        s2 = seq2[j1:j2]
        if tag == "equal" or tag == "replace":
            # Pad the shorter side with gaps so lengths match.
            len1, len2 = len(s1), len(s2)
            max_len = max(len1, len2)
            al1.append(s1 + "-" * (max_len - len1))
            al2.append(s2 + "-" * (max_len - len2))
        elif tag == "insert":
            # seq2 has extra bases; seq1 gets gaps.
            al1.append("-" * len(s2))
            al2.append(s2)
        elif tag == "delete":
            # seq1 has extra bases; seq2 gets gaps.
            al1.append(s1)
            al2.append("-" * len(s1))
    return "".join(al1), "".join(al2)


# ---------------------------------------------------------------------------
# Public API: align_sequences()
# ---------------------------------------------------------------------------

CHUNK_SIZE = 100_000   # 100 kb windows for large-sequence chunking
CHUNK_OVERLAP = 500    # overlap to avoid missing variants at boundaries


def align_sequences(seq1: str, seq2: str) -> "Alignment":
    """Align two sequences and return an Alignment namedtuple.

    For sequences up to CHUNK_SIZE bp the alignment is computed directly.
    Longer sequences (up to ~10 Mb) are chunked into CHUNK_SIZE windows
    with CHUNK_OVERLAP overlap; alignments are concatenated (the overlapping
    tails from all but the last chunk are trimmed).

    Parameters
    ----------
    seq1 : str
        Proxy (reference) sequence — plain nucleotide string.
    seq2 : str
        Target (query) sequence — plain nucleotide string.

    Returns
    -------
    Alignment
        namedtuple with fields aligned_seq1, aligned_seq2, cigar.
    """
    # Short sequences: align directly.
    if len(seq1) <= CHUNK_SIZE and len(seq2) <= CHUNK_SIZE:
        a1, a2 = _difflib_align(seq1, seq2)
        return Alignment(aligned_seq1=a1, aligned_seq2=a2, cigar=_build_cigar(a1, a2))

    # Long sequences: chunk seq1 and find corresponding region in seq2.
    parts1: list[str] = []
    parts2: list[str] = []

    pos1 = 0
    pos2 = 0
    while pos1 < len(seq1):
        end1 = min(pos1 + CHUNK_SIZE, len(seq1))
        # Estimate the corresponding position in seq2 by proportional scaling.
        scale = len(seq2) / max(len(seq1), 1)
        est_end2 = min(int(end1 * scale) + CHUNK_OVERLAP, len(seq2))
        est_start2 = max(int(pos1 * scale) - CHUNK_OVERLAP, 0)

        chunk1 = seq1[pos1:end1]
        chunk2 = seq2[est_start2:est_end2]

        ca1, ca2 = _difflib_align(chunk1, chunk2)

        # Trim overlap from the end of all but the last chunk to avoid
        # double-counting boundary variants.
        if end1 < len(seq1):
            trim = CHUNK_OVERLAP // 2
            # Only trim if we have enough aligned columns.
            if len(ca1) > trim:
                ca1 = ca1[:-trim]
                ca2 = ca2[:-trim]

        parts1.append(ca1)
        parts2.append(ca2)

        pos1 = end1
        pos2 = est_end2

    aligned1 = "".join(parts1)
    aligned2 = "".join(parts2)
    return Alignment(
        aligned_seq1=aligned1,
        aligned_seq2=aligned2,
        cigar=_build_cigar(aligned1, aligned2),
    )


MATCH = 1
MISMATCH = -1
GAP_OPEN = -2


def needleman_wunsch(seq1, seq2):
    n = len(seq1)
    m = len(seq2)
    score_mat = np.zeros((n + 1, m + 1), dtype=np.int32)
    trace_mat = np.zeros((n + 1, m + 1), dtype=np.int8)

    for i in range(1, n + 1):
        score_mat[i, 0] = GAP_OPEN * i
        trace_mat[i, 0] = 1
    for j in range(1, m + 1):
        score_mat[0, j] = GAP_OPEN * j
        trace_mat[0, j] = 2

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match_score = MATCH if seq1[i-1] == seq2[j-1] else MISMATCH
            diag = score_mat[i-1, j-1] + match_score
            up = score_mat[i-1, j] + GAP_OPEN
            left = score_mat[i, j-1] + GAP_OPEN
            best = max(diag, up, left)
            score_mat[i, j] = best
            if best == diag:
                trace_mat[i, j] = 0
            elif best == up:
                trace_mat[i, j] = 1
            else:
                trace_mat[i, j] = 2

    aligned1 = []
    aligned2 = []
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0 and trace_mat[i, j] == 0:
            aligned1.append(seq1[i-1])
            aligned2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif i > 0 and trace_mat[i, j] == 1:
            aligned1.append(seq1[i-1])
            aligned2.append("-")
            i -= 1
        else:
            aligned1.append("-")
            aligned2.append(seq2[j-1])
            j -= 1

    return "".join(reversed(aligned1)), "".join(reversed(aligned2))


def call_variants_from_alignment(aligned_proxy, aligned_target, chrom, window_start):
    variants = []
    proxy_pos = window_start
    target_pos = window_start

    i = 0
    n = len(aligned_proxy)
    while i < n:
        pb = aligned_proxy[i]
        tb = aligned_target[i]

        if pb == "-":
            ins_seq = []
            j = i
            while j < n and aligned_proxy[j] == "-":
                ins_seq.append(aligned_target[j])
                j += 1
            variants.append({
                "chrom": chrom,
                "pos": proxy_pos,
                "ref": ".",
                "alt": "".join(ins_seq),
                "type": "INS",
                "length": len(ins_seq),
            })
            target_pos += len(ins_seq)
            i = j
            continue
        elif tb == "-":
            del_seq = []
            j = i
            while j < n and aligned_target[j] == "-":
                del_seq.append(aligned_proxy[j])
                j += 1
            variants.append({
                "chrom": chrom,
                "pos": proxy_pos,
                "ref": "".join(del_seq),
                "alt": ".",
                "type": "DEL",
                "length": len(del_seq),
            })
            proxy_pos += len(del_seq)
            i = j
            continue
        else:
            if pb != tb:
                variants.append({
                    "chrom": chrom,
                    "pos": proxy_pos,
                    "ref": pb,
                    "alt": tb,
                    "type": "SNV",
                    "length": 1,
                })
            proxy_pos += 1
            target_pos += 1
            i += 1

    return variants


def align_genomes(proxy_seqs, target_seqs, window_size=500, n_samples=1000, rng=None):
    if rng is None:
        rng = np.random.default_rng(42)

    all_variants = []
    chrom_lengths = {}

    for chrom in proxy_seqs:
        if chrom not in target_seqs:
            continue
        pseq = proxy_seqs[chrom]
        tseq = target_seqs[chrom]
        chrom_lengths[chrom] = len(pseq)

        max_start = max(len(pseq) - window_size, 1)
        n_windows = min(n_samples // len(proxy_seqs), max_start)
        n_windows = max(n_windows, 1)

        if max_start <= n_windows:
            starts = list(range(0, max_start, max(1, max_start // n_windows)))
        else:
            starts = sorted(rng.integers(0, max_start, size=n_windows).tolist())

        sampled_fraction = (n_windows * window_size) / max(len(pseq), 1)
        sampled_fraction = min(sampled_fraction, 1.0)

        chrom_variants = []
        for start in starts:
            end_p = min(start + window_size, len(pseq))
            end_t = min(start + window_size, len(tseq))
            pw = pseq[start:end_p]
            tw = tseq[start:end_t]
            if not pw or not tw:
                continue
            ap, at = needleman_wunsch(pw, tw)
            window_vars = call_variants_from_alignment(ap, at, chrom, start)
            chrom_variants.extend(window_vars)

        if sampled_fraction < 1.0 and sampled_fraction > 0:
            scale = 1.0 / sampled_fraction
            scaled = []
            seen_pos = set()
            for v in chrom_variants:
                key = (v["chrom"], v["pos"], v["type"])
                if key not in seen_pos:
                    seen_pos.add(key)
                    scaled.append(v)
            all_variants.extend(scaled)
        else:
            all_variants.extend(chrom_variants)

    return all_variants, chrom_lengths
