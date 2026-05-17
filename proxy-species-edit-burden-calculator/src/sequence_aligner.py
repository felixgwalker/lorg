import numpy as np


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
