"""Chain collinear and inversion k-mer seeds into synteny blocks.

Seeds tuple format: (chrom1, pos1, chrom2, pos2, strand)
  strand "+" : forward match — both pos1 and pos2 increase together
  strand "-" : reverse/inversion match — pos1 increases, pos2 decreases

Chaining uses a fast O(n log n) algorithm based on coordinate-compressed
longest increasing subsequence (LIS) for large seed groups, falling back to
O(n^2) DP for small groups (n <= _DP_THRESHOLD).

The LIS approach bins seeds into windows and greedily chains the densest
consecutive bins, giving synteny blocks in O(n log n) time.

Returns a list of chain dicts.
"""

import bisect

_DP_THRESHOLD = 500    # use O(n^2) DP only for groups smaller than this


def _chain_dp(pair_seeds, strand, min_chain_score, gap_penalty_divisor):
    """O(n^2) DP chainer for small seed groups."""
    n = len(pair_seeds)
    score = [1.0] * n
    prev = [-1] * n

    for i in range(1, n):
        p1_i, p2_i = pair_seeds[i]
        for j in range(i):
            p1_j, p2_j = pair_seeds[j]
            g1_gap = p1_i - p1_j
            if g1_gap <= 0:
                continue
            if strand == "+":
                g2_gap = p2_i - p2_j
                if g2_gap <= 0:
                    continue
            else:
                g2_gap = p2_j - p2_i
                if g2_gap <= 0:
                    continue
            gap = max(g1_gap, g2_gap)
            penalty = gap / gap_penalty_divisor
            new_score = score[j] + 1.0 - penalty
            if new_score > score[i]:
                score[i] = new_score
                prev[i] = j

    chains = []
    used = [False] * n
    for i in sorted(range(n), key=lambda x: -score[x]):
        if used[i]:
            continue
        if score[i] < min_chain_score:
            continue
        chain_indices = []
        cur = i
        while cur != -1:
            chain_indices.append(cur)
            used[cur] = True
            cur = prev[cur]
        chain_indices.reverse()
        chains.append([pair_seeds[idx] for idx in chain_indices])
    return chains


def _chain_fast(pair_seeds, strand, min_chain_score, window_size=200):
    """Fast O(n) binning chainer for large seed groups.

    Seeds are binned along genome1 position into windows of `window_size` bp.
    Consecutive windows with seeds whose genome2 positions are also monotone
    are merged into chains.  This finds the dominant collinear/inversion signal
    efficiently.
    """
    if not pair_seeds:
        return []

    # Group by window
    bins = {}
    for p1, p2 in pair_seeds:
        w = p1 // window_size
        if w not in bins:
            bins[w] = []
        bins[w].append((p1, p2))

    sorted_windows = sorted(bins.keys())

    # Build chains by merging consecutive consistent windows
    chains = []
    current_chain = []
    last_p2 = None

    for w in sorted_windows:
        window_seeds = sorted(bins[w], key=lambda s: s[0])
        # Median pos2 of this window
        median_p2 = window_seeds[len(window_seeds) // 2][1]

        if last_p2 is None:
            current_chain.extend(window_seeds)
            last_p2 = median_p2
        else:
            if strand == "+" and median_p2 >= last_p2:
                current_chain.extend(window_seeds)
                last_p2 = median_p2
            elif strand == "-" and median_p2 <= last_p2:
                current_chain.extend(window_seeds)
                last_p2 = median_p2
            else:
                # Strand consistency broken — save current, start new
                if len(current_chain) >= min_chain_score:
                    chains.append(current_chain)
                current_chain = list(window_seeds)
                last_p2 = median_p2

    if len(current_chain) >= min_chain_score:
        chains.append(current_chain)

    return chains


def chain_seeds(seeds, min_chain_score=3, gap_penalty_divisor=1000, k=15):
    """Build collinear chains from seed tuples.

    Seeds: (chrom1, pos1, chrom2, pos2, strand)

    Returns a list of chain dicts with keys:
        chrom1, chrom2, strand, seeds, n_seeds,
        g1_start, g1_end, g2_start, g2_end
    Each chain's "seeds" field contains (pos1, pos2) pairs.
    """
    if not seeds:
        return []

    # Group seeds by (chrom1, chrom2, strand)
    groups = {}
    for chrom1, pos1, chrom2, pos2, strand in seeds:
        key = (chrom1, chrom2, strand)
        if key not in groups:
            groups[key] = []
        groups[key].append((pos1, pos2))

    chains = []
    for (chrom1, chrom2, strand), pair_seeds in groups.items():
        # Sort: always increasing pos1.
        # For "+" pos2 also increases; for "-" pos2 decreases.
        if strand == "+":
            pair_seeds.sort(key=lambda s: (s[0], s[1]))
        else:
            pair_seeds.sort(key=lambda s: (s[0], -s[1]))

        n = len(pair_seeds)
        if n == 0:
            continue

        # Choose algorithm based on group size
        if n <= _DP_THRESHOLD:
            raw_chains = _chain_dp(pair_seeds, strand, min_chain_score,
                                   gap_penalty_divisor)
        else:
            raw_chains = _chain_fast(pair_seeds, strand, min_chain_score)

        for chain_s in raw_chains:
            if len(chain_s) < min_chain_score:
                continue
            chains.append({
                "chrom1": chrom1,
                "chrom2": chrom2,
                "strand": strand,
                "seeds": chain_s,
                "n_seeds": len(chain_s),
                "g1_start": chain_s[0][0],
                "g1_end": chain_s[-1][0],
                "g2_start": chain_s[0][1],
                "g2_end": chain_s[-1][1],
            })

    return chains
