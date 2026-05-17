def chain_seeds(seeds, min_chain_score=3, gap_penalty_divisor=1000):
    if not seeds:
        return []

    chrom_pairs = {}
    for seed in seeds:
        chrom1, pos1, chrom2, pos2, kmer = seed
        key = (chrom1, chrom2)
        if key not in chrom_pairs:
            chrom_pairs[key] = []
        chrom_pairs[key].append((pos1, pos2, kmer))

    chains = []
    for (chrom1, chrom2), pair_seeds in chrom_pairs.items():
        pair_seeds.sort(key=lambda s: (s[0], s[1]))
        n = len(pair_seeds)
        if n == 0:
            continue

        score = [1] * n
        prev = [-1] * n

        for i in range(1, n):
            for j in range(i):
                p1_j, p2_j, _ = pair_seeds[j]
                p1_i, p2_i, _ = pair_seeds[i]
                g1_gap = p1_i - p1_j
                g2_gap = p2_i - p2_j
                if g1_gap > 0 and g2_gap > 0:
                    penalty = max(g1_gap, g2_gap) / gap_penalty_divisor
                    new_score = score[j] + 1 - penalty
                    if new_score > score[i]:
                        score[i] = new_score
                        prev[i] = j

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
            chain_seeds_list = [pair_seeds[idx] for idx in chain_indices]
            chains.append({
                "chrom1": chrom1,
                "chrom2": chrom2,
                "seeds": chain_seeds_list,
                "n_seeds": len(chain_seeds_list),
                "g1_start": chain_seeds_list[0][0],
                "g1_end": chain_seeds_list[-1][0],
                "g2_start": chain_seeds_list[0][1],
                "g2_end": chain_seeds_list[-1][1],
            })

    return chains


def chain_seeds_with_inversions(seeds, min_chain_score=3, gap_penalty_divisor=1000, k=15):
    forward_seeds = []
    reverse_seeds = []

    comp = str.maketrans("ACGT", "TGCA")

    kmer_set = {s[4] for s in seeds}
    rc_map = {km: km.translate(comp)[::-1] for km in kmer_set}

    seed_by_kmer1 = {}
    for s in seeds:
        chrom1, pos1, chrom2, pos2, kmer = s
        seed_by_kmer1[kmer] = s

    for s in seeds:
        chrom1, pos1, chrom2, pos2, kmer = s
        rc = rc_map.get(kmer, "")
        if rc in seed_by_kmer1:
            rev_s = seed_by_kmer1[rc]
            reverse_seeds.append((chrom1, pos1, chrom2, len_placeholder := pos2, kmer, "rev"))
        forward_seeds.append(s)

    fwd_chains = chain_seeds(seeds, min_chain_score, gap_penalty_divisor)
    return fwd_chains
