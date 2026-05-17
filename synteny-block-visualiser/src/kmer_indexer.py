def build_kmer_index(genome, k=15):
    index = {}
    for chrom, seq in genome.items():
        for pos in range(len(seq) - k + 1):
            kmer = seq[pos:pos + k]
            if kmer not in index:
                index[kmer] = []
            index[kmer].append((chrom, pos))
    return index


def find_unique_seeds(genome2, index1, k=15):
    seeds = []
    seen_kmers = set()
    for chrom2, seq in genome2.items():
        for pos2 in range(len(seq) - k + 1):
            kmer = seq[pos2:pos2 + k]
            if kmer in seen_kmers:
                continue
            hits1 = index1.get(kmer, [])
            if len(hits1) == 1:
                chrom1, pos1 = hits1[0]
                seeds.append((chrom1, pos1, chrom2, pos2, kmer))
                seen_kmers.add(kmer)
    return seeds
