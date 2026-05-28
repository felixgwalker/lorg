"""K-mer indexing and seed finding for synteny detection.

build_kmer_index: builds a {kmer: (chrom, pos)} dict from genome1 where
  each entry holds exactly one position (non-unique k-mers are stored as None).

find_unique_seeds: scans genome2 for k-mers (and their reverse complements)
  that map uniquely to genome1, returning (chrom1, pos1, chrom2, pos2, strand)
  tuples.  strand "+" means forward match; "-" means k-mer in genome2 matches
  the reverse complement of a genome1 k-mer (inversion evidence).
"""

_COMP = str.maketrans("ACGT", "TGCA")


def _rc(seq):
    return seq.translate(_COMP)[::-1]


def build_kmer_index(genome, k=15):
    """Build a k-mer index from genome.

    Returns {kmer: (chrom, pos)} for k-mers appearing exactly once,
    and {kmer: None} for k-mers appearing more than once (non-unique).
    Non-unique entries are kept to allow O(1) uniqueness checks without
    storing all positions.
    """
    index = {}
    for chrom, seq in genome.items():
        n = len(seq)
        for pos in range(n - k + 1):
            kmer = seq[pos:pos + k]
            if kmer in index:
                if index[kmer] is not None:
                    index[kmer] = None   # mark non-unique
            else:
                index[kmer] = (chrom, pos)
    return index


def find_unique_seeds(genome2, index1, k=15):
    """Find unique anchor seeds between genome2 and a pre-built genome1 index.

    Scans every k-mer in genome2 (forward and reverse complement) and emits a
    seed when the k-mer (or its RC) appears exactly once in genome1 and exactly
    once in genome2 on the same strand.

    Seeds: (chrom1, pos1, chrom2, pos2, strand)
      strand "+"  — forward match
      strand "-"  — RC match (inversion evidence)
    """
    # First pass: count occurrences of each forward and RC k-mer in genome2
    # Store as {kmer: (chrom, pos)} if unique, None if multi-hit
    fwd_count = {}   # forward k-mers in genome2
    rev_count = {}   # RC k-mers in genome2

    for chrom2, seq in genome2.items():
        n = len(seq)
        for pos2 in range(n - k + 1):
            kmer = seq[pos2:pos2 + k]
            # Forward
            if kmer in fwd_count:
                if fwd_count[kmer] is not None:
                    fwd_count[kmer] = None
            else:
                fwd_count[kmer] = (chrom2, pos2)
            # Reverse complement
            rc = _rc(kmer)
            if rc in rev_count:
                if rev_count[rc] is not None:
                    rev_count[rc] = None
            else:
                rev_count[rc] = (chrom2, pos2)

    seeds = []

    # Forward seeds: kmer unique in genome2 (fwd) AND unique in genome1
    for kmer, hit2 in fwd_count.items():
        if hit2 is None:
            continue
        hit1 = index1.get(kmer)
        if hit1 is None:
            continue   # absent or non-unique in genome1
        chrom1, pos1 = hit1
        chrom2, pos2 = hit2
        seeds.append((chrom1, pos1, chrom2, pos2, "+"))

    # Reverse seeds: RC of kmer unique in genome2 AND the RC is unique in genome1
    # This finds inversion evidence: a region in genome2 that, read backwards,
    # matches a unique forward k-mer in genome1.
    for rc_kmer, hit2 in rev_count.items():
        if hit2 is None:
            continue
        hit1 = index1.get(rc_kmer)
        if hit1 is None:
            continue   # absent or non-unique in genome1
        # Avoid double-emitting seeds already covered as forward matches
        if fwd_count.get(rc_kmer) is not None:
            # This RC kmer also appears forward in genome2 at same chrom2 — skip
            # to avoid ambiguity (genuine forward seed already emitted)
            fwd_hit = fwd_count[rc_kmer]
            if fwd_hit is not None and fwd_hit[0] == hit2[0]:
                continue
        chrom1, pos1 = hit1
        chrom2, pos2 = hit2
        seeds.append((chrom1, pos1, chrom2, pos2, "-"))

    return seeds
