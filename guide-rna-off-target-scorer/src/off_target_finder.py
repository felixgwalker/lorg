"""Find off-target sites for guide RNAs using k-mer index + mismatch counting.

Strategy: seed-and-extend using the first 12 bases of the guide (PAM-proximal
seed region) as lookup keys into the genome k-mer index, then verify the full
20-mer against the candidate site.
"""

from .genome_indexer import KMER_LEN, build_genome_index, reverse_complement

MAX_MISMATCHES = 3
SEED_LEN = 12  # PAM-proximal seed (last 12 bases of 20-mer guide)


def _count_mismatches(seq_a: str, seq_b: str) -> int:
    return sum(a != b for a, b in zip(seq_a, seq_b))


def _mismatch_positions(guide: str, target: str) -> list[int]:
    """Return 1-based positions (PAM-proximal = 20) where mismatches occur."""
    return [i + 1 for i, (a, b) in enumerate(zip(guide, target)) if a != b]


def _get_pam_fwd(chrom_seq: str, pos: int) -> str:
    """Return 2-base NGG PAM immediately 3' of protospacer on fwd strand."""
    pam_start = pos + KMER_LEN
    if pam_start + 2 <= len(chrom_seq):
        return chrom_seq[pam_start:pam_start + 2]
    return "NN"


def _get_pam_rev(chrom_seq: str, pos: int) -> str:
    """Return 2-base PAM for reverse-strand hit at genomic pos.

    On the reverse strand the protospacer runs right-to-left; its PAM is
    immediately 5' of pos on the forward strand, i.e. bases [pos-2:pos],
    whose reverse-complement gives the PAM as read 5'->3' on the guide strand.
    """
    pam_start = pos - 2
    if pam_start >= 0:
        raw = chrom_seq[pam_start:pam_start + 2]
        if len(raw) == 2:
            return reverse_complement(raw)
    return "NN"


def _is_ngg_pam(pam: str) -> bool:
    """NGG PAM: second and third bases must both be G (pam is 2-char GG part)."""
    return len(pam) == 2 and pam[0] == "G" and pam[1] == "G"


def find_off_targets(
    guide_name: str,
    guide_seq: str,
    genome: dict[str, str],
    max_mismatches: int = MAX_MISMATCHES,
) -> list[dict]:
    """Return a list of off-target hit dicts for a single guide.

    Uses seed-and-extend: queries the genome k-mer index with every possible
    seed window of length SEED_LEN drawn from the guide, then extends each
    candidate to the full 20-mer to count total mismatches.

    Each returned dict has keys:
        guide_name, guide_seq, chrom, pos, strand,
        target_seq, mismatches, mismatch_positions, pam, has_pam
    """
    guide_seq = guide_seq.upper()
    index = build_genome_index(genome)

    # Collect candidate (chrom, pos, strand) locations via seed-and-extend.
    # We generate all sub-strings of the guide of length SEED_LEN and look
    # each one up in the index.  Any hit whose surrounding 20-mer has
    # <= max_mismatches vs the guide is kept.
    seen: set[tuple[str, int, str]] = set()
    results: list[dict] = []

    # Slide a SEED_LEN window across the full guide
    for seed_offset in range(KMER_LEN - SEED_LEN + 1):
        seed = guide_seq[seed_offset:seed_offset + SEED_LEN]
        if seed not in index:
            continue
        for chrom, hit_pos, strand in index[seed]:
            # Reconstruct the start of the 20-mer that contains this seed
            candidate_start = hit_pos - seed_offset
            loc_key = (chrom, candidate_start, strand)
            if loc_key in seen:
                continue
            seen.add(loc_key)

            # Bounds check
            if candidate_start < 0:
                continue
            chrom_seq = genome[chrom].upper()
            if candidate_start + KMER_LEN > len(chrom_seq):
                continue

            # Extract target 20-mer on the appropriate strand
            if strand == "+":
                target_seq = chrom_seq[candidate_start:candidate_start + KMER_LEN]
            else:
                # For reverse-strand hits the index stores the forward genomic
                # position of the reverse-complement k-mer.  We need to
                # extract the RC of the genomic sequence at that position.
                target_seq = reverse_complement(
                    chrom_seq[candidate_start:candidate_start + KMER_LEN]
                )

            if "N" in target_seq:
                continue

            mm = _count_mismatches(guide_seq, target_seq)
            if mm > max_mismatches:
                continue

            # PAM
            if strand == "+":
                pam = _get_pam_fwd(chrom_seq, candidate_start)
            else:
                pam = _get_pam_rev(chrom_seq, candidate_start)

            has_pam = _is_ngg_pam(pam)
            mm_positions = _mismatch_positions(guide_seq, target_seq)

            results.append({
                "guide_name": guide_name,
                "guide_seq": guide_seq,
                "chrom": chrom,
                "pos": candidate_start,
                "strand": strand,
                "target_seq": target_seq,
                "mismatches": mm,
                "mismatch_positions": mm_positions,
                "pam": pam,
                "has_pam": has_pam,
            })

    return results
