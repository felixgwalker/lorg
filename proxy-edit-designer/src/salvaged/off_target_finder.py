"""Find off-target sites for guide RNAs using k-mer index + mismatch counting.

Salvaged from guide-rna-off-target-scorer (deleted stage1f).
"""

from .genome_indexer import KMER_LEN, build_genome_index, reverse_complement

MAX_MISMATCHES = 3
SEED_LEN = 12


def _count_mismatches(seq_a: str, seq_b: str) -> int:
    return sum(a != b for a, b in zip(seq_a, seq_b))


def _mismatch_positions(guide: str, target: str) -> list[int]:
    return [i + 1 for i, (a, b) in enumerate(zip(guide, target)) if a != b]


def _get_pam_fwd(chrom_seq: str, pos: int) -> str:
    pam_start = pos + KMER_LEN
    if pam_start + 2 <= len(chrom_seq):
        return chrom_seq[pam_start:pam_start + 2]
    return "NN"


def _get_pam_rev(chrom_seq: str, pos: int) -> str:
    pam_start = pos - 2
    if pam_start >= 0:
        raw = chrom_seq[pam_start:pam_start + 2]
        if len(raw) == 2:
            return reverse_complement(raw)
    return "NN"


def _is_ngg_pam(pam: str) -> bool:
    return len(pam) == 2 and pam[0] == "G" and pam[1] == "G"


def find_off_targets(
    guide_name: str,
    guide_seq: str,
    genome: dict[str, str],
    max_mismatches: int = MAX_MISMATCHES,
) -> list[dict]:
    """Return off-target hit dicts for a single guide via seed-and-extend."""
    guide_seq = guide_seq.upper()
    index = build_genome_index(genome)
    seen: set[tuple[str, int, str]] = set()
    results: list[dict] = []
    for seed_offset in range(KMER_LEN - SEED_LEN + 1):
        seed = guide_seq[seed_offset:seed_offset + SEED_LEN]
        if seed not in index:
            continue
        for chrom, hit_pos, strand in index[seed]:
            candidate_start = hit_pos - seed_offset
            loc_key = (chrom, candidate_start, strand)
            if loc_key in seen:
                continue
            seen.add(loc_key)
            if candidate_start < 0:
                continue
            chrom_seq = genome[chrom].upper()
            if candidate_start + KMER_LEN > len(chrom_seq):
                continue
            if strand == "+":
                target_seq = chrom_seq[candidate_start:candidate_start + KMER_LEN]
            else:
                target_seq = reverse_complement(
                    chrom_seq[candidate_start:candidate_start + KMER_LEN]
                )
            if "N" in target_seq:
                continue
            mm = _count_mismatches(guide_seq, target_seq)
            if mm > max_mismatches:
                continue
            pam = _get_pam_fwd(chrom_seq, candidate_start) if strand == "+" else _get_pam_rev(chrom_seq, candidate_start)
            results.append({
                "guide_name": guide_name,
                "guide_seq": guide_seq,
                "chrom": chrom,
                "pos": candidate_start,
                "strand": strand,
                "target_seq": target_seq,
                "mismatches": mm,
                "mismatch_positions": _mismatch_positions(guide_seq, target_seq),
                "pam": pam,
                "has_pam": _is_ngg_pam(pam),
            })
    return results
