"""Find off-target sites for guide RNAs using k-mer index + mismatch counting."""

from .genome_indexer import KMER_LEN, reverse_complement

MAX_MISMATCHES = 3


def _count_mismatches(seq_a: str, seq_b: str) -> int:
    return sum(a != b for a, b in zip(seq_a, seq_b))


def _get_pam(chrom_seq: str, pos: int, strand: str) -> str:
    """Return 2-base PAM (positions 20-21 relative to protospacer start)."""
    if strand == "+":
        pam_start = pos + KMER_LEN
        if pam_start + 2 <= len(chrom_seq):
            return chrom_seq[pam_start:pam_start + 2]
    else:
        pam_start = pos - 2
        if pam_start >= 0:
            rc = reverse_complement(chrom_seq[pam_start:pam_start + 2])
            return rc
    return "NN"


def _is_valid_pam(pam: str) -> bool:
    """Check NGG PAM (N=any, GG required)."""
    return len(pam) == 2 and pam[1] == "G" and pam[0] == "G"


def _get_candidates_by_scanning(
    guide: str,
    genome: dict[str, str],
    max_mismatches: int = MAX_MISMATCHES,
) -> list[dict]:
    """Scan genome for sites with <= max_mismatches vs guide (no PAM filter)."""
    hits = []
    guide_upper = guide.upper()

    for chrom, seq in genome.items():
        n = len(seq)
        rc_seq = reverse_complement(seq)

        for strand, scan_seq in (("+", seq), ("-", rc_seq)):
            for i in range(n - KMER_LEN + 1):
                target = scan_seq[i:i + KMER_LEN]
                if "N" in target:
                    continue
                mm = _count_mismatches(guide_upper, target)
                if mm <= max_mismatches:
                    if strand == "+":
                        genomic_pos = i
                        pam_seq = seq[i + KMER_LEN:i + KMER_LEN + 2] if i + KMER_LEN + 2 <= n else "NN"
                    else:
                        genomic_pos = n - i - KMER_LEN
                        pam_raw = seq[genomic_pos - 2:genomic_pos] if genomic_pos >= 2 else "NN"
                        pam_seq = reverse_complement(pam_raw) if len(pam_raw) == 2 else "NN"
                    has_pam = _is_valid_pam(pam_seq)
                    hits.append({
                        "chrom": chrom,
                        "pos": genomic_pos,
                        "strand": strand,
                        "target_seq": target,
                        "mismatches": mm,
                        "pam": pam_seq,
                        "has_pam": has_pam,
                    })
    return hits


def find_off_targets(
    guide_name: str,
    guide_seq: str,
    genome: dict[str, str],
    max_mismatches: int = MAX_MISMATCHES,
) -> list[dict]:
    """Return list of off-target dicts for a single guide."""
    raw_hits = _get_candidates_by_scanning(guide_seq, genome, max_mismatches)
    results = []
    for hit in raw_hits:
        hit["guide_name"] = guide_name
        hit["guide_seq"] = guide_seq
        results.append(hit)
    return results
