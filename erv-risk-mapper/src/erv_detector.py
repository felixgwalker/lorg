"""K-mer-based ERV signature detection across genome windows."""

WINDOW_SIZE = 5000
STEP_SIZE = 2500
MIN_KMER_HITS = 3

# ERV family signatures: each family has a set of consensus k-mers
FAMILY_SIGNATURES: dict[str, list[str]] = {
    "HERV-K": ["TGAAAGAC", "GTCTTCA", "AATAAA", "CCCGGG", "CTGCAG"],
    "HERV-H": ["TGAAAGAC", "AATAAA", "TATTAG", "CCCGGG", "GCCCGG"],
    "HERV-W": ["TGAAAGAC", "GTCTTCA", "AATAAA", "GCCCGG", "TTTAAA"],
    "HERV-E": ["TGAAAGAC", "AATAAA", "CTGCAG", "ACATGT", "CCCGGG"],
    "HERV-L": ["TGAAAGAC", "AATAAA", "GGATCC", "GAGCTC", "GTCTTCA"],
}

LTR_5_MOTIF = "TGAAAGAC"
LTR_3_MOTIF = "GTCTTCA"
POLY_A_SIGNAL = "AATAAA"


def _count_kmer_hits(window: str, kmers: list[str]) -> int:
    """Count how many distinct kmers from the list appear in window."""
    return sum(1 for kmer in kmers if kmer in window)


def _check_ltr(window: str) -> str:
    """Determine LTR type based on motif presence."""
    has_5 = LTR_5_MOTIF in window
    has_3 = LTR_3_MOTIF in window
    if has_5 and has_3:
        return "full"
    elif has_5 or has_3:
        return "partial"
    return "solo"


def _find_longest_orf(window: str, min_aa: int = 500) -> int:
    """Return length of longest ORF in any frame (in bp). min_aa in amino acids."""
    stop_codons = {"TAA", "TAG", "TGA"}
    max_orf_bp = 0
    for frame in range(3):
        orf_start = None
        i = frame
        while i + 2 < len(window):
            codon = window[i:i + 3]
            if codon == "ATG" and orf_start is None:
                orf_start = i
            elif codon in stop_codons and orf_start is not None:
                orf_len = i - orf_start
                max_orf_bp = max(max_orf_bp, orf_len)
                orf_start = None
            i += 3
    return max_orf_bp


def detect_erv_windows(genome: dict[str, str]) -> list[dict]:
    """Scan genome in windows and return ERV hit records."""
    hits = []
    for chrom, seq in genome.items():
        n = len(seq)
        pos = 0
        while pos < n:
            end = min(pos + WINDOW_SIZE, n)
            window = seq[pos:end]

            best_family = None
            best_hits = 0
            for family, kmers in FAMILY_SIGNATURES.items():
                kmer_hits = _count_kmer_hits(window, kmers)
                if kmer_hits >= MIN_KMER_HITS and kmer_hits > best_hits:
                    best_hits = kmer_hits
                    best_family = family

            if best_family is not None:
                ltr_type = _check_ltr(window)
                has_poly_a = POLY_A_SIGNAL in window
                longest_orf_bp = _find_longest_orf(window)
                gc = (window.count("G") + window.count("C")) / max(len(window), 1)
                hits.append({
                    "chrom": chrom,
                    "start": pos,
                    "end": end,
                    "family": best_family,
                    "kmer_hits": best_hits,
                    "ltr_type": ltr_type,
                    "has_poly_a": has_poly_a,
                    "longest_orf_bp": longest_orf_bp,
                    "gc_content": round(gc, 4),
                    "window_len": end - pos,
                })

            pos += STEP_SIZE

    return hits
