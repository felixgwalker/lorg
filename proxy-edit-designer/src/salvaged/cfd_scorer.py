"""CFD scoring for CRISPR off-target sites.

Implements the Cutting Frequency Determination (CFD) score from:
    Doench et al. 2016, Nature Biotechnology 34:184-191.
    doi:10.1038/nbt.3437

Salvaged from guide-rna-off-target-scorer (deleted stage1f).
"""

MM_SCORES: dict[str, list[float]] = {
    "rA:dA": [0.0]*20, "rA:dC": [0.5]*20, "rA:dG": [0.5]*20, "rA:dT": [1.0]*20,
    "rC:dA": [0.5]*20, "rC:dC": [0.0]*20, "rC:dG": [0.5]*20, "rC:dT": [0.5]*20,
    "rG:dA": [0.5]*20, "rG:dC": [0.5]*20, "rG:dG": [0.0]*20, "rG:dT": [0.5]*20,
    "rT:dA": [0.5]*20, "rT:dC": [0.5]*20, "rT:dG": [0.5]*20, "rT:dT": [0.0]*20,
}

PAM_SCORES: dict[str, float] = {
    "AGG": 0.857, "TGG": 0.619, "CGG": 0.555, "GGG": 0.508,
    "AAG": 0.131, "ATG": 0.106, "ACG": 0.070, "AGT": 0.023,
    "AGC": 0.023, "AGA": 0.023, "AAT": 0.023, "AAA": 0.023,
    "AAC": 0.023, "ATC": 0.023, "ATT": 0.023, "ATA": 0.023,
    "ACC": 0.023, "ACT": 0.023, "ACA": 0.023, "TAG": 0.131,
    "TTG": 0.119, "TCG": 0.060, "TAT": 0.023, "TAA": 0.023,
    "TAC": 0.023, "TTT": 0.023, "TTA": 0.023, "TTC": 0.023,
    "TCC": 0.023, "TCT": 0.023, "TCA": 0.023, "TGT": 0.023,
    "TGA": 0.023, "TGC": 0.023, "CAG": 0.131, "CTG": 0.119,
    "CCG": 0.060, "CAT": 0.023, "CAA": 0.023, "CAC": 0.023,
    "CTT": 0.023, "CTA": 0.023, "CTC": 0.023, "CCT": 0.023,
    "CCA": 0.023, "CCC": 0.023, "CGT": 0.023, "CGA": 0.023,
    "CGC": 0.023, "GAG": 0.131, "GTG": 0.119, "GCG": 0.060,
    "GAT": 0.023, "GAA": 0.023, "GAC": 0.023, "GTT": 0.023,
    "GTA": 0.023, "GTC": 0.023, "GCC": 0.023, "GCT": 0.023,
    "GCA": 0.023, "GGT": 0.023, "GGA": 0.023, "GGC": 0.023,
}

KMER_LEN = 20


def _pam_score(pam: str) -> float:
    if len(pam) == 3:
        return PAM_SCORES.get(pam.upper(), 0.023)
    if len(pam) == 2:
        if pam.upper() == "GG":
            ngg = [PAM_SCORES[k] for k in PAM_SCORES if k.endswith("GG")]
            return sum(ngg) / len(ngg) if ngg else 0.5
        return 0.023
    return 0.0


def cfd_score(guide: str, off_target_seq: str, pam: str = "GG") -> float:
    """Compute CFD score between a guide and an off-target sequence (0–1)."""
    guide = guide.upper()
    off_target_seq = off_target_seq.upper()
    if len(guide) != KMER_LEN or len(off_target_seq) != KMER_LEN:
        raise ValueError(
            f"Guide and target must each be {KMER_LEN} nt; "
            f"got guide={len(guide)}, target={len(off_target_seq)}"
        )
    score = 1.0
    for i in range(KMER_LEN):
        g, t = guide[i], off_target_seq[i]
        if g != t:
            key = f"r{g}:d{t}"
            mm_vals = MM_SCORES.get(key)
            score *= mm_vals[i] if mm_vals is not None else 0.5
    score *= _pam_score(pam)
    return round(max(0.0, min(1.0, score)), 6)


def score_hits(hits: list[dict]) -> list[dict]:
    """Add cfd_score to each hit dict in-place and return the list."""
    for hit in hits:
        mm = hit["mismatches"]
        pam = hit.get("pam", "GG")
        if mm == 0:
            hit["cfd_score"] = round(_pam_score(pam), 6)
        else:
            hit["cfd_score"] = cfd_score(hit["guide_seq"], hit["target_seq"], pam)
    return hits
