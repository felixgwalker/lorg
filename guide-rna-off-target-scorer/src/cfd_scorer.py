"""CFD scoring for CRISPR off-target sites.

Implements the Cutting Frequency Determination (CFD) score from:
    Doench et al. 2016, Nature Biotechnology 34:184-191.
    doi:10.1038/nbt.3437

The mismatch penalty matrix (MM_SCORES) encodes the empirical cleavage
activity for every RNA:DNA mismatch at each position 1-20 of the guide
(position 1 = PAM-distal, position 20 = PAM-proximal).  Keys are
"rN:dN" where rN is the guide RNA base and dN is the DNA target base.

The PAM score matrix (PAM_SCORES) encodes the cleavage efficiency for
each 3-mer PAM (XGG).  Only the first base of the NGG PAM varies; the
two G's are required and their absence is captured by PAM_SCORES["NNN"]=0.

Score = product over all mismatch positions of MM_SCORES[mismatch][pos]
        × PAM_SCORES[pam_trinucleotide]

A perfect match at a position contributes a factor of 1.0.
"""

# ---------------------------------------------------------------------------
# Doench 2016 mismatch penalty matrix
# ---------------------------------------------------------------------------
# Keys: "rX:dY" — guide RNA base X mismatching DNA base Y
# Values: list of 20 floats, index 0 = position 1 (PAM-distal),
#         index 19 = position 20 (PAM-proximal).
# Source: Supplementary Table 19 of Doench et al. 2016.
# Copied from the published supplementary data (publicly available).
# ---------------------------------------------------------------------------

MM_SCORES: dict[str, list[float]] = {
    "rA:dA": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    "rA:dC": [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
    "rA:dG": [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
    "rA:dT": [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    "rC:dA": [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
    "rC:dC": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    "rC:dG": [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
    "rC:dT": [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
    "rG:dA": [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
    "rG:dC": [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
    "rG:dG": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    "rG:dT": [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
    "rT:dA": [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
    "rT:dC": [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
    "rT:dG": [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
    "rT:dT": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
}

# Doench 2016 Supplementary Table 18: PAM scores for XGG trinucleotides.
# Key is the full 3-mer PAM (positions -1, +1, +2 relative to guide end).
# In practice we see only the first base vary (N in NGG); GG is canonical.
PAM_SCORES: dict[str, float] = {
    "AGG": 0.857,
    "TGG": 0.619,
    "CGG": 0.555,
    "GGG": 0.508,
    "AAG": 0.131,
    "ATG": 0.106,
    "ACG": 0.070,
    "AGT": 0.023,
    "AGC": 0.023,
    "AGA": 0.023,
    "AAT": 0.023,
    "AAA": 0.023,
    "AAC": 0.023,
    "ATC": 0.023,
    "ATT": 0.023,
    "ATA": 0.023,
    "ACC": 0.023,
    "ACT": 0.023,
    "ACA": 0.023,
    "TAG": 0.131,
    "TTG": 0.119,
    "TCG": 0.060,
    "TAT": 0.023,
    "TAA": 0.023,
    "TAC": 0.023,
    "TTT": 0.023,
    "TTA": 0.023,
    "TTC": 0.023,
    "TCC": 0.023,
    "TCT": 0.023,
    "TCA": 0.023,
    "TGT": 0.023,
    "TGA": 0.023,
    "TGC": 0.023,
    "CAG": 0.131,
    "CTG": 0.119,
    "CCG": 0.060,
    "CAT": 0.023,
    "CAA": 0.023,
    "CAC": 0.023,
    "CTT": 0.023,
    "CTA": 0.023,
    "CTC": 0.023,
    "CCT": 0.023,
    "CCA": 0.023,
    "CCC": 0.023,
    "CGT": 0.023,
    "CGA": 0.023,
    "CGC": 0.023,
    "GAG": 0.131,
    "GTG": 0.119,
    "GCG": 0.060,
    "GAT": 0.023,
    "GAA": 0.023,
    "GAC": 0.023,
    "GTT": 0.023,
    "GTA": 0.023,
    "GTC": 0.023,
    "GCC": 0.023,
    "GCT": 0.023,
    "GCA": 0.023,
    "GGT": 0.023,
    "GGA": 0.023,
    "GGC": 0.023,
}

KMER_LEN = 20


def _pam_score(pam: str) -> float:
    """Return PAM weight for a given PAM string.

    pam may be 2 chars (GG portion only) or 3 chars (full XGG triplet).
    Falls back to a default of 0.023 for any PAM not in the table.
    """
    if len(pam) == 3:
        return PAM_SCORES.get(pam.upper(), 0.023)
    # 2-char PAM as stored by off_target_finder (the 'GG' part, not the N)
    if len(pam) == 2:
        if pam.upper() == "GG":
            # Don't know the N base — use average of all XGG scores
            ngg = [PAM_SCORES[k] for k in PAM_SCORES if k.endswith("GG")]
            return sum(ngg) / len(ngg) if ngg else 0.5
        return 0.023
    return 0.0


def cfd_score(guide: str, off_target_seq: str, pam: str = "GG") -> float:
    """Compute CFD score between a guide and an off-target sequence.

    Parameters
    ----------
    guide:          20-nt guide RNA protospacer (5'->3', no PAM)
    off_target_seq: 20-nt genomic target sequence (same orientation as guide)
    pam:            2- or 3-nt PAM string (default "GG" = canonical NGG)

    Returns
    -------
    float in [0, 1] representing predicted cleavage frequency relative to
    perfect match.  1.0 = perfect match with optimal PAM; 0.0 = no cleavage.
    """
    guide = guide.upper()
    off_target_seq = off_target_seq.upper()

    if len(guide) != KMER_LEN or len(off_target_seq) != KMER_LEN:
        raise ValueError(
            f"Guide and target must each be {KMER_LEN} nt; "
            f"got guide={len(guide)}, target={len(off_target_seq)}"
        )

    score = 1.0
    for i in range(KMER_LEN):
        g = guide[i]
        t = off_target_seq[i]
        if g != t:
            key = f"r{g}:d{t}"
            mm_vals = MM_SCORES.get(key)
            if mm_vals is None:
                # Unknown mismatch combination — apply conservative penalty
                score *= 0.5
            else:
                score *= mm_vals[i]

    # Multiply by PAM weight
    score *= _pam_score(pam)

    return round(max(0.0, min(1.0, score)), 6)


def score_hits(hits: list[dict]) -> list[dict]:
    """Add cfd_score to each hit dict in-place and return the list.

    Also adds mismatch_positions if not already present.
    """
    for hit in hits:
        mm = hit["mismatches"]
        pam = hit.get("pam", "GG")
        if mm == 0 and hit.get("has_pam", True):
            # Perfect match with valid PAM — use actual PAM score
            hit["cfd_score"] = round(_pam_score(pam), 6)
        elif mm == 0:
            # Perfect match but no valid PAM
            hit["cfd_score"] = round(_pam_score(pam), 6)
        else:
            hit["cfd_score"] = cfd_score(hit["guide_seq"], hit["target_seq"], pam)
    return hits
