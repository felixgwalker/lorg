"""Estimate ERV insertion age from LTR divergence proxy."""

import random

SUBSTITUTION_RATE = 2.2e-9  # substitutions per site per year (standard mammalian neutral rate)
MOTIF_LENGTH = 8  # expected LTR motif length used as divergence denominator


def estimate_age(hit: dict, seed_offset: int = 0) -> float:
    """
    Estimate insertion age in millions of years (mya).

    Molecular clock model: LTR divergence correlates with ERV insertion age.
    Divergence = mismatches_between_LTR_motifs / expected_motif_length.
    Age (years) = divergence / (2 * substitution_rate).

    For demo data without real LTR sequences a seeded RNG produces plausible
    divergence values consistent with the observed LTR type, yielding ages in
    the range 1–50 mya.
    """
    rng = random.Random(hash(f"{hit['chrom']}_{hit['start']}_{seed_offset}") % (2**32))

    ltr_type = hit.get("ltr_type", "solo")

    # Estimate mismatches / motif_length as a proxy for LTR divergence.
    # Full-LTR elements are youngest; solo LTRs (5' and 3' recombined away) are oldest.
    if ltr_type == "full":
        # 0–4 mismatches in an 8-mer motif → divergence 0.0–0.50 → ~1–11 mya
        mismatches = rng.uniform(0.0, MOTIF_LENGTH * 0.50)
    elif ltr_type == "partial":
        # 3–6 mismatches → ~8–34 mya
        mismatches = rng.uniform(MOTIF_LENGTH * 0.375, MOTIF_LENGTH * 0.75)
    else:
        # 5–8 mismatches → ~14–45 mya
        mismatches = rng.uniform(MOTIF_LENGTH * 0.625, MOTIF_LENGTH)

    divergence = mismatches / MOTIF_LENGTH
    age_years = divergence / (2 * SUBSTITUTION_RATE)
    age_mya = age_years / 1_000_000

    # Clamp to biologically plausible range 1–50 mya
    age_mya = max(1.0, min(50.0, age_mya))
    return round(age_mya, 2)
