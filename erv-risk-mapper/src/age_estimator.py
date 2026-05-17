"""Estimate ERV insertion age from LTR divergence proxy."""

import random

SUBSTITUTION_RATE = 2e-9  # substitutions per site per year


def estimate_age(hit: dict, seed_offset: int = 0) -> float:
    """
    Estimate insertion age in years.
    For full-LTR elements: use LTR divergence proxy (random in demo).
    For partial/solo: assume older (higher divergence).
    Returns age in millions of years (Mya).
    """
    rng = random.Random(hash(f"{hit['chrom']}_{hit['start']}_{seed_offset}") % (2**32))

    ltr_type = hit.get("ltr_type", "solo")
    if ltr_type == "full":
        divergence = rng.uniform(0.001, 0.05)
    elif ltr_type == "partial":
        divergence = rng.uniform(0.05, 0.15)
    else:
        divergence = rng.uniform(0.15, 0.40)

    age_years = divergence / (2 * SUBSTITUTION_RATE)
    return round(age_years / 1_000_000, 2)
