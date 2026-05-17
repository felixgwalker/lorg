"""
Bayesian per-read classifier for ancient DNA authenticity.

For each read, computes P(ancient | observed terminal mismatches) and
P(contaminated | observed terminal mismatches) using Bayes' theorem with
binomial likelihoods evaluated in log-space to prevent numerical underflow.

The model assumes:
  - Authentic ancient reads exhibit C→T and G→A terminal mismatches at a rate
    equal to the library's fitted deamination rate (from decay_model).
  - Contaminated/modern reads exhibit such mismatches only at the background
    sequencing error rate (BACKGROUND_ERROR_RATE from config).
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass

from scipy.stats import binom

from src.config import AUTH_THRESHOLD, BACKGROUND_ERROR_RATE, CONT_THRESHOLD
from src.damage_profiler import ReadFeatures

logger = logging.getLogger(__name__)

# Minimum probability to avoid log(0)
_LOG_FLOOR = -1e10


@dataclass
class ReadClassification:
    """Classification result for a single read."""

    read_id: str
    classification: str             # "authentic" | "contaminated" | "ambiguous"
    posterior_ancient: float        # P(ancient | data)
    posterior_contaminated: float   # P(contaminated | data)
    ct_terminal: int                # raw C→T count at 5' terminus
    ga_terminal: int                # raw G→A count at 3' terminus
    read_length: int
    uninformative: bool             # True if no C or G reference bases in terminal windows


@dataclass
class ClassificationSummary:
    """Library-level summary of per-read classifications."""

    n_authentic: int
    n_contaminated: int
    n_ambiguous: int
    n_total: int
    fraction_authentic: float
    fraction_contaminated: float
    fraction_ambiguous: float
    mean_posterior_ancient: float


def classify_reads(
    read_features: list[ReadFeatures],
    library_deamination_rate: float,
    prior_ancient: float = 0.9,
    auth_threshold: float = AUTH_THRESHOLD,
    cont_threshold: float = CONT_THRESHOLD,
) -> tuple[list[ReadClassification], ClassificationSummary]:
    """
    Classify each read as authentic, contaminated, or ambiguous.

    Uses log-space Bayesian inference with binomial likelihoods to avoid
    numerical underflow on reads with many or very few terminal mismatches.

    Edge cases handled:
    - No informative terminal bases (n_ct + n_ga == 0): forced "ambiguous".
    - library_deamination_rate <= BACKGROUND_ERROR_RATE: posterior dominated
      by prior; a warning is emitted once.

    Args:
        read_features:            List from damage_profiler.profile_damage().
        library_deamination_rate: Fitted deamination rate (ancient model rate).
        prior_ancient:            Library-level prior P(ancient).
        auth_threshold:           Posterior >= this -> "authentic".
        cont_threshold:           Posterior <= this -> "contaminated".

    Returns:
        Tuple of (list[ReadClassification], ClassificationSummary).
    """
    p_anc = library_deamination_rate
    p_mod = BACKGROUND_ERROR_RATE

    if p_anc <= p_mod * 1.5:
        logger.warning(
            "Library deamination rate (%.5f) is close to background error rate "
            "(%.5f). Classifications will be dominated by the prior (%.2f).",
            p_anc, p_mod, prior_ancient,
        )

    classifications: list[ReadClassification] = []
    for rf in read_features:
        rc = _classify_single(rf, p_anc, p_mod, prior_ancient, auth_threshold, cont_threshold)
        classifications.append(rc)

    # ── summary statistics ──────────────────────────────────────────────────
    n_total = len(classifications)
    n_authentic    = sum(1 for c in classifications if c.classification == "authentic")
    n_contaminated = sum(1 for c in classifications if c.classification == "contaminated")
    n_ambiguous    = sum(1 for c in classifications if c.classification == "ambiguous")

    if n_total > 0:
        frac_auth = n_authentic / n_total
        frac_cont = n_contaminated / n_total
        frac_amb  = n_ambiguous / n_total
        mean_post = sum(c.posterior_ancient for c in classifications) / n_total
    else:
        frac_auth = frac_cont = frac_amb = mean_post = 0.0

    summary = ClassificationSummary(
        n_authentic=n_authentic,
        n_contaminated=n_contaminated,
        n_ambiguous=n_ambiguous,
        n_total=n_total,
        fraction_authentic=frac_auth,
        fraction_contaminated=frac_cont,
        fraction_ambiguous=frac_amb,
        mean_posterior_ancient=mean_post,
    )

    return classifications, summary


def _classify_single(
    rf: ReadFeatures,
    p_anc: float,
    p_mod: float,
    prior_anc: float,
    auth_threshold: float,
    cont_threshold: float,
) -> ReadClassification:
    """
    Classify a single read using log-space Bayesian inference.

    Log-space algorithm:
        log_L_anc = logpmf(k_ct, n_ct, p_anc) + logpmf(k_ga, n_ga, p_anc)
        log_L_mod = logpmf(k_ct, n_ct, p_mod) + logpmf(k_ga, n_ga, p_mod)
        log_post_anc_unnorm = log_L_anc + log(prior_anc)
        log_post_mod_unnorm = log_L_mod + log(1 - prior_anc)
        log_Z = log-sum-exp([log_post_anc_unnorm, log_post_mod_unnorm])
        posterior_ancient = exp(log_post_anc_unnorm - log_Z)

    If n_ct + n_ga == 0 (no informative bases), forces "ambiguous" without
    performing the Bayes update.
    """
    k_ct = rf.ct_terminal_count
    n_ct = rf.ct_terminal_opportunities
    k_ga = rf.ga_terminal_count
    n_ga = rf.ga_terminal_opportunities

    uninformative = (n_ct + n_ga) == 0

    if uninformative:
        # No reference C or G bases in the terminal window; no information.
        posterior_ancient = prior_anc
        posterior_contaminated = 1.0 - prior_anc
        classification = "ambiguous"
    else:
        # Log-likelihoods: sum of two binomial log-PMFs
        log_L_anc = (
            _safe_logpmf(k_ct, n_ct, p_anc) +
            _safe_logpmf(k_ga, n_ga, p_anc)
        )
        log_L_mod = (
            _safe_logpmf(k_ct, n_ct, p_mod) +
            _safe_logpmf(k_ga, n_ga, p_mod)
        )

        log_prior_anc = math.log(prior_anc)
        log_prior_mod = math.log(1.0 - prior_anc)

        log_post_anc_unnorm = log_L_anc + log_prior_anc
        log_post_mod_unnorm = log_L_mod + log_prior_mod

        # Numerically stable log-sum-exp
        log_Z = _log_sum_exp(log_post_anc_unnorm, log_post_mod_unnorm)

        posterior_ancient       = math.exp(log_post_anc_unnorm - log_Z)
        posterior_contaminated  = math.exp(log_post_mod_unnorm - log_Z)

        # Clamp floating-point rounding artefacts
        posterior_ancient      = max(0.0, min(1.0, posterior_ancient))
        posterior_contaminated = max(0.0, min(1.0, posterior_contaminated))

        if posterior_ancient >= auth_threshold:
            classification = "authentic"
        elif posterior_ancient <= cont_threshold:
            classification = "contaminated"
        else:
            classification = "ambiguous"

    return ReadClassification(
        read_id=rf.read_id,
        classification=classification,
        posterior_ancient=posterior_ancient,
        posterior_contaminated=posterior_contaminated,
        ct_terminal=k_ct,
        ga_terminal=k_ga,
        read_length=rf.read_length,
        uninformative=uninformative,
    )


def _safe_logpmf(k: int, n: int, p: float) -> float:
    """
    Compute binom.logpmf(k, n, p) with a floor to avoid -inf propagation.

    When n == 0, scipy returns 0.0 (probability 1 of drawing 0 from 0 trials),
    which is correct and requires no special handling.  We only floor the result
    to avoid -inf when the event is technically impossible but n > 0.

    Args:
        k: Number of observed successes.
        n: Number of trials.
        p: Success probability.

    Returns:
        Log probability, floored at _LOG_FLOOR.
    """
    if n == 0:
        return 0.0
    val = float(binom.logpmf(k, n, p))
    return max(val, _LOG_FLOOR)


def _log_sum_exp(a: float, b: float) -> float:
    """
    Compute log(exp(a) + exp(b)) in a numerically stable way.

    Args:
        a: First log-value.
        b: Second log-value.

    Returns:
        log(exp(a) + exp(b)).
    """
    if a >= b:
        return a + math.log1p(math.exp(b - a))
    return b + math.log1p(math.exp(a - b))
