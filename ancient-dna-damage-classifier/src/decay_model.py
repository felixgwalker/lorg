"""
Geometric decay model fitting for ancient DNA damage profiles.

Fits f(x) = amplitude * (1 - rate)^x + background to observed C→T and G→A
substitution rates using non-linear least squares (scipy.optimize.curve_fit).

The fitted amplitude represents the damage signal above background, the decay
rate captures how quickly damage frequency drops off from the terminus, and the
background captures the residual error rate at internal positions.
"""

from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass

import numpy as np
from scipy.optimize import OptimizeWarning, curve_fit

from src.config import DECAY_PARAM_BOUNDS, DECAY_PARAM_P0, GRADE_THRESHOLDS
from src.damage_profiler import DamageProfile

logger = logging.getLogger(__name__)


@dataclass
class DecayFit:
    """Result of fitting the geometric decay model to one terminus."""

    terminus: str           # "5prime" or "3prime"
    amplitude: float        # fitted damage magnitude above background
    rate: float             # per-position geometric decay rate
    background: float       # asymptotic substitution rate at internal positions
    r_squared: float        # coefficient of determination
    fitted_values: np.ndarray
    converged: bool         # False if curve_fit raised RuntimeError/OptimizeWarning
    signal_quality: str     # "strong" | "moderate" | "weak" | "absent"


@dataclass
class ModelResult:
    """Fitted decay models for both termini plus library-level deamination rate."""

    five_prime: DecayFit
    three_prime: DecayFit
    library_deamination_rate: float     # used as "ancient" rate in Bayesian classifier
    overall_signal_quality: str         # worst of the two terminus quality grades


def _geometric_decay(
    x: np.ndarray, amplitude: float, rate: float, background: float
) -> np.ndarray:
    """
    Geometric decay function for damage frequency.

    f(x) = amplitude * (1 - rate)^x + background

    Args:
        x:          Array of position indices (0 = most terminal).
        amplitude:  Damage signal above background.
        rate:       Per-position decay rate (0 < rate < 1).
        background: Asymptotic rate at internal positions.

    Returns:
        Predicted substitution rates at each position.
    """
    return amplitude * (1.0 - rate) ** x + background


def _compute_r_squared(observed: np.ndarray, fitted: np.ndarray) -> float:
    """
    Compute the coefficient of determination R².

    Returns 0.0 when variance in observed is zero (flat signal) or fitted has
    no variation relative to observed mean.

    Args:
        observed: Observed substitution rates.
        fitted:   Model-predicted rates at the same positions.

    Returns:
        R² clamped to [0.0, 1.0].
    """
    ss_tot = np.sum((observed - np.mean(observed)) ** 2)
    if ss_tot == 0.0:
        return 0.0
    ss_res = np.sum((observed - fitted) ** 2)
    r2 = 1.0 - ss_res / ss_tot
    return float(np.clip(r2, 0.0, 1.0))


def _assess_signal_quality(r_squared: float, amplitude: float) -> str:
    """
    Classify damage signal quality using GRADE_THRESHOLDS from config.

    Args:
        r_squared: Coefficient of determination from decay fit.
        amplitude: Fitted damage amplitude above background.

    Returns:
        One of: "strong", "moderate", "weak", "absent".
    """
    for min_r2, min_amp, label in GRADE_THRESHOLDS:
        if r_squared >= min_r2 and amplitude >= min_amp:
            return label
    return "absent"


def fit_decay(
    rates: np.ndarray,
    terminus: str,
) -> DecayFit:
    """
    Fit the geometric decay model to observed substitution rates.

    Positions with rate == 0.0 where the denominator was also 0 are excluded
    before fitting (they carry no information).  If curve_fit fails to converge,
    a flat fallback model is returned with converged=False.

    Args:
        rates:   Observed substitution rates, shape (n_terminal,).
        terminus: "5prime" or "3prime" (for labelling only).

    Returns:
        A DecayFit dataclass with fitted parameters and quality assessment.
    """
    positions = np.arange(len(rates), dtype=np.float64)

    # Use all positions for fitting; curve_fit handles near-zero values fine.
    # If all rates are zero, fall back immediately.
    if np.all(rates == 0.0):
        fitted_values = np.zeros_like(rates)
        return DecayFit(
            terminus=terminus,
            amplitude=0.0,
            rate=0.0,
            background=0.0,
            r_squared=0.0,
            fitted_values=fitted_values,
            converged=False,
            signal_quality="absent",
        )

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error", OptimizeWarning)
            popt, _ = curve_fit(
                _geometric_decay,
                positions,
                rates,
                p0=DECAY_PARAM_P0,
                bounds=DECAY_PARAM_BOUNDS,
                maxfev=10_000,
            )
        amplitude, rate, background = float(popt[0]), float(popt[1]), float(popt[2])
        fitted_values = _geometric_decay(positions, amplitude, rate, background)
        r_squared = _compute_r_squared(rates, fitted_values)
        converged = True

    except (RuntimeError, OptimizeWarning, ValueError) as exc:
        logger.warning(
            "Decay model did not converge for %s terminus (%s). "
            "Using flat fallback model.",
            terminus, exc,
        )
        background = float(np.mean(rates))
        amplitude, rate = 0.0, 0.0
        fitted_values = np.full_like(rates, background)
        r_squared = 0.0
        converged = False

    signal_quality = _assess_signal_quality(r_squared, amplitude)

    return DecayFit(
        terminus=terminus,
        amplitude=amplitude,
        rate=rate,
        background=background,
        r_squared=r_squared,
        fitted_values=fitted_values,
        converged=converged,
        signal_quality=signal_quality,
    )


def _quality_rank(quality: str) -> int:
    """Return a numeric rank for comparing signal quality (lower = worse)."""
    order = {"absent": 0, "weak": 1, "moderate": 2, "strong": 3}
    return order.get(quality, 0)


def fit_model(profile: DamageProfile) -> ModelResult:
    """
    Fit decay models to both 5' C→T and 3' G→A damage profiles.

    The library_deamination_rate is derived as the mean of the two fitted
    (amplitude + background) values, representing the expected per-base
    deamination probability at the most terminal position.  It is clamped to a
    minimum equal to BACKGROUND_ERROR_RATE so the Bayesian classifier always
    has a non-degenerate ancient model.

    Args:
        profile: The DamageProfile returned by damage_profiler.profile_damage().

    Returns:
        A ModelResult with fitted DecayFit objects and library-level stats.
    """
    from src.config import BACKGROUND_ERROR_RATE

    five_fit = fit_decay(profile.ct_rate, "5prime")
    three_fit = fit_decay(profile.ga_rate, "3prime")

    # Library deamination rate: mean peak damage at position 0 for both termini
    peak_ct = five_fit.amplitude + five_fit.background
    peak_ga = three_fit.amplitude + three_fit.background
    raw_rate = (peak_ct + peak_ga) / 2.0
    library_deamination_rate = max(raw_rate, BACKGROUND_ERROR_RATE)

    # Overall quality = worse of the two
    if _quality_rank(five_fit.signal_quality) <= _quality_rank(three_fit.signal_quality):
        overall_signal_quality = five_fit.signal_quality
    else:
        overall_signal_quality = three_fit.signal_quality

    if library_deamination_rate <= BACKGROUND_ERROR_RATE * 1.5:
        logger.warning(
            "Library deamination rate (%.4f) is at or near background error rate "
            "(%.4f). Damage signal may be absent; classification posteriors will be "
            "dominated by the prior.",
            library_deamination_rate,
            BACKGROUND_ERROR_RATE,
        )

    return ModelResult(
        five_prime=five_fit,
        three_prime=three_fit,
        library_deamination_rate=library_deamination_rate,
        overall_signal_quality=overall_signal_quality,
    )
