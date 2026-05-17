"""
Configuration for the CRISPR Base Editor Window Visualiser.

Contains all base editor profiles, activity window definitions, published
positional efficiency data, and output constants.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List

import numpy as np


# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

DEFAULT_BYSTANDER_THRESHOLD: float = 0.10   # absolute editing freq; above → HIGH risk
DEFAULT_OUTPUT_PREFIX: str = "be_output"


# ---------------------------------------------------------------------------
# Base editor profile
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class BaseEditorProfile:
    """Describes a base editor variant: window, target base, and efficiency."""

    name: str
    editor_class: str               # 'ABE' or 'CBE'
    window_start: int               # 1-indexed, PAM-distal
    window_end: int                 # 1-indexed, PAM-distal
    target_base: str                # 'A' (ABE) or 'C' (CBE)
    product_base: str               # 'G' (ABE) or 'T' (CBE)
    efficiency_profile: tuple       # 20 relative efficiencies, index 0 = pos 1
    max_absolute_efficiency: float  # published peak absolute editing frequency (0–1)
    description: str


def _gaussian_profile(
    center: float,
    sigma: float,
    window_start: int,
    window_end: int,
) -> tuple:
    """
    Build a length-20 relative efficiency profile using a Gaussian curve.

    Values outside the declared window are zeroed.  Peak is normalised to 1.0.
    """
    raw = []
    for pos in range(1, 21):
        if window_start <= pos <= window_end:
            val = np.exp(-0.5 * ((pos - center) / sigma) ** 2)
            raw.append(float(val))
        else:
            raw.append(0.0)

    peak = max(raw) or 1.0
    return tuple(round(v / peak, 4) for v in raw)


# ---------------------------------------------------------------------------
# Built-in editor profiles
# ---------------------------------------------------------------------------
#
# Efficiency profiles are Gaussian approximations calibrated to published
# positional preference data.  Absolute editing frequencies (max_absolute_efficiency)
# are representative values from the primary discovery papers.
#
# Sources:
#   ABE7.10  — Gaudelli et al. 2017, Nature 551:464
#   ABE8e    — Richter et al. 2020, Nature Biotechnology 38:883
#   BE3      — Komor et al. 2016, Nature 533:420
#   BE4max   — Koblan et al. 2018, Nature Biotechnology 36:843
#   evoAPOBEC — Thuronyi et al. 2019, Nature Biotechnology 37:1070
#   AncBE4max — Koblan et al. 2018, Nature Biotechnology 36:843

EDITOR_PROFILES: Dict[str, BaseEditorProfile] = {
    "ABE7.10": BaseEditorProfile(
        name="ABE7.10",
        editor_class="ABE",
        window_start=4,
        window_end=7,
        target_base="A",
        product_base="G",
        efficiency_profile=_gaussian_profile(5.5, 0.90, 4, 7),
        max_absolute_efficiency=0.50,
        description="Adenine base editor v7.10 (Gaudelli et al. 2017), narrow window 4–7",
    ),
    "ABE8e": BaseEditorProfile(
        name="ABE8e",
        editor_class="ABE",
        window_start=4,
        window_end=8,
        target_base="A",
        product_base="G",
        efficiency_profile=_gaussian_profile(5.5, 1.20, 4, 8),
        max_absolute_efficiency=0.75,
        description="Evolved adenine base editor 8e (Richter et al. 2020), broad window 4–8",
    ),
    "BE3": BaseEditorProfile(
        name="BE3",
        editor_class="CBE",
        window_start=4,
        window_end=8,
        target_base="C",
        product_base="T",
        efficiency_profile=_gaussian_profile(5.5, 1.10, 4, 8),
        max_absolute_efficiency=0.55,
        description="Cytosine base editor 3, rAPOBEC1–nCas9–UGI (Komor et al. 2016)",
    ),
    "BE4max": BaseEditorProfile(
        name="BE4max",
        editor_class="CBE",
        window_start=4,
        window_end=8,
        target_base="C",
        product_base="T",
        efficiency_profile=_gaussian_profile(5.5, 1.15, 4, 8),
        max_absolute_efficiency=0.70,
        description="Optimised CBE with enhanced deaminase (Koblan et al. 2018)",
    ),
    "evoAPOBEC": BaseEditorProfile(
        name="evoAPOBEC",
        editor_class="CBE",
        window_start=4,
        window_end=7,
        target_base="C",
        product_base="T",
        efficiency_profile=_gaussian_profile(5.0, 0.85, 4, 7),
        max_absolute_efficiency=0.48,
        description="Evolved APOBEC1-based CBE, narrow window 4–7 (Thuronyi et al. 2019)",
    ),
    "AncBE4max": BaseEditorProfile(
        name="AncBE4max",
        editor_class="CBE",
        window_start=4,
        window_end=8,
        target_base="C",
        product_base="T",
        efficiency_profile=_gaussian_profile(4.8, 0.90, 4, 8),
        max_absolute_efficiency=0.68,
        description="Ancestral APOBEC-based BE4max, tight PAM-distal bias (Koblan et al. 2018)",
    ),
}


def build_custom_profile(
    editor_class: str,
    window_start: int,
    window_end: int,
    max_absolute_efficiency: float = 0.60,
) -> BaseEditorProfile:
    """
    Construct a custom editor profile from user-specified window parameters.

    Args:
        editor_class:            'ABE' or 'CBE'.
        window_start:            PAM-distal window boundary (1-indexed).
        window_end:              PAM-proximal window boundary (1-indexed).
        max_absolute_efficiency: Ceiling absolute editing frequency (default 0.60).

    Returns:
        A BaseEditorProfile with a Gaussian efficiency centred on the window.
    """
    cls = editor_class.upper()
    if cls not in ("ABE", "CBE"):
        raise ValueError(f"editor_class must be 'ABE' or 'CBE'; got '{editor_class}'")
    if not (1 <= window_start <= window_end <= 20):
        raise ValueError(
            f"Window positions must satisfy 1 ≤ window_start ≤ window_end ≤ 20; "
            f"got {window_start}–{window_end}"
        )

    center = (window_start + window_end) / 2
    sigma = max((window_end - window_start) / 2.5, 0.5)

    return BaseEditorProfile(
        name=f"Custom-{cls}({window_start}-{window_end})",
        editor_class=cls,
        window_start=window_start,
        window_end=window_end,
        target_base="A" if cls == "ABE" else "C",
        product_base="G" if cls == "ABE" else "T",
        efficiency_profile=_gaussian_profile(center, sigma, window_start, window_end),
        max_absolute_efficiency=max_absolute_efficiency,
        description=f"Custom {cls} editor, window {window_start}–{window_end}",
    )


# ---------------------------------------------------------------------------
# Output filenames
# ---------------------------------------------------------------------------

OUTFILE_DUPLEX_PNG: str    = "{prefix}_duplex.png"
OUTFILE_DUPLEX_SVG: str    = "{prefix}_duplex.svg"
OUTFILE_EDITABILITY: str   = "{prefix}_editability.csv"
OUTFILE_OUTCOMES: str      = "{prefix}_outcomes.csv"
OUTFILE_BYSTANDERS: str    = "{prefix}_bystander_warnings.txt"
