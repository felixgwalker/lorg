"""Published base editor profiles (window, target base, efficiency).

Salvaged from crispr-base-editor-window-visualiser (deleted stage1f).
Sources: Gaudelli 2017, Richter 2020, Komor 2016, Koblan 2018, Thuronyi 2019.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Dict
import numpy as np


DEFAULT_BYSTANDER_THRESHOLD: float = 0.10


@dataclass(frozen=True)
class BaseEditorProfile:
    name: str
    editor_class: str
    window_start: int
    window_end: int
    target_base: str
    product_base: str
    efficiency_profile: tuple
    max_absolute_efficiency: float
    description: str


def _gaussian_profile(center: float, sigma: float, window_start: int, window_end: int) -> tuple:
    raw = []
    for pos in range(1, 21):
        if window_start <= pos <= window_end:
            raw.append(float(np.exp(-0.5 * ((pos - center) / sigma) ** 2)))
        else:
            raw.append(0.0)
    peak = max(raw) or 1.0
    return tuple(round(v / peak, 4) for v in raw)


EDITOR_PROFILES: Dict[str, BaseEditorProfile] = {
    "ABE7.10": BaseEditorProfile("ABE7.10", "ABE", 4, 7, "A", "G",
        _gaussian_profile(5.5, 0.90, 4, 7), 0.50, "Gaudelli et al. 2017, window 4-7"),
    "ABE8e": BaseEditorProfile("ABE8e", "ABE", 4, 8, "A", "G",
        _gaussian_profile(5.5, 1.20, 4, 8), 0.75, "Richter et al. 2020, window 4-8"),
    "BE3": BaseEditorProfile("BE3", "CBE", 4, 8, "C", "T",
        _gaussian_profile(5.5, 1.10, 4, 8), 0.55, "Komor et al. 2016"),
    "BE4max": BaseEditorProfile("BE4max", "CBE", 4, 8, "C", "T",
        _gaussian_profile(5.5, 1.15, 4, 8), 0.70, "Koblan et al. 2018"),
    "evoAPOBEC": BaseEditorProfile("evoAPOBEC", "CBE", 4, 7, "C", "T",
        _gaussian_profile(5.0, 0.85, 4, 7), 0.48, "Thuronyi et al. 2019"),
    "AncBE4max": BaseEditorProfile("AncBE4max", "CBE", 4, 8, "C", "T",
        _gaussian_profile(4.8, 0.90, 4, 8), 0.68, "Koblan et al. 2018"),
    "CBE-NG": BaseEditorProfile("CBE-NG", "CBE", 4, 8, "C", "T",
        _gaussian_profile(5.5, 1.15, 4, 8), 0.65, "Nishimasu et al. 2018, NG PAM"),
    "ABE-NG": BaseEditorProfile("ABE-NG", "ABE", 4, 8, "A", "G",
        _gaussian_profile(5.5, 1.20, 4, 8), 0.70, "Nishimasu et al. 2018, NG PAM"),
}

EDITOR_PAM_REQUIREMENT: Dict[str, str] = {
    "ABE7.10": "NGG", "ABE8e": "NGG", "BE3": "NGG", "BE4max": "NGG",
    "evoAPOBEC": "NGG", "AncBE4max": "NGG", "CBE-NG": "NG", "ABE-NG": "NG",
}


def get_pam_requirement(editor_name: str) -> str:
    return EDITOR_PAM_REQUIREMENT.get(editor_name, "NGG")


def pam_matches(pam_seq: str, pattern: str) -> bool:
    seq, pat = pam_seq.upper(), pattern.upper()
    if len(seq) < len(pat):
        return False
    iupac = {"N": set("ACGT"), "R": set("AG"), "Y": set("CT"),
             "W": set("AT"), "S": set("CG"), "K": set("GT"), "M": set("AC")}
    for s, p in zip(seq, pat):
        allowed = iupac.get(p, {p})
        if s not in allowed:
            return False
    return True
