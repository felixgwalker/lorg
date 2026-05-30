"""Report generation for proxy-species edit-burden calculator.

Outputs:
  - edit_burden_summary.csv  — counts, weights, weighted contribution per class
                                plus estimated CRISPR cycles (10 edits / cycle)
  - prioritized_edits.csv    — variants sorted by descending impact score
  - edit_burden.json         — machine-readable burden dict
"""

from __future__ import annotations

import json
import os
import math

import pandas as pd


# Assumed maximum edits deliverable in one CRISPR engineering cycle.
EDITS_PER_CYCLE = 10


def _crispr_cycles(n_edits: int, edits_per_cycle: int = EDITS_PER_CYCLE) -> int:
    """Return minimum number of CRISPR cycles required for *n_edits* edits."""
    if n_edits <= 0:
        return 0
    return math.ceil(n_edits / edits_per_cycle)


def write_burden_summary(burden: dict, output_dir: str) -> str:
    """Write edit burden summary CSV with per-class counts and CRISPR cycle estimates.

    Parameters
    ----------
    burden : dict
        Output of compute_burden() — must contain 'class_counts', 'total_edits',
        'weighted_burden', 'normalized_burden_per_mb'.
    output_dir : str
        Directory in which to write ``edit_burden_summary.csv``.

    Returns
    -------
    str
        Path to the written CSV file.
    """
    path = os.path.join(output_dir, "edit_burden_summary.csv")
    weights = {
        "SNV": 1,
        "SMALL_INS": 3,
        "SMALL_DEL": 3,
        "LARGE_INS": 10,
        "LARGE_DEL": 10,
        "SV_INS": 50,
        "SV_DEL": 50,
    }

    rows = []
    for vc, count in burden["class_counts"].items():
        w = weights.get(vc, 1)
        rows.append({
            "variant_class": vc,
            "count": count,
            "weight": w,
            "weighted_contribution": count * w,
            "estimated_crispr_cycles": _crispr_cycles(count),
        })

    total_edits = burden["total_edits"]
    rows.append({
        "variant_class": "TOTAL",
        "count": total_edits,
        "weight": "",
        "weighted_contribution": burden["weighted_burden"],
        "estimated_crispr_cycles": _crispr_cycles(total_edits),
    })

    df = pd.DataFrame(rows)
    df.to_csv(path, index=False)
    return path


def write_prioritized_edits(variants, output_dir: str, top_n=None) -> str:
    """Write variants sorted by descending impact score to CSV.

    Parameters
    ----------
    variants : list
        List of variant dicts or Variant dataclass instances.
    output_dir : str
        Output directory.
    top_n : int or None
        If set, only write the top *top_n* variants.

    Returns
    -------
    str
        Path to the written CSV file.
    """
    path = os.path.join(output_dir, "prioritized_edits.csv")

    def _score(v):
        return v.impact_score if hasattr(v, "impact_score") else v.get("impact_score", 0)

    sorted_vars = sorted(variants, key=lambda v: -_score(v))
    if top_n:
        sorted_vars = sorted_vars[:top_n]

    rows = []
    for v in sorted_vars:
        _is_dc = hasattr(v, "__dataclass_fields__")
        if _is_dc:
            rows.append({
                "chrom": v.chrom,
                "pos": v.position,
                "ref": v.ref_allele,
                "alt": v.alt_allele,
                "type": v.type,
                "variant_class": v.variant_class,
                "impact_category": v.impact_category,
                "impact_score": v.impact_score,
                "weight": v.weight,
            })
        else:
            rows.append({
                "chrom": v.get("chrom", ""),
                "pos": v.get("pos", 0),
                "ref": v.get("ref", ""),
                "alt": v.get("alt", ""),
                "type": v.get("type", ""),
                "variant_class": v.get("variant_class", ""),
                "impact_category": v.get("impact_category", ""),
                "impact_score": v.get("impact_score", 0),
                "weight": v.get("weight", 1),
            })

    df = (
        pd.DataFrame(rows)
        if rows
        else pd.DataFrame(
            columns=[
                "chrom", "pos", "ref", "alt", "type", "variant_class",
                "impact_category", "impact_score", "weight",
            ]
        )
    )
    df.to_csv(path, index=False)
    return path


def write_burden_json(burden: dict, output_dir: str) -> str:
    """Write the burden dict as JSON, enriched with CRISPR cycle estimates.

    Parameters
    ----------
    burden : dict
        Output of compute_burden().
    output_dir : str
        Output directory.

    Returns
    -------
    str
        Path to the written JSON file.
    """
    path = os.path.join(output_dir, "edit_burden.json")

    # Augment the burden dict with CRISPR cycle data without mutating caller's copy.
    enriched = dict(burden)
    total = burden.get("total_edits", 0)
    enriched["estimated_crispr_cycles"] = _crispr_cycles(total)
    enriched["edits_per_cycle_assumption"] = EDITS_PER_CYCLE

    # Per-class cycle estimates.
    enriched["crispr_cycles_by_class"] = {
        vc: _crispr_cycles(count)
        for vc, count in burden.get("class_counts", {}).items()
    }

    with open(path, "w") as fh:
        json.dump(enriched, fh, indent=2)
    return path
