"""
FR6 · Pocket Ranking
Scores and ranks candidate pockets by persistence, volume, druggability,
conservation, and enclosure.

Usage:
    python scripts/06_pocket_ranking.py [--demo]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))
from scripts.utils import (
    DATA_DIR, RESULTS_DIR, DemoData, load_config, load_json,
    save_json, setup_logger,
)

log = setup_logger("06_ranking")

POCKETS_DIR = DATA_DIR / "pockets"


# Weights for the composite pocket score
WEIGHTS = {
    "persistence":      0.30,
    "volume_norm":      0.25,
    "druggability":     0.20,
    "conservation":     0.15,
    "enclosure":        0.10,
}


def _volume_score(vol: float, v_min: float = 100, v_max: float = 700) -> float:
    """Map pocket volume to [0, 1], capped at v_max."""
    if vol <= v_min:
        return 0.0
    return min(1.0, (vol - v_min) / (v_max - v_min))


def _druggability_score(pocket: dict) -> float:
    """
    Estimate druggability from hydrophobic fraction, charged fraction,
    and volume.  Crude proxy; replace with fpocket dScore or DoGSiteScorer.
    """
    hydro   = pocket.get("hydrophobic_fraction", 0.5)
    charged = pocket.get("charged_fraction", 0.2)
    vol     = pocket.get("mean_volume_A3", 300)
    vol_s   = _volume_score(vol)
    # Preferred: moderate hydrophobicity, some charge, decent volume
    score = 0.5 * hydro + 0.2 * (1 - abs(charged - 0.3)) + 0.3 * vol_s
    return round(float(np.clip(score, 0.0, 1.0)), 3)


def _conservation_score(pocket: dict, conservation: dict | None) -> float:
    """Average ConSurf conservation score (1–9) for pocket residues, normalised."""
    if conservation is None:
        return 0.5  # neutral if unavailable
    scores = conservation.get("scores", [])
    residues = pocket.get("residues", [])
    if not residues or not scores:
        return 0.5
    vals = []
    for r in residues:
        idx = r - 1
        if 0 <= idx < len(scores):
            vals.append(scores[idx])
    if not vals:
        return 0.5
    return round(float(np.mean(vals) / 9.0), 3)


def rank_pockets(pockets: list[dict], conservation: dict | None = None) -> pd.DataFrame:
    rows = []
    for p in pockets:
        persist  = p.get("persistence", 0.0)
        vol_n    = _volume_score(p.get("mean_volume_A3", 0))
        drug     = p.get("druggability_score") or _druggability_score(p)
        cons     = _conservation_score(p, conservation)
        encl     = p.get("enclosure", 0.5)

        score = (
            WEIGHTS["persistence"]  * persist  +
            WEIGHTS["volume_norm"]  * vol_n    +
            WEIGHTS["druggability"] * drug     +
            WEIGHTS["conservation"] * cons     +
            WEIGHTS["enclosure"]    * encl
        )

        rows.append({
            "id":                  p["id"],
            "name":                p.get("name", ""),
            "pocket_score":        round(float(score), 3),
            "persistence":         round(float(persist), 3),
            "mean_volume_A3":      p.get("mean_volume_A3", 0),
            "max_volume_A3":       p.get("max_volume_A3", 0),
            "druggability_score":  round(float(drug), 3),
            "conservation_score":  round(float(cons), 3),
            "enclosure":           round(float(encl), 3),
            "n_residues":          len(p.get("residues", [])),
            "first_appearance_ps": p.get("first_appearance_ps", 0),
            "description":         p.get("description", ""),
            "color":               p.get("color", "#888888"),
        })

    df = pd.DataFrame(rows).sort_values("pocket_score", ascending=False).reset_index(drop=True)
    df.insert(0, "rank", df.index + 1)
    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--demo", action="store_true")
    args   = parser.parse_args()
    config = load_config()
    primary = config["fto"]["primary_structure"]

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # ── Load pocket and conservation data ─────────────────────────────────────
    pocket_file = POCKETS_DIR / f"{primary}_pockets.json"
    cons_file   = DATA_DIR / "conservation" / "fto_conservation.json"

    if args.demo or not pocket_file.exists():
        log.info("Using demo pocket data …")
        dd      = DemoData()
        pockets = dd.pockets()
        cons    = dd.conservation_scores()
    else:
        pockets = load_json(pocket_file)
        cons    = load_json(cons_file) if cons_file.exists() else None

    # ── Rank ──────────────────────────────────────────────────────────────────
    log.info(f"Ranking {len(pockets)} pockets …")
    df = rank_pockets(pockets, cons)

    csv_path  = RESULTS_DIR / "pocket_ranking.csv"
    json_path = RESULTS_DIR / "pocket_scores.json"
    df.to_csv(csv_path, index=False)
    save_json(df.to_dict(orient="records"), json_path)

    log.info(f"\nTop pockets:\n{df[['rank','id','name','pocket_score','persistence','mean_volume_A3']].to_string(index=False)}")
    log.info(f"Saved → {csv_path.name}  {json_path.name}")


if __name__ == "__main__":
    main()
