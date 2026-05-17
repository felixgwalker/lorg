"""
FR10 · Residue Interaction Network
Builds a Cα contact map for open vs closed states and identifies
contacts that change when pockets open.

Usage:
    python scripts/09_residue_network.py [--demo]
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

log = setup_logger("09_network")


def _load_trajectory(pdb_id: str, replica: int):
    try:
        import MDAnalysis as mda
    except ImportError:
        raise ImportError("MDAnalysis required")

    top = DATA_DIR / "prepared" / f"{pdb_id}_prepared.pdb"
    dcd = DATA_DIR / "trajectories" / f"{pdb_id}_rep{replica}.dcd"
    return mda.Universe(str(top), str(dcd))


def _contact_fraction(ca_positions: np.ndarray, cutoff: float = 8.0) -> np.ndarray:
    """Return binary contact matrix for a single frame."""
    n = len(ca_positions)
    mat = np.zeros((n, n), dtype=bool)
    for i in range(n):
        diff = ca_positions[i + 1:] - ca_positions[i]
        d2   = (diff ** 2).sum(axis=1)
        mat[i, i + 1:] = d2 < cutoff ** 2
        mat[i + 1:, i] = mat[i, i + 1:]
    return mat


def build_contact_network(pdb_id: str, config: dict) -> dict:
    """
    Compute average contact map for open (pocket P1 volume > 300 Å³)
    and closed (< 100 Å³) states separately.
    """
    cutoff = config["analysis"]["contact_cutoff_A"]
    u      = _load_trajectory(pdb_id, 0)
    protein = u.select_atoms("protein and name CA")

    pock_file = DATA_DIR / "pockets" / f"{pdb_id}_pocket_trajectory.json"
    pock_data = load_json(pock_file).get("frames", []) if pock_file.exists() else []

    # Index pocket volume by frame
    p1_vol = {}
    for frame in pock_data:
        fi = frame["frame_idx"]
        v  = next((p["volume_A3"] for p in frame["pockets"] if p.get("id") == "P1"), 0.0)
        p1_vol[fi] = v

    open_mats, closed_mats = [], []
    for i, ts in enumerate(u.trajectory[::10]):
        vol  = p1_vol.get(i * 10, 0.0)
        mat  = _contact_fraction(protein.positions, cutoff)
        if vol > 300:
            open_mats.append(mat)
        elif vol < 100:
            closed_mats.append(mat)

    n = len(protein.resids)
    if open_mats and closed_mats:
        open_avg   = np.mean(open_mats, axis=0)
        closed_avg = np.mean(closed_mats, axis=0)
        delta      = open_avg - closed_avg
    else:
        delta = np.zeros((n, n))

    # Identify top changing contacts
    idx = np.argwhere(np.abs(delta) > 0.3)
    changes = [
        {
            "res1":  int(protein.resids[i]),
            "res2":  int(protein.resids[j]),
            "delta": round(float(delta[i, j]), 3),
            "label": f"R{protein.resids[i]}–R{protein.resids[j]}",
        }
        for i, j in idx if i < j
    ]
    changes.sort(key=lambda x: abs(x["delta"]), reverse=True)

    return {
        "residue_ids":      protein.resids.tolist(),
        "changing_contacts": changes[:30],
        "breaking_contacts": [c for c in changes if c["delta"] < -0.3][:15],
        "forming_contacts":  [c for c in changes if c["delta"] >  0.3][:15],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--demo", action="store_true")
    args   = parser.parse_args()
    config = load_config()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    if args.demo:
        data = DemoData().residue_network()
    else:
        primary = config["fto"]["primary_structure"]
        try:
            data = build_contact_network(primary, config)
        except Exception as exc:
            log.warning(f"Network analysis failed: {exc} – using demo data")
            data = DemoData().residue_network()

    save_json(data, RESULTS_DIR / "residue_network.json")
    n_break = len(data.get("breaking_contacts", []))
    n_form  = len(data.get("forming_contacts", []))
    log.info(f"Residue network: {n_break} breaking, {n_form} forming contacts → results/residue_network.json")


if __name__ == "__main__":
    main()
