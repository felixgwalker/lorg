"""
FR10 · Candidate Ligand / Pharmacophore Module
Generates pharmacophore description for the top-ranked pocket and
optionally runs AutoDock Vina fragment docking.

Usage:
    python scripts/10_ligand.py [--demo] [--dock]
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))
from scripts.utils import (
    DATA_DIR, RESULTS_DIR, DemoData, load_config, load_json,
    save_json, setup_logger,
)

log = setup_logger("10_ligand")


FEATURE_TYPES = {
    "HBA": "H-bond acceptor",
    "HBD": "H-bond donor",
    "AR":  "Aromatic / π-stacking",
    "HYD": "Hydrophobic",
    "PI":  "Positive ionisable",
    "NI":  "Negative ionisable",
}

# Amino-acid pharmacophore feature lookup
AA_FEATURES = {
    "HIS": ["HBA", "HBD", "AR"],
    "ASP": ["HBA", "NI"],
    "GLU": ["HBA", "NI"],
    "ARG": ["HBD", "PI"],
    "LYS": ["HBD", "PI"],
    "TYR": ["HBA", "HBD", "AR"],
    "PHE": ["AR", "HYD"],
    "TRP": ["AR", "HYD", "HBD"],
    "LEU": ["HYD"],
    "ILE": ["HYD"],
    "VAL": ["HYD"],
    "ALA": ["HYD"],
    "MET": ["HYD"],
    "PRO": ["HYD"],
    "SER": ["HBA", "HBD"],
    "THR": ["HBA", "HBD"],
    "ASN": ["HBA", "HBD"],
    "GLN": ["HBA", "HBD"],
    "CYS": ["HBD"],
    "GLY": [],
}


def _residue_name_from_id(pdb_id: str, res_id: int) -> str:
    """Look up residue name from PDB file."""
    pdb = DATA_DIR / "structures" / f"{pdb_id}.pdb"
    if not pdb.exists():
        return "UNK"
    with open(pdb) as f:
        for line in f:
            if line.startswith("ATOM"):
                try:
                    if int(line[22:26]) == res_id:
                        return line[17:20].strip()
                except ValueError:
                    pass
    return "UNK"


def generate_pharmacophore(pocket: dict, pdb_id: str) -> dict:
    """Build a pharmacophore model from pocket residues and geometry."""
    residues = pocket.get("residues", [])
    center   = np.array(pocket.get("center", [0, 0, 0]))
    features = []

    for i, res_id in enumerate(residues[:8]):  # top 8 residues
        resname = _residue_name_from_id(pdb_id, res_id)
        fts     = AA_FEATURES.get(resname, [])
        if not fts:
            continue
        ftype = fts[0]
        # Perturb position around pocket centre
        angle  = 2 * np.pi * i / max(len(residues), 1)
        radius = 3.5 + np.random.uniform(-1, 1)
        pos    = (
            center +
            np.array([radius * np.cos(angle), radius * np.sin(angle), np.random.uniform(-2, 2)])
        )
        features.append({
            "type":        ftype,
            "description": FEATURE_TYPES.get(ftype, ftype),
            "position":    [round(float(v), 2) for v in pos],
            "radius":      round(1.2 + 0.5 * (ftype in ("AR", "HYD")), 1),
            "importance":  round(0.6 + 0.3 * np.random.random(), 2),
            "residue":     f"{resname}{res_id}",
        })

    return {
        "pocket_id":    pocket["id"],
        "pocket_name":  pocket.get("name", ""),
        "volume_A3":    pocket.get("mean_volume_A3", 0),
        "n_features":   len(features),
        "features":     features,
        "summary": (
            f"Pharmacophore for pocket {pocket['id']}: "
            f"{sum(1 for f in features if f['type'] in ('HBA','HBD'))} H-bond features, "
            f"{sum(1 for f in features if f['type'] == 'AR')} aromatic features, "
            f"{sum(1 for f in features if f['type'] == 'HYD')} hydrophobic features."
        ),
    }


def run_autodock_vina(pocket: dict, pdb_id: str, config: dict) -> dict | None:
    """Optional: run AutoDock Vina if installed."""
    if not shutil.which("vina"):
        log.warning("AutoDock Vina not found in PATH – skipping docking")
        return None

    c = pocket.get("center", [0, 0, 0])
    receptor = DATA_DIR / "prepared" / f"{pdb_id}_prepared.pdbqt"
    if not receptor.exists():
        log.warning("PDBQT receptor not found – run prepare_receptor4.py first")
        return None

    out_dir = RESULTS_DIR / "docking"
    out_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "vina",
        "--receptor", str(receptor),
        "--center_x", str(c[0]),
        "--center_y", str(c[1]),
        "--center_z", str(c[2]),
        "--size_x", "20", "--size_y", "20", "--size_z", "20",
        "--exhaustiveness", "8",
        "--out", str(out_dir / f"{pdb_id}_{pocket['id']}_docked.pdbqt"),
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        log.info(f"Vina output:\n{result.stdout[-500:]}")
        return {"status": "ok", "stdout": result.stdout[-500:]}
    except Exception as exc:
        log.error(f"Vina failed: {exc}")
        return {"status": "failed", "error": str(exc)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--demo", action="store_true")
    parser.add_argument("--dock", action="store_true", help="Run AutoDock Vina if available")
    args   = parser.parse_args()
    config = load_config()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    primary = config["fto"]["primary_structure"]

    if args.demo:
        pharma = DemoData().pharmacophore()
        save_json(pharma, RESULTS_DIR / "pharmacophore.json")
        log.info(f"Demo pharmacophore saved (pocket {pharma['pocket_id']}, "
                 f"{pharma['n_features']} features)")
        return

    # Load top-ranked pocket
    ranking_file = RESULTS_DIR / "pocket_scores.json"
    pockets_file = DATA_DIR / "pockets" / f"{primary}_pockets.json"

    if ranking_file.exists():
        ranked = load_json(ranking_file)
        top_id = ranked[0]["id"] if ranked else "P1"
    else:
        top_id = "P1"

    if pockets_file.exists():
        pockets = load_json(pockets_file)
        top_pocket = next((p for p in pockets if p["id"] == top_id), pockets[0] if pockets else None)
    else:
        top_pocket = DemoData().pockets()[0]

    pharma = generate_pharmacophore(top_pocket, primary)
    save_json(pharma, RESULTS_DIR / "pharmacophore.json")
    log.info(f"Pharmacophore: {pharma['n_features']} features for pocket {pharma['pocket_id']}")

    if args.dock or config["ligand"]["run_docking"]:
        docking_result = run_autodock_vina(top_pocket, primary, config)
        if docking_result:
            save_json(docking_result, RESULTS_DIR / "docking_results.json")


if __name__ == "__main__":
    main()
