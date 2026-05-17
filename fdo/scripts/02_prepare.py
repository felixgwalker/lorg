"""
FR2 · Structure Preparation
Repairs, protonates, solvates, and parameterises FTO structures.

Requires: pdbfixer, openmm
Usage:
    python scripts/02_prepare.py [--pdb 3LFM]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from scripts.utils import DATA_DIR, load_config, save_json, setup_logger

log = setup_logger("02_prepare")

STRUCTURES_DIR = DATA_DIR / "structures"
PREPARED_DIR   = DATA_DIR / "prepared"


def prepare_structure(pdb_id: str, config: dict) -> dict:
    """Fix, solvate, and serialise one FTO structure."""
    PREPARED_DIR.mkdir(parents=True, exist_ok=True)

    pdb_path = STRUCTURES_DIR / f"{pdb_id}.pdb"
    if not pdb_path.exists():
        raise FileNotFoundError(f"{pdb_path} not found – run 01_acquire.py first")

    try:
        from pdbfixer import PDBFixer
        import openmm
        import openmm.app as app
        import openmm.unit as unit
    except ImportError as exc:
        log.error(f"OpenMM / PDBFixer not installed: {exc}")
        log.error("Install with:  conda install -c conda-forge openmm pdbfixer")
        raise

    sim_cfg = config["simulation"]
    log.info(f"{pdb_id}: loading PDB …")
    fixer = PDBFixer(filename=str(pdb_path))

    log.info(f"{pdb_id}: finding and adding missing residues …")
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.4)

    log.info(f"{pdb_id}: solvating …")
    padding      = sim_cfg["padding_nm"] * unit.nanometer
    ionic_str    = sim_cfg["ionic_strength_molar"] * unit.molar
    fixer.addSolvent(padding=padding, ionicStrength=ionic_str)

    prepared_pdb = PREPARED_DIR / f"{pdb_id}_prepared.pdb"
    with open(prepared_pdb, "w") as fh:
        app.PDBFile.writeFile(fixer.topology, fixer.positions, fh)
    log.info(f"{pdb_id}: wrote {prepared_pdb.name}")

    log.info(f"{pdb_id}: building OpenMM system …")
    ff = app.ForceField(sim_cfg["forcefield"], sim_cfg["water_model"])
    system = ff.createSystem(
        fixer.topology,
        nonbondedMethod  = app.PME,
        nonbondedCutoff  = sim_cfg["nonbonded_cutoff_nm"] * unit.nanometer,
        constraints      = app.HBonds,
        rigidWater       = True,
        hydrogenMass     = 1.5 * unit.amu,
    )
    import openmm.openmm as mm
    system_xml = PREPARED_DIR / f"{pdb_id}_system.xml"
    with open(system_xml, "w") as fh:
        fh.write(mm.XmlSerializer.serialize(system))
    log.info(f"{pdb_id}: serialised system → {system_xml.name}")

    n_atoms = fixer.topology.getNumAtoms()
    n_res   = fixer.topology.getNumResidues()
    return {
        "pdb_id":        pdb_id,
        "prepared_pdb":  str(prepared_pdb),
        "system_xml":    str(system_xml),
        "n_atoms":       n_atoms,
        "n_residues":    n_res,
        "status":        "ok",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", help="Single PDB ID to prepare (default: all)")
    args = parser.parse_args()

    config  = load_config()
    pdb_ids = [args.pdb] if args.pdb else config["fto"]["pdb_ids"]

    results = {}
    for pdb_id in pdb_ids:
        log.info(f"=== Preparing {pdb_id} ===")
        try:
            results[pdb_id] = prepare_structure(pdb_id, config)
        except Exception as exc:
            log.error(f"{pdb_id}: FAILED — {exc}")
            results[pdb_id] = {"pdb_id": pdb_id, "status": "failed", "error": str(exc)}

    save_json(results, PREPARED_DIR / "preparation_log.json")
    ok = sum(1 for v in results.values() if v["status"] == "ok")
    log.info(f"Prepared {ok}/{len(results)} structures")


if __name__ == "__main__":
    main()
