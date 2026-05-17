"""
FR3 · Molecular Dynamics Simulation
Runs OpenMM NPT replicate simulations for prepared FTO systems.

Requires: openmm
Usage:
    python scripts/03_simulate.py [--pdb 3LFM] [--replicas 3]
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from scripts.utils import DATA_DIR, RESULTS_DIR, load_config, save_json, setup_logger

log = setup_logger("03_simulate")

PREPARED_DIR    = DATA_DIR / "prepared"
TRAJECTORIES_DIR = DATA_DIR / "trajectories"


def run_simulation(pdb_id: str, replica: int, config: dict) -> dict:
    """Run one NPT replica and save trajectory."""
    TRAJECTORIES_DIR.mkdir(parents=True, exist_ok=True)

    try:
        import openmm
        import openmm.app as app
        import openmm.unit as unit
        import openmm.openmm as mm
    except ImportError as exc:
        raise ImportError(
            "OpenMM not found. Install with: conda install -c conda-forge openmm"
        ) from exc

    sim_cfg = config["simulation"]
    prefix  = f"{pdb_id}_rep{replica}"

    # ── Load prepared system ──────────────────────────────────────────────────
    prepared_pdb = PREPARED_DIR / f"{pdb_id}_prepared.pdb"
    system_xml   = PREPARED_DIR / f"{pdb_id}_system.xml"
    if not prepared_pdb.exists() or not system_xml.exists():
        raise FileNotFoundError(f"Prepared files missing for {pdb_id} – run 02_prepare.py first")

    pdb    = app.PDBFile(str(prepared_pdb))
    with open(system_xml) as f:
        system = mm.XmlSerializer.deserialize(f.read())

    # ── Integrator & barostat ─────────────────────────────────────────────────
    dt          = sim_cfg["timestep_fs"] * unit.femtosecond
    temperature = sim_cfg["temperature_K"] * unit.kelvin
    pressure    = sim_cfg["pressure_atm"] * unit.atmosphere

    integrator = mm.LangevinMiddleIntegrator(temperature, 1.0 / unit.picosecond, dt)
    integrator.setRandomNumberSeed(42 + replica * 1_000)
    system.addForce(mm.MonteCarloBarostat(pressure, temperature, sim_cfg["barostat_frequency"]))

    # ── Platform ──────────────────────────────────────────────────────────────
    platform = None
    for name in ("CUDA", "OpenCL", "CPU"):
        try:
            platform = mm.Platform.getPlatformByName(name)
            log.info(f"  {prefix}: using {name} platform")
            break
        except Exception:
            continue

    simulation = app.Simulation(pdb.topology, system, integrator, platform)
    simulation.context.setPositions(pdb.positions)

    # ── Energy minimisation ───────────────────────────────────────────────────
    log.info(f"  {prefix}: minimising energy ({sim_cfg['n_minimization_steps']} steps) …")
    simulation.minimizeEnergy(maxIterations=sim_cfg["n_minimization_steps"])

    # ── NVT equilibration ─────────────────────────────────────────────────────
    equil_steps = int(sim_cfg["equilibration_ns"] * 1e6 / sim_cfg["timestep_fs"])
    log.info(f"  {prefix}: equilibrating {sim_cfg['equilibration_ns']} ns ({equil_steps:,} steps) …")
    simulation.step(equil_steps)

    # ── Production NPT ────────────────────────────────────────────────────────
    prod_steps    = int(sim_cfg["production_ns"] * 1e6 / sim_cfg["timestep_fs"])
    save_interval = int(sim_cfg["save_interval_ps"] * 1_000 / sim_cfg["timestep_fs"])

    dcd_path = TRAJECTORIES_DIR / f"{prefix}.dcd"
    log.path = TRAJECTORIES_DIR / f"{prefix}.log"

    simulation.reporters.append(app.DCDReporter(str(dcd_path), save_interval))
    simulation.reporters.append(
        app.StateDataReporter(
            str(TRAJECTORIES_DIR / f"{prefix}.log"),
            save_interval,
            step=True, time=True,
            potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
            temperature=True, volume=True, density=True,
            progress=True, remainingTime=True, speed=True,
            totalSteps=prod_steps,
        )
    )

    log.info(f"  {prefix}: production run {sim_cfg['production_ns']} ns ({prod_steps:,} steps) …")
    t0 = time.perf_counter()
    simulation.step(prod_steps)
    elapsed = time.perf_counter() - t0

    log.info(f"  {prefix}: done in {elapsed/60:.1f} min → {dcd_path.name}")
    return {
        "pdb_id":    pdb_id,
        "replica":   replica,
        "dcd":       str(dcd_path),
        "log":       str(TRAJECTORIES_DIR / f"{prefix}.log"),
        "prod_ns":   sim_cfg["production_ns"],
        "elapsed_s": round(elapsed, 1),
        "status":    "ok",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb",      default=None, help="Single PDB ID")
    parser.add_argument("--replicas", type=int, default=None, help="Number of replicas")
    args = parser.parse_args()

    config   = load_config()
    pdb_ids  = [args.pdb] if args.pdb else config["fto"]["pdb_ids"]
    n_rep    = args.replicas or config["simulation"]["n_replicas"]

    results = {}
    for pdb_id in pdb_ids:
        for rep in range(n_rep):
            key = f"{pdb_id}_rep{rep}"
            log.info(f"=== {key} ===")
            try:
                results[key] = run_simulation(pdb_id, rep, config)
            except Exception as exc:
                log.error(f"{key}: FAILED — {exc}")
                results[key] = {"key": key, "status": "failed", "error": str(exc)}

    save_json(results, TRAJECTORIES_DIR / "simulation_log.json")
    ok = sum(1 for v in results.values() if v["status"] == "ok")
    log.info(f"Completed {ok}/{len(results)} simulations")


if __name__ == "__main__":
    main()
