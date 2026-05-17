"""
FR4 · Motion Analysis
Computes RMSD, RMSF, interdomain distance/angle, PCA, and clustering.

Requires: MDAnalysis, scikit-learn
Usage:
    python scripts/04_motion_analysis.py [--pdb 3LFM]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))
from scripts.utils import (
    DATA_DIR, RESULTS_DIR, DemoData, load_config, save_json, setup_logger,
)

log = setup_logger("04_motion")


def _check_mda():
    try:
        import MDAnalysis as mda
        return mda
    except ImportError:
        return None


def analyse_trajectory(pdb_id: str, replica: int, config: dict) -> dict:
    mda = _check_mda()
    if mda is None:
        raise ImportError("MDAnalysis not installed: pip install MDAnalysis")

    from sklearn.decomposition import PCA
    from sklearn.cluster import KMeans

    cfg = config["analysis"]
    traj_dir = DATA_DIR / "trajectories"
    prefix   = f"{pdb_id}_rep{replica}"
    top_path = DATA_DIR / "prepared" / f"{pdb_id}_prepared.pdb"
    dcd_path = traj_dir / f"{prefix}.dcd"

    if not top_path.exists() or not dcd_path.exists():
        raise FileNotFoundError(f"Missing topology or trajectory for {prefix}")

    log.info(f"  {prefix}: loading universe …")
    u = mda.Universe(str(top_path), str(dcd_path))
    protein = u.select_atoms("protein")

    # ── RMSD ─────────────────────────────────────────────────────────────────
    log.info(f"  {prefix}: computing RMSD …")
    from MDAnalysis.analysis import rms, align
    align.AlignTraj(u, u, select="backbone", in_memory=True).run()
    rmsd_analysis = rms.RMSD(protein, select=cfg["rmsd_selection"])
    rmsd_analysis.run()
    rmsd_df = pd.DataFrame({
        "time_ps": rmsd_analysis.times,
        "rmsd_A":  rmsd_analysis.rmsd[:, 2],
        "replica": replica,
    })

    # ── RMSF ─────────────────────────────────────────────────────────────────
    log.info(f"  {prefix}: computing RMSF …")
    rmsf_analysis = rms.RMSF(protein.select_atoms(cfg["rmsf_selection"]))
    rmsf_analysis.run()
    rmsf_df = pd.DataFrame({
        "residue": rmsf_analysis.atomgroup.resids,
        "rmsf_A":  rmsf_analysis.rmsf,
    })

    # ── Interdomain distance & angle ──────────────────────────────────────────
    log.info(f"  {prefix}: computing interdomain motion …")
    cat_res = config["fto"]["domains"]["catalytic"]
    cterm_res = config["fto"]["domains"]["cterminal"]
    cat_sel   = f"protein and backbone and resid {cat_res['start']}:{cat_res['end']}"
    cterm_sel = f"protein and backbone and resid {cterm_res['start']}:{cterm_res['end']}"

    dists, angles = [], []
    for ts in u.trajectory:
        cat_cg   = u.select_atoms(cat_sel).center_of_geometry()
        cterm_cg = u.select_atoms(cterm_sel).center_of_geometry()
        diff     = cterm_cg - cat_cg
        dists.append(float(np.linalg.norm(diff)))
        # Angle relative to initial orientation (simple proxy)
        angles.append(float(np.degrees(np.arctan2(diff[1], diff[0]))))

    interdomain_df = pd.DataFrame({
        "time_ps":               u.trajectory.times,
        "interdomain_dist_A":    dists,
        "interdomain_angle_deg": angles,
    })

    # ── PCA ───────────────────────────────────────────────────────────────────
    log.info(f"  {prefix}: PCA …")
    ca_sel = protein.select_atoms("name CA")
    n_comp = cfg["n_pca_components"]
    coords = []
    for ts in u.trajectory:
        coords.append(ca_sel.positions.flatten())
    coords_arr = np.array(coords)
    pca = PCA(n_components=n_comp)
    projections = pca.fit_transform(coords_arr)

    # ── Clustering ────────────────────────────────────────────────────────────
    log.info(f"  {prefix}: clustering …")
    k = cfg["n_clusters"]
    km = KMeans(n_clusters=k, random_state=42, n_init=10)
    clusters = km.fit_predict(projections[:, :3])

    pca_data = {
        "pc1": projections[:, 0].tolist(),
        "pc2": projections[:, 1].tolist(),
        "pc3": projections[:, 2].tolist(),
        "explained_variance": (pca.explained_variance_ratio_ * 100).tolist(),
        "clusters": clusters.tolist(),
        "time_ps": u.trajectory.times.tolist(),
    }

    # ── Save results ──────────────────────────────────────────────────────────
    out = RESULTS_DIR
    out.mkdir(parents=True, exist_ok=True)
    rmsd_df.to_csv(out / f"{prefix}_rmsd.csv", index=False)
    rmsf_df.to_csv(out / f"{prefix}_rmsf.csv", index=False)
    interdomain_df.to_csv(out / f"{prefix}_interdomain.csv", index=False)
    save_json(pca_data, out / f"{prefix}_pca.json")

    log.info(f"  {prefix}: results saved to results/")
    return {"pdb_id": pdb_id, "replica": replica, "status": "ok"}


def generate_demo(config: dict) -> None:
    """Write synthetic analysis results so the dashboard works without MD."""
    log.info("Generating demo motion-analysis results …")
    dd  = DemoData()
    out = RESULTS_DIR
    out.mkdir(parents=True, exist_ok=True)
    primary = config["fto"]["primary_structure"]
    prefix  = f"{primary}_rep0"

    dd.rmsd().to_csv(out / f"{prefix}_rmsd.csv", index=False)
    dd.rmsf().to_csv(out / f"{prefix}_rmsf.csv", index=False)
    dd.interdomain_motion().to_csv(out / f"{prefix}_interdomain.csv", index=False)
    save_json(dd.pca(), out / f"{prefix}_pca.json")
    log.info(f"Demo results written → results/{prefix}_*.csv / *.json")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb",     default=None)
    parser.add_argument("--replica", type=int, default=0)
    parser.add_argument("--demo",    action="store_true",
                        help="Generate synthetic demo data instead of running analysis")
    args = parser.parse_args()

    config = load_config()
    if args.demo:
        generate_demo(config)
        return

    pdb_ids = [args.pdb] if args.pdb else config["fto"]["pdb_ids"]
    for pdb_id in pdb_ids:
        log.info(f"=== {pdb_id} ===")
        try:
            analyse_trajectory(pdb_id, args.replica, config)
        except Exception as exc:
            log.error(f"{pdb_id}: FAILED — {exc}")


if __name__ == "__main__":
    main()
