"""
FR5 · Pocket Detection
Detects transient cavities using fpocket (preferred) or a grid-based fallback.

Usage:
    python scripts/05_pocket_detection.py [--pdb 3LFM] [--method grid|fpocket]
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))
from scripts.utils import (
    DATA_DIR, RESULTS_DIR, DemoData, load_config, save_json, setup_logger,
)

log = setup_logger("05_pockets")

POCKETS_DIR = DATA_DIR / "pockets"


# ── Grid-based fallback ───────────────────────────────────────────────────────

def _read_atom_coords(pdb_text: str) -> np.ndarray:
    coords = []
    for line in pdb_text.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])
            except ValueError:
                continue
    return np.array(coords) if coords else np.zeros((0, 3))


def _grid_pockets(coords: np.ndarray, spacing: float = 1.5,
                  probe_r: float = 1.4, burial_cutoff: int = 4) -> list[dict]:
    """
    Minimal grid-based pocket detection.
    Returns a list of pocket dicts sorted by volume descending.
    """
    if len(coords) == 0:
        return []

    lo = coords.min(axis=0) - 5.0
    hi = coords.max(axis=0) + 5.0
    gx = np.arange(lo[0], hi[0], spacing)
    gy = np.arange(lo[1], hi[1], spacing)
    gz = np.arange(lo[2], hi[2], spacing)

    # For speed, randomly sample grid points rather than dense enumeration
    rng = np.random.default_rng(0)
    n_sample = min(20_000, int(len(gx) * len(gy) * len(gz)))
    pts = rng.uniform(lo, hi, (n_sample, 3))

    # Burial: count nearby protein atoms
    from scipy.spatial import cKDTree  # type: ignore
    tree = cKDTree(coords)
    counts = tree.query_ball_point(pts, r=8.0, return_length=True)  # type: ignore
    buried_mask = counts >= burial_cutoff

    # Remove points inside VdW radius
    dists, _ = tree.query(pts)
    buried_mask &= dists > probe_r

    buried_pts = pts[buried_mask]
    if len(buried_pts) == 0:
        return []

    # Cluster buried points into pockets
    from sklearn.cluster import DBSCAN
    labels = DBSCAN(eps=3.0, min_samples=5).fit_predict(buried_pts)

    pockets = []
    for lbl in sorted(set(labels)):
        if lbl == -1:
            continue
        cluster_pts = buried_pts[labels == lbl]
        centre = cluster_pts.mean(axis=0)
        vol    = len(cluster_pts) * spacing ** 3
        pockets.append({
            "center": centre.tolist(),
            "volume_A3": round(float(vol), 1),
            "n_grid_points": int(len(cluster_pts)),
        })

    pockets.sort(key=lambda p: p["volume_A3"], reverse=True)
    return pockets[:10]


def detect_pockets_grid(pdb_id: str, config: dict) -> list[dict]:
    """Run grid detection on every Nth frame of a trajectory."""
    traj_dir = DATA_DIR / "trajectories"
    dcd_path = traj_dir / f"{pdb_id}_rep0.dcd"
    top_path = DATA_DIR / "prepared" / f"{pdb_id}_prepared.pdb"

    try:
        import MDAnalysis as mda
    except ImportError:
        raise ImportError("MDAnalysis required for trajectory pocket detection")

    cfg    = config["pocket_detection"]
    stride = 20  # sample every 20th frame for speed

    log.info(f"  {pdb_id}: loading trajectory …")
    u        = mda.Universe(str(top_path), str(dcd_path))
    protein  = u.select_atoms("protein")

    frame_pockets = []
    for i, ts in enumerate(u.trajectory[::stride]):
        pdb_txt = ""
        for atom in protein.atoms:
            pdb_txt += (
                f"ATOM  {atom.id:5d} {atom.name:<4s} {atom.resname:<3s} "
                f"{atom.segid:1s}{atom.resid:4d}    "
                f"{atom.position[0]:8.3f}{atom.position[1]:8.3f}{atom.position[2]:8.3f}"
                f"  1.00  0.00\n"
            )
        pockets = _grid_pockets(
            protein.positions,
            spacing=cfg["grid_spacing_A"],
            probe_r=cfg["probe_radius_A"],
            burial_cutoff=cfg["burial_cutoff"],
        )
        frame_pockets.append({
            "frame_idx": i * stride,
            "time_ps":   float(ts.time),
            "pockets":   pockets[: cfg["max_pockets_per_frame"]],
        })
        if i % 20 == 0:
            log.info(f"  {pdb_id}: frame {i * stride} → {len(pockets)} pockets found")

    return frame_pockets


# ── fpocket integration ───────────────────────────────────────────────────────

def _fpocket_available() -> bool:
    return shutil.which("fpocket") is not None


def detect_pockets_fpocket(pdb_id: str, config: dict) -> list[dict]:
    """Run fpocket on the primary structure."""
    pdb_path = DATA_DIR / "structures" / f"{pdb_id}.pdb"
    if not pdb_path.exists():
        raise FileNotFoundError(f"{pdb_path} not found")

    POCKETS_DIR.mkdir(parents=True, exist_ok=True)
    out_dir = POCKETS_DIR / f"{pdb_id}_fpocket"

    cmd = ["fpocket", "-f", str(pdb_path), "-o", str(out_dir)]
    log.info(f"  {pdb_id}: running fpocket …")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    if result.returncode != 0:
        raise RuntimeError(f"fpocket failed: {result.stderr[:500]}")

    # Parse fpocket info file
    info_file = out_dir / f"{pdb_id}_info.txt"
    pockets = _parse_fpocket_info(info_file) if info_file.exists() else []
    log.info(f"  {pdb_id}: fpocket found {len(pockets)} pockets")
    return pockets


def _parse_fpocket_info(info_path: Path) -> list[dict]:
    pockets, current = [], {}
    with open(info_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("Pocket"):
                if current:
                    pockets.append(current)
                current = {"id": line}
            elif ":" in line:
                k, _, v = line.partition(":")
                try:
                    current[k.strip()] = float(v.strip())
                except ValueError:
                    current[k.strip()] = v.strip()
    if current:
        pockets.append(current)
    return pockets


# ── Demo ──────────────────────────────────────────────────────────────────────

def generate_demo(config: dict) -> None:
    log.info("Generating demo pocket-detection results …")
    POCKETS_DIR.mkdir(parents=True, exist_ok=True)
    dd = DemoData()
    primary = config["fto"]["primary_structure"]
    pockets = dd.pockets()

    # Build fake per-frame data from volume trajectories
    n_frames  = DemoData.N_FRAMES
    time_ns   = np.linspace(0, 20, n_frames)
    stride    = 20
    frame_data = []

    for fi in range(0, n_frames, stride):
        frame_pockets = []
        for p in pockets:
            vol = p["volume_trajectory"][fi]
            if vol > 50:
                frame_pockets.append({
                    "id":        p["id"],
                    "volume_A3": vol,
                    "center":    p["center"],
                })
        frame_data.append({
            "frame_idx": fi,
            "time_ps":   float(time_ns[fi] * 1_000),
            "pockets":   frame_pockets,
        })

    save_json(
        {"pdb_id": primary, "n_frames": n_frames, "frames": frame_data},
        POCKETS_DIR / f"{primary}_pocket_trajectory.json",
    )
    save_json(pockets, POCKETS_DIR / f"{primary}_pockets.json")
    log.info("Demo pocket data written → data/pockets/")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb",    default=None)
    parser.add_argument("--method", choices=["grid", "fpocket", "auto"], default="auto")
    parser.add_argument("--demo",   action="store_true")
    args = parser.parse_args()

    config = load_config()
    if args.demo:
        generate_demo(config)
        return

    method  = args.method
    pdb_ids = [args.pdb] if args.pdb else config["fto"]["pdb_ids"]

    for pdb_id in pdb_ids:
        log.info(f"=== {pdb_id} ===")
        try:
            if method in ("auto", "fpocket") and _fpocket_available():
                pockets = detect_pockets_fpocket(pdb_id, config)
                save_json(pockets, POCKETS_DIR / f"{pdb_id}_pockets.json")
            elif method in ("auto", "grid"):
                log.info(f"  fpocket not found — using grid method")
                frame_pockets = detect_pockets_grid(pdb_id, config)
                save_json(
                    {"pdb_id": pdb_id, "frames": frame_pockets},
                    POCKETS_DIR / f"{pdb_id}_pocket_trajectory.json",
                )
            else:
                log.error(f"Method {method} unavailable")
        except Exception as exc:
            log.error(f"{pdb_id}: FAILED — {exc}")


if __name__ == "__main__":
    main()
