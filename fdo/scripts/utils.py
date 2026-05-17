"""
Shared utilities for the FTO Pocket Discovery pipeline.
"""
from __future__ import annotations

import json
import logging
import sys
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd
import yaml

# ── Project paths ─────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).parent.parent
DATA_DIR     = PROJECT_ROOT / "data"
RESULTS_DIR  = PROJECT_ROOT / "results"
CONFIG_PATH  = PROJECT_ROOT / "config.yaml"


# ── Config ────────────────────────────────────────────────────────────────────

def load_config(path: Optional[Path] = None) -> dict:
    cfg_path = path or CONFIG_PATH
    with open(cfg_path) as f:
        return yaml.safe_load(f)


# ── Logging ───────────────────────────────────────────────────────────────────

def setup_logger(name: str, level: int = logging.INFO) -> logging.Logger:
    log = logging.getLogger(name)
    log.setLevel(level)
    if not log.handlers:
        handler = logging.StreamHandler(sys.stdout)
        fmt = logging.Formatter(
            "%(asctime)s | %(name)-22s | %(levelname)-8s | %(message)s",
            datefmt="%H:%M:%S",
        )
        handler.setFormatter(fmt)
        log.addHandler(handler)
    return log


# ── I/O helpers ───────────────────────────────────────────────────────────────

def save_json(data: Any, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(data, f, indent=2, default=str)


def load_json(path: Path) -> Any:
    with open(path) as f:
        return json.load(f)


def ensure_dirs() -> None:
    for sub in ["structures", "prepared", "trajectories", "pockets", "conservation"]:
        (DATA_DIR / sub).mkdir(parents=True, exist_ok=True)
    for sub in ["figures", "animations", "reports"]:
        (RESULTS_DIR / sub).mkdir(parents=True, exist_ok=True)


def results_path(name: str, subdir: str = "") -> Path:
    out = RESULTS_DIR / subdir if subdir else RESULTS_DIR
    out.mkdir(parents=True, exist_ok=True)
    return out / name


# ── Demo data generator ───────────────────────────────────────────────────────

class DemoData:
    """
    Generates scientifically plausible synthetic data for dashboard demo mode.
    Used whenever real pipeline results are absent.
    """

    N_FRAMES   = 2000
    N_RESIDUES = 498
    N_REPLICAS = 3

    def __init__(self, seed: int = 42):
        self.rng = np.random.default_rng(seed)
        self.time_ps = np.linspace(0, 20_000, self.N_FRAMES)

    # ── Trajectory metrics ────────────────────────────────────────────────────

    def rmsd(self, replica: int = 0) -> pd.DataFrame:
        t = self.time_ps
        core  = 1.5 + 1.2 * (1 - np.exp(-t / 1_000))
        noise = self.rng.standard_normal(len(t)) * 0.25
        drift = np.cumsum(self.rng.standard_normal(len(t)) * 0.008)
        vals  = np.clip(core + noise + drift, 0.3, 5.5)
        return pd.DataFrame({"time_ps": t, "rmsd_A": vals, "replica": replica})

    def rmsf(self) -> pd.DataFrame:
        r  = np.arange(1, self.N_RESIDUES + 1)
        b  = 1.2 + self.rng.exponential(0.4, self.N_RESIDUES)
        b[30:326]  *= 0.65   # catalytic core — rigid
        b[326:]    *= 1.25   # CTD — slightly flexible
        b[249:270] *= 1.9    # interdomain linker
        b[229:235] *= 0.35   # active site — very rigid
        b[:30]     *= 1.5    # N-terminus
        b[-20:]    *= 1.4    # C-terminus
        return pd.DataFrame({"residue": r, "rmsf_A": np.abs(b)})

    def interdomain_motion(self) -> pd.DataFrame:
        t = self.time_ps
        n = len(t)
        dist = 25.0 + 2.0 * np.sin(t / 5_000) + 0.4 * self.rng.standard_normal(n)
        dist = np.convolve(dist, np.ones(40) / 40, mode="same")

        angle = 4.5 + 2.5 * np.sin(t / 7_000 + 0.5)
        # Large opening event centred at frame 800
        opening = np.zeros(n)
        s, e = 750, 1_150
        opening[s:e] = 12.0 * np.sin(np.pi * np.arange(e - s) / (e - s))
        # Secondary event at frame 1 200
        s2, e2 = 1_180, 1_380
        opening[s2:e2] = 7.5 * np.sin(np.pi * np.arange(e2 - s2) / (e2 - s2))
        angle = np.convolve(angle + opening + 0.4 * self.rng.standard_normal(n),
                            np.ones(25) / 25, mode="same")
        return pd.DataFrame({
            "time_ps": t,
            "interdomain_dist_A":  dist,
            "interdomain_angle_deg": angle,
        })

    def pca(self) -> dict:
        n = self.N_FRAMES
        t = self.time_ps
        pc1 = np.cumsum(self.rng.standard_normal(n) * 0.08)
        pc1 -= pc1.mean()
        pc2 = 1.8 * np.sin(t / 3_500) + self.rng.standard_normal(n) * 0.4
        pc3 = self.rng.standard_normal(n) * 0.6

        clusters = np.zeros(n, dtype=int)
        clusters[(pc1 > 1.0) & (pc2 > 0.5)]  = 1
        clusters[pc1 < -1.2]                   = 2
        clusters[(pc1 > 0.3) & (pc2 < -1.0)]  = 3
        clusters[(clusters == 0) & (pc1 > 0.5)] = 4

        return {
            "pc1": pc1.tolist(),
            "pc2": pc2.tolist(),
            "pc3": pc3.tolist(),
            "explained_variance": [35.2, 18.7, 9.4, 6.1, 4.2, 3.1, 2.4, 1.9, 1.5, 1.2],
            "clusters": clusters.tolist(),
            "time_ps": t.tolist(),
        }

    # ── Pockets ───────────────────────────────────────────────────────────────

    def pockets(self) -> list[dict]:
        n, t = self.N_FRAMES, self.time_ps
        idx  = np.arange(n)

        def _envelope(centers_widths):
            env = np.zeros(n)
            for c, w in centers_widths:
                env += np.exp(-0.5 * ((idx - c) / (w / 3)) ** 2)
            return env

        # P1 – transient interdomain pocket
        p1_raw = _envelope([(800, 200), (1_190, 160), (1_600, 240), (1_820, 120)])
        p1_vol = np.clip(p1_raw * 460 + self.rng.standard_normal(n) * 28, 0, 720)

        # P2 – persistent active site pocket
        p2_vol = np.clip(
            320 + 45 * np.sin(t / 3_000) + self.rng.standard_normal(n) * 18, 200, 470
        )

        # P3 – transient CTD groove
        p3_raw = _envelope([(400, 110), (1_010, 85), (1_430, 115)])
        p3_vol = np.clip(p3_raw * 210 + self.rng.standard_normal(n) * 14, 0, 360)

        def _persist(vol, thr):
            return round(float(np.sum(vol > thr) / n), 3)

        def _first(vol, thr):
            mask = vol > thr
            return float(t[np.argmax(mask)]) if mask.any() else float(t[-1])

        return [
            {
                "id": "P1",
                "name": "Interdomain Interface Pocket",
                "rank": 1,
                "persistence": _persist(p1_vol, 100),
                "mean_volume_A3": round(float(p1_vol[p1_vol > 100].mean()), 1),
                "max_volume_A3":  round(float(p1_vol.max()), 1),
                "first_appearance_ps": _first(p1_vol, 100),
                "volume_trajectory": [round(float(v), 1) for v in p1_vol],
                "residues": [105, 106, 109, 130, 131, 203, 205, 240, 242, 330, 332, 335, 360, 362],
                "center": [12.4, 8.7, 22.1],
                "druggability_score": 0.74,
                "hydrophobic_fraction": 0.64,
                "charged_fraction": 0.21,
                "pocket_score": 0.82,
                "description": (
                    "Transient pocket at the catalytic–CTD interface, "
                    "exposed by interdomain rotation. Adjacent to H231/D233 "
                    "iron-binding triad."
                ),
                "color": "#FF6B35",
            },
            {
                "id": "P2",
                "name": "Catalytic Active Site",
                "rank": 2,
                "persistence": _persist(p2_vol, 150),
                "mean_volume_A3": round(float(p2_vol.mean()), 1),
                "max_volume_A3":  round(float(p2_vol.max()), 1),
                "first_appearance_ps": 0.0,
                "volume_trajectory": [round(float(v), 1) for v in p2_vol],
                "residues": [53, 131, 190, 205, 231, 233, 307],
                "center": [8.2, 12.1, 18.5],
                "druggability_score": 0.91,
                "hydrophobic_fraction": 0.43,
                "charged_fraction": 0.43,
                "pocket_score": 0.88,
                "description": (
                    "Canonical Fe²⁺ / α-KG binding pocket. Highly persistent "
                    "across all frames. Known drug target."
                ),
                "color": "#4A90D9",
            },
            {
                "id": "P3",
                "name": "CTD Surface Groove",
                "rank": 3,
                "persistence": _persist(p3_vol, 80),
                "mean_volume_A3": round(float(p3_vol[p3_vol > 80].mean()), 1) if np.any(p3_vol > 80) else 0.0,
                "max_volume_A3":  round(float(p3_vol.max()), 1),
                "first_appearance_ps": _first(p3_vol, 80),
                "volume_trajectory": [round(float(v), 1) for v in p3_vol],
                "residues": [340, 343, 346, 378, 381, 392, 394, 420, 421],
                "center": [18.3, 22.6, 5.4],
                "druggability_score": 0.48,
                "hydrophobic_fraction": 0.56,
                "charged_fraction": 0.22,
                "pocket_score": 0.51,
                "description": "Shallow surface groove on the C-terminal domain; transiently accessible.",
                "color": "#50C878",
            },
        ]

    # ── Conservation ──────────────────────────────────────────────────────────

    def conservation_scores(self) -> dict:
        scores = self.rng.integers(3, 10, self.N_RESIDUES).tolist()
        for r in [53, 131, 190, 205, 231, 233, 307]:
            scores[r - 1] = 9
        for i in range(30):
            scores[i] = max(1, scores[i] - 3)
        for i in range(self.N_RESIDUES - 20, self.N_RESIDUES):
            scores[i] = max(1, scores[i] - 2)
        return {
            "residues": list(range(1, self.N_RESIDUES + 1)),
            "scores": scores,
            "method": "ConSurf (demo)",
            "colormap": {
                "1": "#00CCFF", "2": "#33DDFF", "3": "#66EEFF",
                "4": "#99FFEE", "5": "#CCFFCC", "6": "#FFFF99",
                "7": "#FFCC66", "8": "#FF9933", "9": "#FF6600",
            },
        }

    # ── Events ────────────────────────────────────────────────────────────────

    def events(self) -> list[dict]:
        return [
            {
                "id": "E1", "type": "domain_motion",
                "time_ps": 7_500, "duration_ps": 800,
                "description": "Major interdomain rotation: catalytic domain tilts 9.3° relative to CTD",
                "frame": 750, "metrics": {"angle_change_deg": 9.3, "dist_change_A": 3.2},
            },
            {
                "id": "E2", "type": "contact_break",
                "time_ps": 7_600, "duration_ps": 1_800,
                "description": "Salt bridge E244–R365 breaks, correlating with P1 opening",
                "frame": 760, "metrics": {"dist_before_A": 3.2, "dist_after_A": 9.7},
            },
            {
                "id": "E3", "type": "pocket_open",
                "time_ps": 7_850, "duration_ps": 2_100,
                "description": "P1 interdomain pocket opens (volume > 300 Å³) after 9.3° domain rotation",
                "frame": 785, "pocket_id": "P1",
                "metrics": {"max_volume_A3": 487.2, "trigger_angle_deg": 9.3},
            },
            {
                "id": "E4", "type": "cluster_transition",
                "time_ps": 7_200, "duration_ps": 200,
                "description": "Trajectory shifts from cluster C1 (closed) to cluster C2 (open interdomain)",
                "frame": 720, "metrics": {"from_cluster": 0, "to_cluster": 1},
            },
            {
                "id": "E5", "type": "contact_form",
                "time_ps": 9_800, "duration_ps": 3_200,
                "description": "New H-bond Q108–N338 forms across domain interface",
                "frame": 980, "metrics": {"distance_A": 2.9},
            },
            {
                "id": "E6", "type": "pocket_open",
                "time_ps": 11_900, "duration_ps": 1_500,
                "description": "P1 second opening event: volume exceeds 400 Å³",
                "frame": 1_190, "pocket_id": "P1",
                "metrics": {"max_volume_A3": 412.5, "trigger_angle_deg": 7.8},
            },
            {
                "id": "E7", "type": "pocket_open",
                "time_ps": 4_000, "duration_ps": 900,
                "description": "P3 CTD groove transiently opens",
                "frame": 400, "pocket_id": "P3",
                "metrics": {"max_volume_A3": 198.3},
            },
        ]

    # ── Residue network ───────────────────────────────────────────────────────

    def residue_network(self) -> dict:
        return {
            "stable_contacts": [
                (50, 130), (51, 131), (100, 200), (150, 250),
                (200, 300), (231, 233), (231, 307), (233, 307),
            ],
            "breaking_contacts": [
                {"res1": 244, "res2": 365, "delta": -1.0, "label": "E244–R365 salt bridge"},
                {"res1": 108, "res2": 338, "delta": -0.8, "label": "Q108–N338 H-bond"},
                {"res1": 205, "res2": 332, "delta": -0.7, "label": "Y205–T332 H-bond"},
            ],
            "forming_contacts": [
                {"res1": 108, "res2": 340, "delta": 0.9, "label": "Q108–D340 new H-bond"},
                {"res1": 130, "res2": 335, "delta": 0.7, "label": "F130–L335 new contact"},
            ],
            "key_residues": [108, 130, 131, 205, 231, 233, 244, 307, 330, 332, 335, 338, 340, 360, 365],
        }

    # ── Structure metadata ────────────────────────────────────────────────────

    def structure_metadata(self) -> list[dict]:
        return [
            {
                "pdb_id": "3LFM", "resolution": 2.5, "method": "X-ray",
                "ligands": ["3DT", "FE2", "OGA"], "chains": ["A"], "n_residues": 420,
                "title": "Crystal structure of the fat mass and obesity associated (FTO) "
                         "protein reveals basis for its substrate specificity",
            },
            {
                "pdb_id": "4IE4", "resolution": 2.5, "method": "X-ray",
                "ligands": ["8XQ", "GOL", "ZN"], "chains": ["A"], "n_residues": 438,
                "title": "Crystal structure of the human FTO in complex with "
                         "5-carboxy-8-hydroxyquinoline (IOX1)",
            },
            {
                "pdb_id": "4ZS3", "resolution": 2.45, "method": "X-ray",
                "ligands": ["A4F", "AKG", "MN"], "chains": ["A"], "n_residues": 432,
                "title": "Structural complex of 5-aminofluorescein bound to the FTO protein",
            },
        ]

    # ── Pharmacophore ─────────────────────────────────────────────────────────

    def pharmacophore(self) -> dict:
        return {
            "pocket_id": "P1",
            "pocket_name": "Interdomain Interface Pocket",
            "volume_A3": 450,
            "surface_A2": 620,
            "enclosure": 0.68,
            "summary": (
                "Pharmacophore model for the P1 interdomain interface pocket, "
                "derived from the most-open trajectory frame (volume 450 Å³). "
                "Five binding features are predicted: 2 hydrogen-bond interactions "
                "(His231, Asp233), 1 aromatic centre (Tyr106), 1 hydrophobic patch "
                "(Leu109), and 1 positive-ionisable region (Arg96). A ligand "
                "complementary to these features could selectively occupy P1 and "
                "allosterically modulate FTO demethylase activity."
            ),
            "features": [
                {"type": "HBA", "position": [11.2, 9.1, 23.4], "radius": 1.5,
                 "importance": 0.90, "residue": "H231",
                 "description": "H-bond acceptor – His231 imidazole"},
                {"type": "HBD", "position": [13.8, 7.2, 21.0], "radius": 1.5,
                 "importance": 0.80, "residue": "D233",
                 "description": "H-bond donor – Asp233 carboxylate"},
                {"type": "AR",  "position": [10.5, 11.3, 24.7], "radius": 2.0,
                 "importance": 0.70, "residue": "Y106",
                 "description": "Aromatic ring – Tyr106 π-stacking"},
                {"type": "HYD", "position": [15.1, 8.8, 19.2], "radius": 2.5,
                 "importance": 0.75, "residue": "L109",
                 "description": "Hydrophobic – Leu109"},
                {"type": "PI",  "position": [9.8, 6.5, 22.8],  "radius": 1.8,
                 "importance": 0.65, "residue": "R96",
                 "description": "Positive ionisable – Arg96"},
            ],
        }
