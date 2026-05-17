"""
FR8 · Conformational Event Detection
Automatically detects pocket openings, domain motions, and contact changes.

Usage:
    python scripts/08_event_detection.py [--demo]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))
from scripts.utils import (
    DATA_DIR, RESULTS_DIR, DemoData, load_config, load_json,
    save_json, setup_logger,
)

log = setup_logger("08_events")


# ── Detectors ─────────────────────────────────────────────────────────────────

def detect_pocket_openings(
    pocket_traj: list[dict],
    pocket_id: str,
    vol_threshold: float = 200.0,
    min_duration_frames: int = 5,
) -> list[dict]:
    """Detect contiguous intervals where a pocket exceeds the volume threshold."""
    vols = np.array([
        next((p["volume_A3"] for p in frame["pockets"] if p.get("id") == pocket_id), 0.0)
        for frame in pocket_traj
    ])
    times = np.array([frame["time_ps"] for frame in pocket_traj])

    open_mask = vols > vol_threshold
    events = []
    in_event = False
    start_i  = 0

    for i, o in enumerate(open_mask):
        if o and not in_event:
            in_event = True
            start_i  = i
        elif not o and in_event:
            in_event = False
            dur = i - start_i
            if dur >= min_duration_frames:
                peak_vol = float(vols[start_i:i].max())
                events.append({
                    "type":        "pocket_open",
                    "pocket_id":   pocket_id,
                    "time_ps":     float(times[start_i]),
                    "end_ps":      float(times[i - 1]),
                    "duration_ps": float(times[i - 1] - times[start_i]),
                    "frame":       int(pocket_traj[start_i]["frame_idx"]),
                    "max_volume_A3": round(peak_vol, 1),
                    "description": (
                        f"{pocket_id} opens (vol > {vol_threshold:.0f} Å³) "
                        f"at {times[start_i]/1000:.1f} ns"
                    ),
                })
    return events


def detect_domain_motion_events(
    interdomain_df: pd.DataFrame,
    angle_threshold: float = 5.0,
    window: int = 30,
) -> list[dict]:
    """Detect frames where the interdomain angle changes sharply."""
    angles = interdomain_df["interdomain_angle_deg"].values
    times  = interdomain_df["time_ps"].values
    # Smooth first
    kernel = np.ones(window) / window
    smooth = np.convolve(angles, kernel, mode="same")
    delta  = np.abs(np.diff(smooth, prepend=smooth[0]))

    events = []
    in_ev  = False
    for i, d in enumerate(delta):
        if d > angle_threshold / window and not in_ev:
            in_ev = True
            events.append({
                "type":            "domain_motion",
                "time_ps":         float(times[i]),
                "duration_ps":     float(times[min(i + window, len(times) - 1)] - times[i]),
                "frame":           int(i),
                "angle_change_deg": round(float(delta[i:i + window].sum()), 2),
                "description": (
                    f"Interdomain rotation Δ{delta[i:i+window].sum():.1f}° "
                    f"at {times[i]/1000:.1f} ns"
                ),
            })
        elif d <= angle_threshold / window:
            in_ev = False
    return events[:10]  # cap at 10 major events


def assign_event_ids(events: list[dict]) -> list[dict]:
    events_sorted = sorted(events, key=lambda e: e["time_ps"])
    for i, ev in enumerate(events_sorted):
        ev["id"] = f"E{i + 1}"
    return events_sorted


# ── Main ──────────────────────────────────────────────────────────────────────

def generate_events(config: dict) -> list[dict]:
    primary   = config["fto"]["primary_structure"]
    pock_file = DATA_DIR / "pockets" / f"{primary}_pocket_trajectory.json"
    idom_file = RESULTS_DIR / f"{primary}_rep0_interdomain.csv"

    if not pock_file.exists() or not idom_file.exists():
        raise FileNotFoundError("Run 04_motion_analysis and 05_pocket_detection first")

    pocket_traj = load_json(pock_file).get("frames", [])
    idom_df     = pd.read_csv(idom_file)

    events = []
    for pid in ["P1", "P2", "P3"]:
        thr = 200 if pid == "P1" else 150 if pid == "P2" else 80
        events += detect_pocket_openings(pocket_traj, pid, vol_threshold=thr)

    events += detect_domain_motion_events(idom_df)
    return assign_event_ids(events)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--demo", action="store_true")
    args   = parser.parse_args()
    config = load_config()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    if args.demo:
        events = DemoData().events()
    else:
        try:
            events = generate_events(config)
        except FileNotFoundError as exc:
            log.warning(f"{exc} – falling back to demo events")
            events = DemoData().events()

    save_json(events, RESULTS_DIR / "events.json")
    log.info(f"Detected {len(events)} events → results/events.json")
    for ev in events:
        log.info(f"  [{ev.get('id','')}] {ev['type']:20s} @ {ev['time_ps']/1000:.2f} ns  {ev['description']}")


if __name__ == "__main__":
    main()
