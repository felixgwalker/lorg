"""
[Trajectory] tab – RMSD, interdomain motion, frame playback.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from scripts.utils import DATA_DIR, RESULTS_DIR, DemoData, load_config
from dashboard.components.viewer import render_structure

STRUCTURES_DIR = DATA_DIR / "structures"


def _load_or_demo(path: Path, loader, demo_fn):
    if path.exists():
        return loader(path)
    return demo_fn()


def _plotly_theme() -> dict:
    return dict(
        paper_bgcolor="#0d1117",
        plot_bgcolor="#161b22",
        font=dict(color="#f0f6fc", size=12),
        xaxis=dict(gridcolor="#21262d", linecolor="#30363d", zerolinecolor="#30363d"),
        yaxis=dict(gridcolor="#21262d", linecolor="#30363d", zerolinecolor="#30363d"),
        margin=dict(l=50, r=20, t=40, b=50),
    )


def render() -> None:
    config  = load_config()
    primary = config["fto"]["primary_structure"]
    n_rep   = config["simulation"]["n_replicas"]
    dd      = DemoData()

    st.markdown(
        """
        <h2 style="color:#f0f6fc;margin:0 0 4px 0">Trajectory Analysis</h2>
        <p style="color:#8b949e;font-size:14px;margin:0 0 20px 0">
        Structural dynamics across the MD simulation · RMSD · RMSF · Interdomain motion
        </p>
        """,
        unsafe_allow_html=True,
    )

    with st.expander("ℹ️ How to read these plots", expanded=False):
        st.markdown(
            """
            **RMSD (Root Mean Square Deviation)** measures how far the protein has
            moved from its starting structure over time (in Ångströms, Å).
            A rising RMSD indicates the protein is sampling new conformations.
            Plateau behaviour suggests the simulation has reached equilibrium.
            Values of 1–3 Å are typical for a well-folded domain.

            **RMSF (Root Mean Square Fluctuation)** shows the *average* flexibility
            of each residue over the whole simulation. High RMSF = mobile; low RMSF
            = rigid. The active site (yellow markers) is expected to be rigid.
            The interdomain linker region is expected to be more flexible.

            **Interdomain motion** tracks how the catalytic and C-terminal domains
            move relative to each other. Large angle changes indicate domain
            opening events that expose the P1 allosteric pocket.
            """
        )

    # ── Replica selector (must come before data loading) ─────────────────────
    rep = st.select_slider(
        "Replica",
        options=list(range(n_rep)),
        value=0,
        help="Each replica is an independent MD simulation run from the same starting structure. "
             "Consistent findings across replicas increase confidence.",
    )

    prefix = f"{primary}_rep{rep}"

    # ── Load data (falls back to demo if replica file absent) ─────────────────
    rmsd_path  = RESULTS_DIR / f"{prefix}_rmsd.csv"
    rmsf_path  = RESULTS_DIR / f"{prefix}_rmsf.csv"
    idom_path  = RESULTS_DIR / f"{prefix}_interdomain.csv"

    using_fallback = not rmsd_path.exists() and rep > 0 and (RESULTS_DIR / f"{primary}_rep0_rmsd.csv").exists()

    rmsd_df = _load_or_demo(rmsd_path, pd.read_csv, dd.rmsd)
    rmsf_df = _load_or_demo(rmsf_path, pd.read_csv, dd.rmsf)
    idom_df = _load_or_demo(idom_path, pd.read_csv, dd.interdomain_motion)

    if using_fallback:
        st.info(
            f"Results for replica {rep} not found — showing replica 0 data. "
            "Run the full pipeline to generate all replicas."
        )

    # Convert time to ns for display
    rmsd_df = rmsd_df.copy()
    idom_df = idom_df.copy()
    rmsd_df["time_ns"] = rmsd_df["time_ps"] / 1_000
    idom_df["time_ns"] = idom_df["time_ps"] / 1_000

    # ── Top metrics ───────────────────────────────────────────────────────────
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Simulation length", f"{config['simulation']['production_ns']:.0f} ns")
    m2.metric("Mean RMSD",         f"{rmsd_df['rmsd_A'].mean():.2f} Å")
    m3.metric("Peak domain angle", f"{idom_df['interdomain_angle_deg'].max():.1f}°")
    m4.metric("Frames analysed",   f"{len(rmsd_df):,}")

    # ── RMSD plot ─────────────────────────────────────────────────────────────
    st.markdown("#### Backbone RMSD")
    fig_rmsd = go.Figure()
    fig_rmsd.add_trace(go.Scatter(
        x=rmsd_df["time_ns"], y=rmsd_df["rmsd_A"],
        mode="lines", line=dict(color="#58a6ff", width=1.5),
        name="Backbone RMSD", fill="tozeroy",
        fillcolor="rgba(88,166,255,0.08)",
    ))
    w = max(1, len(rmsd_df) // 50)
    smoothed = rmsd_df["rmsd_A"].rolling(w, center=True).mean()
    fig_rmsd.add_trace(go.Scatter(
        x=rmsd_df["time_ns"], y=smoothed,
        mode="lines", line=dict(color="#FFD700", width=2),
        name=f"Rolling average ({w} frames)",
    ))
    fig_rmsd.update_layout(
        **_plotly_theme(),
        height=280,
        xaxis_title="Time (ns)",
        yaxis_title="RMSD (Å)",
        showlegend=True,
        legend=dict(bgcolor="rgba(0,0,0,0)"),
    )
    st.plotly_chart(fig_rmsd, use_container_width=True)
    st.caption(
        "Backbone RMSD relative to the starting structure. The yellow line is a "
        "rolling average that filters out fast thermal noise and reveals slower "
        "conformational changes."
    )

    # ── RMSF plot ─────────────────────────────────────────────────────────────
    st.markdown("#### Per-residue Flexibility (RMSF)")
    cat   = config["fto"]["domains"]["catalytic"]
    cterm = config["fto"]["domains"]["cterminal"]

    fig_rmsf = go.Figure()
    fig_rmsf.add_vrect(
        x0=cat["start"], x1=cat["end"],
        fillcolor="#4A90D9", opacity=0.06,
        annotation_text="Catalytic domain", annotation_position="top left",
        annotation_font_color="#4A90D9",
    )
    fig_rmsf.add_vrect(
        x0=cterm["start"], x1=cterm["end"],
        fillcolor="#50C878", opacity=0.06,
        annotation_text="C-terminal domain", annotation_position="top left",
        annotation_font_color="#50C878",
    )
    fig_rmsf.add_trace(go.Scatter(
        x=rmsf_df["residue"], y=rmsf_df["rmsf_A"],
        mode="lines", line=dict(color="#ff7b72", width=1.2),
        fill="tozeroy", fillcolor="rgba(255,123,114,0.08)",
        name="RMSF",
    ))
    for r in config["fto"]["active_site_residues"]:
        fig_rmsf.add_vline(
            x=r, line_dash="dot", line_color="#FFD700", line_width=1,
            annotation_text=f"{r}", annotation_font_size=9,
            annotation_font_color="#FFD700",
        )
    fig_rmsf.update_layout(
        **_plotly_theme(),
        height=280,
        xaxis_title="Residue number",
        yaxis_title="RMSF (Å)",
        showlegend=False,
    )
    st.plotly_chart(fig_rmsf, use_container_width=True)
    st.caption(
        "Per-residue flexibility. Yellow dotted lines mark active-site residues. "
        "High peaks in the inter-domain linker region (~315–330) indicate the hinge "
        "motion that opens the P1 pocket."
    )

    # ── Interdomain motion ────────────────────────────────────────────────────
    st.markdown("#### Interdomain Motion")
    st.caption(
        "The catalytic and C-terminal domains are not rigidly connected. "
        "These plots track how they move relative to each other. "
        "Large angle changes (> 5°) correlate with transient pocket opening."
    )
    c_left, c_right = st.columns(2)

    with c_left:
        fig_dist = go.Figure()
        fig_dist.add_trace(go.Scatter(
            x=idom_df["time_ns"],
            y=idom_df["interdomain_dist_A"],
            mode="lines", line=dict(color="#50C878", width=1.4),
            name="CoG distance",
        ))
        fig_dist.update_layout(
            **_plotly_theme(), height=240,
            xaxis_title="Time (ns)",
            yaxis_title="Distance (Å)",
            title=dict(text="Domain Centre-of-Geometry Distance", font_color="#8b949e"),
        )
        st.plotly_chart(fig_dist, use_container_width=True)

    with c_right:
        fig_angle = go.Figure()
        fig_angle.add_trace(go.Scatter(
            x=idom_df["time_ns"],
            y=idom_df["interdomain_angle_deg"],
            mode="lines", line=dict(color="#FF6B35", width=1.4),
            fill="tozeroy", fillcolor="rgba(255,107,53,0.08)",
            name="Angle",
        ))
        peak_idx = idom_df["interdomain_angle_deg"].idxmax()
        fig_angle.add_annotation(
            x=idom_df.loc[peak_idx, "time_ns"],
            y=idom_df.loc[peak_idx, "interdomain_angle_deg"],
            text=f"Peak {idom_df.loc[peak_idx,'interdomain_angle_deg']:.1f}°",
            showarrow=True, arrowcolor="#FFD700",
            font=dict(color="#FFD700", size=11),
        )
        fig_angle.update_layout(
            **_plotly_theme(), height=240,
            xaxis_title="Time (ns)",
            yaxis_title="Angle (°)",
            title=dict(text="Interdomain Opening Angle", font_color="#8b949e"),
        )
        st.plotly_chart(fig_angle, use_container_width=True)

    # ── Frame explorer ────────────────────────────────────────────────────────
    with st.expander("Frame explorer — view a specific simulation snapshot"):
        st.caption(
            "Select a time point to inspect the protein conformation at that moment. "
            "This uses the static crystal structure as a placeholder; "
            "actual trajectory frames require running the full pipeline."
        )
        n_frames  = len(rmsd_df)
        frame_idx = st.slider("Frame", 0, n_frames - 1, 0, key="traj_frame")
        t_ns      = float(rmsd_df.iloc[frame_idx]["time_ns"])
        rmsd_val  = float(rmsd_df.iloc[frame_idx]["rmsd_A"])
        angle_val = float(idom_df.iloc[min(frame_idx, len(idom_df) - 1)]["interdomain_angle_deg"])

        fc1, fc2, fc3 = st.columns(3)
        fc1.metric("Time", f"{t_ns:.2f} ns")
        fc2.metric("RMSD", f"{rmsd_val:.2f} Å")
        fc3.metric("Domain angle", f"{angle_val:.1f}°")

        pdb_path = STRUCTURES_DIR / f"{primary}.pdb"
        if pdb_path.exists():
            pdb_content = pdb_path.read_text(encoding="utf-8", errors="ignore")
            render_structure(
                pdb_content, height=420,
                cat_range=(cat["start"], cat["end"]),
                cterm_range=(cterm["start"], cterm["end"]),
                highlight_residues=config["fto"]["active_site_residues"],
                label=f"Frame {frame_idx} · {t_ns:.2f} ns · RMSD {rmsd_val:.2f} Å",
            )
        else:
            st.info(
                "3D frame viewing requires the PDB structure file. "
                "Run `python scripts/01_acquire.py` to download it."
            )
