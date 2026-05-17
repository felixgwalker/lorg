"""
[Pockets] tab – transient pocket explorer, ranking, volumes.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from scripts.utils import DATA_DIR, RESULTS_DIR, DemoData, load_config
from dashboard.components.viewer import render_pocket_sphere

STRUCTURES_DIR = DATA_DIR / "structures"
POCKETS_DIR    = DATA_DIR / "pockets"


def _load_pockets() -> list[dict]:
    for p in [
        RESULTS_DIR / "pocket_scores.json",
        POCKETS_DIR / "3LFM_pockets.json",
    ]:
        if p.exists():
            with open(p) as f:
                return json.load(f)
    return DemoData().pockets()


def _plotly_theme() -> dict:
    return dict(
        paper_bgcolor="#0d1117", plot_bgcolor="#161b22",
        font=dict(color="#f0f6fc", size=12),
        xaxis=dict(gridcolor="#21262d", linecolor="#30363d"),
        yaxis=dict(gridcolor="#21262d", linecolor="#30363d"),
        margin=dict(l=50, r=20, t=40, b=50),
    )


def _hex_to_rgba(hex_color: str, alpha: float = 0.10) -> str:
    """Convert a #RRGGBB hex color to rgba(...) string."""
    h = hex_color.lstrip("#")
    if len(h) == 6:
        r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
        return f"rgba({r},{g},{b},{alpha})"
    return f"rgba(255,107,53,{alpha})"


def render() -> None:
    config  = load_config()
    primary = config["fto"]["primary_structure"]
    pockets = _load_pockets()
    dd      = DemoData()

    st.markdown(
        """
        <h2 style="color:#f0f6fc;margin:0 0 4px 0">Pocket Explorer</h2>
        <p style="color:#8b949e;font-size:14px;margin:0 0 20px 0">
        Transient cavity detection · composite ranking · volume dynamics
        </p>
        """,
        unsafe_allow_html=True,
    )

    with st.expander("ℹ️ How pockets are detected and scored", expanded=False):
        st.markdown(
            """
            **What is a pocket?**
            A pocket (or cavity) is a concave region on the protein surface
            large enough to accommodate a small drug-like molecule. Most drug
            targets are defined by their pockets.

            **Transient pockets** are particularly valuable: they are hidden in
            the static crystal structure but open during molecular dynamics
            simulation, revealing *cryptic* binding sites invisible to traditional
            structure-based drug design.

            **How the composite score is calculated:**

            | Component | Weight | What it means |
            |-----------|--------|---------------|
            | Persistence | 30 % | Fraction of simulation frames in which the pocket is open |
            | Volume (normalised) | 25 % | Mean open-state volume; larger = more room for a ligand |
            | Druggability | 20 % | Hydrophobicity, charge balance, and volume together |
            | Conservation | 15 % | ConSurf evolutionary score of lining residues |
            | Enclosure | 10 % | How buried the pocket is; more enclosed = better binding |

            **Persistence** is the most important single metric: a pocket that
            opens rarely (< 5 % of frames) is unlikely to be pharmacologically useful.
            Pockets persisting in > 20 % of frames are considered viable lead targets.

            **Druggability score (0–1):** A score ≥ 0.7 is considered well-druggable.
            Below 0.4 the pocket is unlikely to bind drug-like compounds tightly.
            """
        )

    # ── Ranking table ─────────────────────────────────────────────────────────
    st.markdown("#### Ranked Pockets")

    df = pd.DataFrame([{
        "Rank":           p.get("rank", i + 1),
        "ID":             p["id"],
        "Name":           p.get("name", ""),
        "Score":          p.get("pocket_score", 0),
        "Persistence":    p.get("persistence", 0),
        "Volume (Å³)":    p.get("mean_volume_A3", 0),
        "Druggability":   p.get("druggability_score", 0),
        "Description":    p.get("description", ""),
    } for i, p in enumerate(pockets)])

    st.dataframe(
        df.drop(columns=["Description"]),
        use_container_width=True,
        hide_index=True,
        column_config={
            "Score":       st.column_config.ProgressColumn(
                "Score", min_value=0, max_value=1, format="%.2f",
                help="Composite druggability score (0–1). ≥ 0.7 = high priority.",
            ),
            "Persistence": st.column_config.ProgressColumn(
                "Persistence", min_value=0, max_value=1, format="%.1%",
                help="Fraction of MD frames in which this pocket is open.",
            ),
            "Druggability": st.column_config.ProgressColumn(
                "Druggability", min_value=0, max_value=1, format="%.2f",
                help="Estimate of ligand-binding suitability based on shape and chemistry.",
            ),
        },
    )

    # ── Pocket selector ───────────────────────────────────────────────────────
    pocket_ids = [p["id"] for p in pockets]
    sel_id     = st.radio("Select pocket to explore", pocket_ids, horizontal=True)
    sel_pocket = next(p for p in pockets if p["id"] == sel_id)
    pock_color = sel_pocket.get("color", "#FF6B35")

    st.markdown(
        f"""
        <div style="
            background:#161b22;border:1px solid #30363d;
            border-left: 4px solid {pock_color};
            border-radius:0 8px 8px 0;padding:16px;margin:12px 0;
        ">
            <b style="color:#f0f6fc;font-size:16px">{sel_pocket['id']} · {sel_pocket.get('name','')}</b><br>
            <span style="color:#8b949e;font-size:13px">{sel_pocket.get('description','')}</span>
        </div>
        """,
        unsafe_allow_html=True,
    )

    # ── Pocket metrics ────────────────────────────────────────────────────────
    m1, m2, m3, m4, m5 = st.columns(5)
    m1.metric("Score",        f"{sel_pocket.get('pocket_score', 0):.2f}",
              help="Composite druggability score (0–1)")
    m2.metric("Persistence",  f"{sel_pocket.get('persistence', 0):.1%}",
              help="Fraction of frames where pocket volume exceeds threshold")
    m3.metric("Max volume",   f"{sel_pocket.get('max_volume_A3', 0):.0f} Å³",
              help="Peak pocket volume across all frames")
    m4.metric("Druggability", f"{sel_pocket.get('druggability_score', 0):.2f}",
              help="Ligand-binding suitability estimate")
    m5.metric("First seen",   f"{sel_pocket.get('first_appearance_ps', 0) / 1000:.1f} ns",
              help="Time at which the pocket first opens in the simulation")

    # ── Volume trajectory ─────────────────────────────────────────────────────
    vol_traj = sel_pocket.get("volume_trajectory")
    if not vol_traj:
        vol_traj = dd.pockets()[0]["volume_trajectory"]

    n   = len(vol_traj)
    t   = np.linspace(0, config["simulation"]["production_ns"], n)
    thr = 100  # Å³ open/closed threshold

    fill_rgba = _hex_to_rgba(pock_color, alpha=0.10)

    fig_vol = go.Figure()
    fig_vol.add_trace(go.Scatter(
        x=t, y=vol_traj,
        mode="lines", line=dict(color=pock_color, width=1.3),
        fill="tozeroy", fillcolor=fill_rgba,
        name="Volume",
    ))
    fig_vol.add_hline(
        y=thr, line_dash="dash", line_color="#FFD700", line_width=1.2,
        annotation_text=f"Open threshold ({thr} Å³)",
        annotation_font_color="#FFD700",
    )
    fig_vol.update_layout(
        **_plotly_theme(), height=260,
        xaxis_title="Time (ns)", yaxis_title="Pocket volume (Å³)",
        title=dict(text=f"Pocket {sel_id} — Volume over Simulation Time", font_color="#8b949e"),
    )
    st.plotly_chart(fig_vol, use_container_width=True)
    st.caption(
        f"The gold dashed line marks the open/closed threshold ({thr} Å³). "
        "Regions above the line represent frames where this pocket is considered open "
        "and accessible to a ligand. Persistence is the fraction of time spent above this threshold."
    )

    # ── 3D viewer + residue composition ──────────────────────────────────────
    c_viewer, c_chart = st.columns([3, 2])

    with c_chart:
        st.markdown("**Chemical character of pocket lining**")
        st.caption(
            "The pocket's chemical environment determines what types of ligands can bind. "
            "Hydrophobic pockets prefer lipophilic ligands; charged pockets suit ionic interactions."
        )
        hydro   = sel_pocket.get("hydrophobic_fraction", 0.5)
        charged = sel_pocket.get("charged_fraction", 0.2)
        polar   = max(0, 1 - hydro - charged)
        prop_df = pd.DataFrame({
            "Type":     ["Hydrophobic", "Charged", "Polar / other"],
            "Fraction": [hydro, charged, polar],
        })
        fig_pie = px.pie(
            prop_df, names="Type", values="Fraction",
            color_discrete_sequence=["#4A90D9", "#FF6B35", "#50C878"],
            hole=0.45,
        )
        fig_pie.update_layout(
            paper_bgcolor="#0d1117",
            font_color="#f0f6fc",
            margin=dict(l=10, r=10, t=10, b=10),
            showlegend=True,
            legend=dict(bgcolor="rgba(0,0,0,0)", font_size=11),
            height=220,
        )
        st.plotly_chart(fig_pie, use_container_width=True)

        st.markdown("**Lining residues**")
        st.caption("Residue numbers that line the pocket interior.")
        res_list = sel_pocket.get("residues", [])
        if res_list:
            st.code(", ".join(map(str, sorted(res_list))), language=None)
        else:
            st.caption("No residue data available.")

    with c_viewer:
        pdb_path = STRUCTURES_DIR / f"{primary}.pdb"
        if pdb_path.exists():
            cat   = config["fto"]["domains"]["catalytic"]
            cterm = config["fto"]["domains"]["cterminal"]
            render_pocket_sphere(
                pdb_path.read_text(encoding="utf-8", errors="ignore"),
                sel_pocket,
                height=460,
                cat_range=(cat["start"], cat["end"]),
                cterm_range=(cterm["start"], cterm["end"]),
            )
            st.caption(
                "Orange residues = pocket lining. The semi-transparent sphere "
                "indicates the approximate pocket volume and position."
            )
        else:
            st.info(
                "3D pocket view requires the PDB structure. "
                "Run `python scripts/01_acquire.py` to download it."
            )

    # ── Persistence comparison ────────────────────────────────────────────────
    st.markdown("#### Pocket Persistence Comparison")
    st.caption(
        "A direct comparison of how often each pocket is open across the simulation. "
        "Higher persistence means greater opportunity for a ligand to bind."
    )
    fig_bar = go.Figure()
    for p in pockets:
        fig_bar.add_trace(go.Bar(
            name=p["id"],
            x=[p.get("name", p["id"])],
            y=[p.get("persistence", 0)],
            marker_color=p.get("color", "#58a6ff"),
            text=[f"{p.get('persistence', 0):.1%}"],
            textposition="outside",
        ))
    fig_bar.update_layout(
        **_plotly_theme(), height=280,
        yaxis_title="Persistence (fraction of frames open)",
        yaxis_tickformat=".0%",
        yaxis_range=[0, 1.15],
        barmode="group",
        showlegend=False,
        title=dict(text="Fraction of Simulation Frames Where Each Pocket Is Open", font_color="#8b949e"),
    )
    st.plotly_chart(fig_bar, use_container_width=True)
