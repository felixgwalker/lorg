"""
[Analysis] tab – PCA, clustering, conservation overlay, RMSD multi-replica.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import streamlit as st

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from scripts.utils import DATA_DIR, RESULTS_DIR, DemoData, load_config

CONSERVATION_DIR = DATA_DIR / "conservation"

CLUSTER_COLORS = ["#58a6ff", "#FF6B35", "#50C878", "#FFD700", "#c778dd"]


def _plotly_theme() -> dict:
    return dict(
        paper_bgcolor="#0d1117", plot_bgcolor="#161b22",
        font=dict(color="#f0f6fc", size=12),
        xaxis=dict(gridcolor="#21262d", linecolor="#30363d", zerolinecolor="#30363d"),
        yaxis=dict(gridcolor="#21262d", linecolor="#30363d", zerolinecolor="#30363d"),
        margin=dict(l=50, r=20, t=45, b=50),
    )


def _load_pca() -> dict:
    primary = load_config()["fto"]["primary_structure"]
    path    = RESULTS_DIR / f"{primary}_rep0_pca.json"
    if path.exists():
        with open(path) as f:
            return json.load(f)
    return DemoData().pca()


def _load_rmsf() -> pd.DataFrame:
    primary = load_config()["fto"]["primary_structure"]
    path    = RESULTS_DIR / f"{primary}_rep0_rmsf.csv"
    if path.exists():
        return pd.read_csv(path)
    return DemoData().rmsf()


def _load_conservation() -> dict:
    path = CONSERVATION_DIR / "fto_conservation.json"
    if path.exists():
        with open(path) as f:
            return json.load(f)
    return DemoData().conservation_scores()


def render() -> None:
    config  = load_config()
    primary = config["fto"]["primary_structure"]
    dd      = DemoData()

    st.markdown(
        """
        <h2 style="color:#f0f6fc;margin:0 0 4px 0">Analysis</h2>
        <p style="color:#8b949e;font-size:14px;margin:0 0 20px 0">
        Conformational landscape · PCA · clustering · evolutionary conservation
        </p>
        """,
        unsafe_allow_html=True,
    )

    pca_data  = _load_pca()
    rmsf_df   = _load_rmsf()
    cons_data = _load_conservation()

    # ── PCA ───────────────────────────────────────────────────────────────────
    st.markdown("#### Conformational Landscape (PCA)")

    with st.expander("ℹ️ What is PCA and how should I read this plot?", expanded=False):
        st.markdown(
            """
            **Principal Component Analysis (PCA)** reduces the enormous complexity
            of a molecular dynamics trajectory — millions of atomic coordinates —
            down to a small number of dominant motions.

            Each point on the plot is one simulation frame (snapshot). Points that
            cluster together represent similar protein conformations. Points far apart
            represent structurally distinct states.

            - **PC1 (x-axis):** The single largest motion in the simulation — typically
              a large-scale domain rearrangement or hinge motion.
            - **PC2 (y-axis):** The second-largest independent motion.
            - **Explained variance:** How much of the total structural variation is
              captured by each component. PC1 + PC2 together often capture 40–60 %
              of all motion.

            **What to look for:**
            - Distinct, well-separated clusters indicate metastable conformational
              states (e.g., open vs closed domain orientations).
            - Trajectories that drift continuously suggest the protein is still
              exploring new conformations; those that revisit the same regions
              indicate convergence.

            **Colour by Cluster** groups frames by k-means clustering.
            **Colour by Time** shows how the simulation progresses through
            conformational space.
            """
        )

    c_pca, c_ev = st.columns([3, 2])

    with c_pca:
        pc1  = np.array(pca_data["pc1"])
        pc2  = np.array(pca_data["pc2"])
        clus = np.array(pca_data["clusters"])
        t_ns = np.array(pca_data["time_ps"]) / 1_000

        colour_by = st.radio(
            "Colour by", ["Cluster", "Time"], horizontal=True, key="pca_colour",
        )
        if colour_by == "Cluster":
            n_clus = int(clus.max()) + 1
            fig_pca = go.Figure()
            for c in range(n_clus):
                mask = clus == c
                fig_pca.add_trace(go.Scatter(
                    x=pc1[mask], y=pc2[mask],
                    mode="markers",
                    marker=dict(color=CLUSTER_COLORS[c % len(CLUSTER_COLORS)],
                                size=3, opacity=0.7),
                    name=f"Cluster C{c}",
                ))
        else:
            fig_pca = go.Figure(go.Scatter(
                x=pc1, y=pc2, mode="markers",
                marker=dict(
                    color=t_ns, colorscale="Viridis", size=3, opacity=0.7,
                    colorbar=dict(title="Time (ns)", thickness=12, len=0.7),
                ),
                name="Frame",
            ))

        ev = pca_data["explained_variance"]
        fig_pca.update_layout(
            **_plotly_theme(), height=380,
            xaxis_title=f"PC1 — {ev[0]:.1f}% of variance",
            yaxis_title=f"PC2 — {ev[1]:.1f}% of variance",
            title=dict(text="Conformational Landscape", font_color="#8b949e"),
            legend=dict(bgcolor="rgba(0,0,0,0)", font_size=11),
        )
        st.plotly_chart(fig_pca, use_container_width=True)

    with c_ev:
        n_comp = min(10, len(ev))
        fig_ev = go.Figure(go.Bar(
            x=[f"PC{i+1}" for i in range(n_comp)],
            y=ev[:n_comp],
            marker_color=["#58a6ff" if i < 2 else "#30363d" for i in range(n_comp)],
            text=[f"{v:.1f}%" for v in ev[:n_comp]],
            textposition="outside",
        ))
        fig_ev.update_layout(
            **_plotly_theme(), height=380,
            yaxis_title="Variance explained (%)",
            title=dict(text="Scree Plot — Variance per Component", font_color="#8b949e"),
            showlegend=False,
        )
        st.plotly_chart(fig_ev, use_container_width=True)
        st.caption(
            "The scree plot shows how much structural information each principal "
            "component captures. Blue bars (PC1, PC2) are plotted in the landscape above."
        )

    # ── Cluster populations ───────────────────────────────────────────────────
    st.markdown("#### Cluster Populations")
    st.caption(
        "K-means clustering groups trajectory frames into distinct conformational states. "
        "Each cluster represents a basin in the protein's energy landscape. "
        "Large clusters are thermodynamically stable states; small clusters may be "
        "transient or high-energy intermediates."
    )
    n_clus  = int(clus.max()) + 1
    counts  = [int(np.sum(clus == c)) for c in range(n_clus)]
    total   = len(clus)
    cols    = st.columns(n_clus)
    for i, (count, col) in enumerate(zip(counts, cols)):
        col.metric(
            f"Cluster C{i}",
            f"{count:,} frames",
            f"{count / total:.1%} of trajectory",
        )

    # ── RMSF + Conservation ───────────────────────────────────────────────────
    st.markdown("#### Flexibility vs Evolutionary Conservation")

    with st.expander("ℹ️ How to read the conservation overlay", expanded=False):
        st.markdown(
            """
            This plot overlays two independent data sources on the same residue axis:

            **Red line — RMSF:** Per-residue flexibility from the MD simulation.
            High peaks = mobile regions. Low regions = rigid regions.

            **Coloured dots — ConSurf conservation score (1–9):**
            How well conserved each residue is across ~50 vertebrate FTO orthologues.
            - Score **9** (dark orange) = invariant across species → functionally critical
            - Score **1** (blue) = variable → less essential

            **The key insight:** Residues that are *both* rigid in MD *and* highly
            conserved across evolution are likely to be structurally and functionally
            important. The active-site residues (H231, D233, H307) score 9/9 and show
            low RMSF — consistent with their essential catalytic role.

            Pocket-lining residues that are also highly conserved suggest the pocket
            geometry is evolutionarily maintained, strengthening the case for its
            functional significance.
            """
        )

    cat   = config["fto"]["domains"]["catalytic"]
    cterm = config["fto"]["domains"]["cterminal"]
    cons_residues = cons_data.get("residues", [])
    cons_scores   = cons_data.get("scores", [])

    fig_cons = go.Figure()
    fig_cons.add_vrect(
        x0=cat["start"], x1=cat["end"],
        fillcolor="#4A90D9", opacity=0.06,
        annotation_text="Catalytic", annotation_position="top left",
        annotation_font_color="#4A90D9",
    )
    fig_cons.add_vrect(
        x0=cterm["start"], x1=cterm["end"],
        fillcolor="#50C878", opacity=0.06,
        annotation_text="C-terminal", annotation_position="top left",
        annotation_font_color="#50C878",
    )
    fig_cons.add_trace(go.Scatter(
        x=rmsf_df["residue"], y=rmsf_df["rmsf_A"],
        mode="lines", line=dict(color="#ff7b72", width=1.2),
        name="RMSF (Å)",
    ))
    if cons_residues:
        fig_cons.add_trace(go.Scatter(
            x=cons_residues, y=cons_scores,
            mode="markers",
            marker=dict(
                symbol="circle",
                color=cons_scores,
                colorscale=[[0, "#00CCFF"], [0.5, "#FFFF99"], [1, "#FF6600"]],
                size=4, opacity=0.85,
                cmin=1, cmax=9,
                colorbar=dict(
                    title="ConSurf<br>score",
                    thickness=12, len=0.6, y=0.5,
                ),
            ),
            name="Conservation (1–9)",
            yaxis="y2",
        ))

    _theme = {k: v for k, v in _plotly_theme().items() if k != "yaxis"}
    fig_cons.update_layout(
        **_theme, height=320,
        xaxis_title="Residue number",
        yaxis=dict(title="RMSF (Å)", gridcolor="#21262d"),
        yaxis2=dict(
            title="Conservation score (1–9)",
            overlaying="y", side="right",
            range=[1, 9], gridcolor="rgba(0,0,0,0)",
        ),
        legend=dict(bgcolor="rgba(0,0,0,0)"),
    )
    st.plotly_chart(fig_cons, use_container_width=True)

    # ── Conservation heatmap ──────────────────────────────────────────────────
    with st.expander("Full-sequence conservation heatmap"):
        st.caption(
            "Each cell shows the ConSurf conservation score for one residue. "
            "Blue = variable (low conservation), orange = invariant (high conservation). "
            "Rows represent blocks of 50 consecutive residues."
        )
        if cons_residues:
            width  = 50
            n      = len(cons_residues)
            n_rows = (n + width - 1) // width
            padded = cons_scores + [0] * (n_rows * width - n)
            mat    = np.array(padded).reshape(n_rows, width)

            fig_hm = go.Figure(go.Heatmap(
                z=mat,
                colorscale=[
                    [0,   "#00CCFF"], [0.25, "#66EEFF"],
                    [0.5, "#FFFF99"], [0.75, "#FF9933"],
                    [1.0, "#FF6600"],
                ],
                zmin=1, zmax=9,
                colorbar=dict(title="ConSurf<br>score (1–9)"),
                xgap=0.5, ygap=0.5,
            ))
            fig_hm.update_layout(
                **_plotly_theme(), height=max(200, n_rows * 12),
                xaxis_title="Position within block (residues)",
                yaxis_title="Block (50 residues each)",
                title=dict(
                    text="Evolutionary Conservation — 1 (variable) → 9 (invariant)",
                    font_color="#8b949e",
                ),
            )
            st.plotly_chart(fig_hm, use_container_width=True)
        else:
            st.info("Conservation data not available. Run `python scripts/07_conservation.py`.")
