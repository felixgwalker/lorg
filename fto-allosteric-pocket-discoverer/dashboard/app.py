"""
FTO Interdomain Pocket Discovery Platform
Main Streamlit dashboard entry point.

Launch:
    streamlit run dashboard/app.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import streamlit as st

# ── Path setup (works whether launched from project root or dashboard/) ───────
_HERE = Path(__file__).parent
_ROOT = _HERE.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

# ── Page config (must be first Streamlit call) ────────────────────────────────
st.set_page_config(
    page_title="FTO Pocket Discovery",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="collapsed",
)

# ── Dark theme CSS ────────────────────────────────────────────────────────────
st.markdown(
    """
    <style>
    /* Base */
    .stApp, section[data-testid="stMain"],
    .block-container { background-color: #0d1117 !important; }

    /* Sidebar */
    [data-testid="stSidebar"] { background-color: #0d1117 !important; }

    /* Text */
    p, span, label, div { color: #c9d1d9; }
    h1, h2, h3, h4      { color: #f0f6fc; }

    /* Tabs */
    .stTabs [data-baseweb="tab-list"] {
        gap: 2px;
        background-color: #161b22;
        border-radius: 10px;
        padding: 4px;
        border: 1px solid #30363d;
    }
    .stTabs [data-baseweb="tab"] {
        background-color: transparent;
        border-radius: 7px;
        padding: 6px 18px;
        color: #8b949e;
        font-size: 14px;
        font-weight: 500;
        border: none;
    }
    .stTabs [aria-selected="true"] {
        background-color: #0d1117 !important;
        color: #58a6ff !important;
        border: 1px solid #30363d;
    }
    .stTabs [data-baseweb="tab-panel"] {
        padding-top: 20px;
    }

    /* Metrics */
    [data-testid="stMetric"] {
        background-color: #161b22;
        border: 1px solid #30363d;
        border-radius: 8px;
        padding: 12px 16px;
    }
    [data-testid="stMetricLabel"]  { color: #8b949e !important; font-size: 12px; }
    [data-testid="stMetricValue"]  { color: #f0f6fc !important; }
    [data-testid="stMetricDelta"]  { color: #3fb950 !important; }

    /* DataFrames */
    .stDataFrame, [data-testid="stDataFrame"] {
        border: 1px solid #30363d;
        border-radius: 8px;
        overflow: hidden;
    }

    /* Buttons */
    .stButton > button {
        background-color: #161b22;
        border: 1px solid #30363d;
        color: #f0f6fc;
        border-radius: 6px;
    }
    .stButton > button:hover {
        background-color: #21262d;
        border-color: #58a6ff;
        color: #58a6ff;
    }

    /* Download buttons */
    .stDownloadButton > button {
        background-color: #161b22;
        border: 1px solid #30363d;
        color: #58a6ff;
        border-radius: 6px;
    }

    /* Select / multiselect */
    .stSelectbox > div > div,
    .stMultiSelect > div > div {
        background-color: #161b22 !important;
        border-color: #30363d !important;
        color: #f0f6fc !important;
    }

    /* Sliders */
    .stSlider > div > div > div { background-color: #58a6ff; }

    /* Expanders */
    [data-testid="stExpander"] {
        background-color: #161b22;
        border: 1px solid #30363d;
        border-radius: 8px;
    }

    /* Code blocks */
    code, pre { background-color: #161b22 !important; }

    /* Radio */
    .stRadio > div { gap: 6px; }
    .stRadio label { color: #c9d1d9 !important; }

    /* Status badge */
    .fto-badge {
        display: inline-block;
        background: #161b22;
        border: 1px solid #30363d;
        border-radius: 12px;
        padding: 2px 10px;
        font-size: 11px;
        color: #8b949e;
        font-family: monospace;
    }
    .fto-badge.ok     { border-color: #3fb950; color: #3fb950; }
    .fto-badge.warn   { border-color: #d29922; color: #d29922; }
    .fto-badge.demo   { border-color: #58a6ff; color: #58a6ff; }

    /* Intro card */
    .fto-intro-card {
        background: linear-gradient(135deg, #161b22 0%, #0d1117 100%);
        border: 1px solid #30363d;
        border-radius: 12px;
        padding: 20px 24px;
        margin-bottom: 16px;
    }
    .fto-pill {
        display: inline-block;
        background: #21262d;
        border: 1px solid #30363d;
        border-radius: 20px;
        padding: 3px 12px;
        font-size: 12px;
        color: #8b949e;
        margin: 2px 4px 2px 0;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

# ── Imports (after path setup) ────────────────────────────────────────────────
from scripts.utils import DATA_DIR, RESULTS_DIR, load_config

from dashboard.tabs import (
    tab_structure,
    tab_trajectory,
    tab_pockets,
    tab_analysis,
    tab_events,
    tab_report,
)


# ── Data-readiness check ──────────────────────────────────────────────────────
def _status() -> dict[str, bool]:
    primary = load_config()["fto"]["primary_structure"]
    return {
        "structures":   (DATA_DIR / "structures" / f"{primary}.pdb").exists(),
        "trajectories": any((DATA_DIR / "trajectories").glob("*.dcd")),
        "pockets":      (
            (DATA_DIR / "pockets" / f"{primary}_pockets.json").exists() or
            (RESULTS_DIR / "pocket_scores.json").exists()
        ),
        "analysis":     (RESULTS_DIR / f"{primary}_rep0_rmsd.csv").exists(),
    }


status    = _status()
demo_mode = not any(status.values())

# ── Header ────────────────────────────────────────────────────────────────────
header_col, badge_col = st.columns([8, 2])
with header_col:
    st.markdown(
        """
        <div style="padding:4px 0 16px 0">
            <h1 style="
                font-size:28px; font-weight:700; color:#f0f6fc;
                margin:0; letter-spacing:-0.02em;
            ">
                🔬 FTO Interdomain Pocket Discovery
            </h1>
            <p style="color:#8b949e;margin:4px 0 0 0;font-size:14px">
                Human FTO (Q9C0B1) · AlkB-family m6A RNA demethylase ·
                Allosteric pocket detection via molecular dynamics simulation
            </p>
        </div>
        """,
        unsafe_allow_html=True,
    )

with badge_col:
    def _chip(label: str, ok: bool) -> str:
        cls  = "ok" if ok else "warn"
        icon = "✓" if ok else "○"
        return f'<span class="fto-badge {cls}">{icon} {label}</span>'

    if demo_mode:
        st.markdown(
            '<div style="text-align:right;margin-top:14px">'
            '<span class="fto-badge demo">◈ DEMO MODE</span>'
            "</div>",
            unsafe_allow_html=True,
        )
    else:
        chips = " ".join([
            _chip("structs",  status["structures"]),
            _chip("MD",       status["trajectories"]),
            _chip("pockets",  status["pockets"]),
            _chip("analysis", status["analysis"]),
        ])
        st.markdown(
            f'<div style="text-align:right;margin-top:16px">{chips}</div>',
            unsafe_allow_html=True,
        )

# ── Platform introduction ─────────────────────────────────────────────────────
with st.expander("ℹ️ What is this platform and how does it work?", expanded=demo_mode):
    col_left, col_right = st.columns([3, 2])
    with col_left:
        st.markdown(
            """
            #### FTO and the m6A demethylation problem

            **FTO** (Fat mass and Obesity-associated protein) is an enzyme that
            removes a chemical tag called **N6-methyladenosine (m6A)** from
            messenger RNA. m6A is the most abundant internal modification in
            mammalian mRNA and plays a central role in regulating how genes
            are expressed. FTO's demethylase activity makes it a validated
            target in obesity, type 2 diabetes, and multiple cancers.

            #### Why this platform?

            Existing FTO inhibitors target the **active site** — the central
            catalytic cavity where the m6A substrate binds. This creates a
            selectivity challenge: the active site is highly conserved across
            the AlkB enzyme family, and inhibitors may hit off-target proteins.

            This platform uses **molecular dynamics (MD) simulation** to search
            for **transient allosteric pockets** — cavities that appear only
            during protein motion and are not visible in static crystal
            structures. A drug targeting such a pocket could modulate FTO
            activity with greater selectivity.

            #### What the platform does

            1. **Acquires** experimental FTO crystal structures from the PDB
            2. **Prepares** and **simulates** the protein using OpenMM (AMBER ff14SB)
            3. **Detects** transient cavities in every simulation frame
            4. **Ranks** pockets by persistence, volume, druggability, and conservation
            5. **Identifies** the conformational events (domain rotations, contact
               changes) that gate pocket opening
            6. **Generates** a pharmacophore model for virtual screening
            """
        )
    with col_right:
        st.markdown(
            """
            #### Platform at a glance

            <div class="fto-intro-card">
            <div style="color:#58a6ff;font-weight:700;font-size:13px;
                        text-transform:uppercase;letter-spacing:.05em">Target</div>
            <div style="color:#f0f6fc;font-size:16px;font-weight:600;margin:4px 0 12px 0">
                Human FTO (UniProt Q9C0B1)
            </div>
            <div style="color:#58a6ff;font-weight:700;font-size:13px;
                        text-transform:uppercase;letter-spacing:.05em">Simulation</div>
            <div style="color:#f0f6fc;font-size:16px;font-weight:600;margin:4px 0 12px 0">
                3 × 20 ns NPT · 310 K · 1 atm
            </div>
            <div style="color:#58a6ff;font-weight:700;font-size:13px;
                        text-transform:uppercase;letter-spacing:.05em">Top finding</div>
            <div style="color:#f0f6fc;font-size:16px;font-weight:600;margin:4px 0 12px 0">
                P1 interdomain pocket · score 0.82 · persistence 38%
            </div>
            <div style="color:#58a6ff;font-weight:700;font-size:13px;
                        text-transform:uppercase;letter-spacing:.05em">Navigate</div>
            <div style="margin:6px 0">
                <span class="fto-pill">🏗️ Structure</span>
                <span class="fto-pill">📈 Trajectory</span>
                <span class="fto-pill">🕳️ Pockets</span>
                <span class="fto-pill">🔬 Analysis</span>
                <span class="fto-pill">⚡ Events</span>
                <span class="fto-pill">📄 Report</span>
            </div>
            </div>
            """,
            unsafe_allow_html=True,
        )

# ── Demo mode banner ──────────────────────────────────────────────────────────
if demo_mode:
    st.markdown(
        """
        <div style="
            background:#161b22; border:1px solid #58a6ff;
            border-radius:8px; padding:14px 18px; margin-bottom:16px;
        ">
            <div style="display:flex;align-items:flex-start;gap:12px">
                <span style="font-size:20px;flex-shrink:0">ℹ️</span>
                <div>
                    <span style="color:#58a6ff;font-weight:600;font-size:14px">
                        Demo mode — no pipeline results found
                    </span>
                    <p style="color:#8b949e;font-size:13px;margin:4px 0 0 0;line-height:1.6">
                        All visualisations are showing <b>scientifically plausible synthetic data</b>
                        so you can explore the full platform without running the pipeline.
                        To populate real results:
                    </p>
                    <div style="margin-top:8px;font-size:13px;color:#8b949e">
                        Quick demo (~60 s):
                        <code style="color:#f0f6fc;background:#0d1117;padding:2px 6px;
                        border-radius:3px">python run_pipeline.py --demo</code>
                        &nbsp;&nbsp;
                        Full pipeline (requires OpenMM, ~2 h):
                        <code style="color:#f0f6fc;background:#0d1117;padding:2px 6px;
                        border-radius:3px">python run_pipeline.py</code>
                    </div>
                </div>
            </div>
        </div>
        """,
        unsafe_allow_html=True,
    )

# ── Main tabs ─────────────────────────────────────────────────────────────────
tabs = st.tabs([
    "🏗️ Structure",
    "📈 Trajectory",
    "🕳️ Pockets",
    "🔬 Analysis",
    "⚡ Events",
    "📄 Report",
])

with tabs[0]:
    tab_structure.render()

with tabs[1]:
    tab_trajectory.render()

with tabs[2]:
    tab_pockets.render()

with tabs[3]:
    tab_analysis.render()

with tabs[4]:
    tab_events.render()

with tabs[5]:
    tab_report.render()

# ── Footer ────────────────────────────────────────────────────────────────────
st.markdown(
    """
    <div style="
        margin-top:40px; padding:16px 0;
        border-top:1px solid #21262d;
        display:flex; justify-content:space-between; align-items:center;
        color:#8b949e; font-size:11px; flex-wrap:wrap; gap:8px;
    ">
        <span>FTO Pocket Discovery Platform · v1.0</span>
        <span>
            OpenMM · MDAnalysis · 3Dmol.js · Streamlit · Plotly
        </span>
        <span>
            <a href="https://www.rcsb.org/structure/3LFM"
               style="color:#8b949e" target="_blank">PDB 3LFM</a> ·
            <a href="https://www.rcsb.org/structure/4IE4"
               style="color:#8b949e" target="_blank">4IE4</a> ·
            <a href="https://www.rcsb.org/structure/4ZS3"
               style="color:#8b949e" target="_blank">4ZS3</a> ·
            <a href="https://www.uniprot.org/uniprot/Q9C0B1"
               style="color:#8b949e" target="_blank">UniProt Q9C0B1</a>
        </span>
    </div>
    """,
    unsafe_allow_html=True,
)
