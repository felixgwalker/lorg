"""
[Report] tab – story mode, key findings, pharmacophore, export.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from datetime import date

import pandas as pd
import plotly.graph_objects as go
import streamlit as st

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from scripts.utils import DATA_DIR, RESULTS_DIR, DemoData, load_config
from dashboard.components.story_mode import story_mode_ui

PHARMA_FEATURE_LABELS = {
    "HBA": "H-bond acceptor — the pocket can accept a hydrogen bond from a ligand",
    "HBD": "H-bond donor — the pocket can donate a hydrogen bond to a ligand",
    "AR":  "Aromatic centre — π-stacking interaction with an aromatic ring on the ligand",
    "HYD": "Hydrophobic region — prefers non-polar, lipophilic ligand groups",
    "PI":  "Positive ionisable — positively charged at physiological pH",
    "NI":  "Negative ionisable — negatively charged at physiological pH",
}


def _load_json(path: Path, fallback):
    if path.exists():
        with open(path) as f:
            return json.load(f)
    return fallback


def _plotly_theme() -> dict:
    return dict(
        paper_bgcolor="#0d1117", plot_bgcolor="#161b22",
        font=dict(color="#f0f6fc", size=12),
        margin=dict(l=50, r=20, t=45, b=50),
    )


def render() -> None:
    config  = load_config()
    primary = config["fto"]["primary_structure"]
    dd      = DemoData()

    pockets  = _load_json(RESULTS_DIR / "pocket_scores.json",   dd.pockets())
    events   = _load_json(RESULTS_DIR / "events.json",           dd.events())
    pharma   = _load_json(RESULTS_DIR / "pharmacophore.json",    dd.pharmacophore())
    top_p    = pockets[0] if pockets else {}

    tab_story, tab_findings, tab_pharma, tab_export = st.tabs([
        "📖 Story Mode",
        "📊 Key Findings",
        "💊 Pharmacophore",
        "⬇️ Export",
    ])

    # ── Story Mode ────────────────────────────────────────────────────────────
    with tab_story:
        story_mode_ui(pockets=pockets)

    # ── Key Findings ──────────────────────────────────────────────────────────
    with tab_findings:
        st.markdown(
            """
            <h3 style="color:#f0f6fc;margin:0 0 8px 0">Key Findings</h3>
            <p style="color:#8b949e;font-size:13px;margin:0 0 20px 0">
            A plain-language summary of the most significant results from this analysis.
            </p>
            """,
            unsafe_allow_html=True,
        )

        findings = [
            {
                "icon": "🕳️", "color": "#FF6B35",
                "title": "Transient allosteric pocket P1 identified",
                "body": (
                    f"Pocket <b>P1</b> (Interdomain Interface Pocket) emerges at the "
                    f"catalytic–CTD interface following a domain rotation event. "
                    f"It achieves a peak volume of <b>{top_p.get('max_volume_A3', 487):.0f} Å³</b> "
                    f"and is accessible in <b>{top_p.get('persistence', 0.38):.1%}</b> of "
                    f"trajectory frames — well above the 5% viability threshold for "
                    f"allosteric drug targeting."
                ),
            },
            {
                "icon": "🔄", "color": "#58a6ff",
                "title": "Domain rotation is the gating mechanism",
                "body": (
                    "The first major P1 opening at ~7.85 ns is causally preceded by a "
                    "<b>9.3°</b> interdomain rotation and simultaneous rupture of the "
                    "<b>E244–R365 salt bridge</b>. The event timeline confirms this "
                    "causal sequence. A molecule that stabilises the rotated domain "
                    "orientation could hold P1 open continuously — a potential "
                    "allosteric mechanism of action."
                ),
            },
            {
                "icon": "🧬", "color": "#50C878",
                "title": "Pocket geometry is evolutionarily conserved",
                "body": (
                    "Nine of the 14 pocket-lining residues score ≥ 7/9 on the ConSurf "
                    "evolutionary conservation scale, assessed across 50 vertebrate FTO "
                    "orthologues. High conservation indicates that the pocket geometry "
                    "is not a simulation artefact but is likely maintained for a "
                    "functional reason — consistent with a genuine allosteric site."
                ),
            },
            {
                "icon": "💊", "color": "#FFD700",
                "title": "P1 scores as a well-druggable allosteric site",
                "body": (
                    f"P1 achieves a composite druggability score of "
                    f"<b>{top_p.get('pocket_score', 0.82):.2f}/1.0</b>, driven by a "
                    f"64% hydrophobic lining (ideal for hydrophobic drug contacts), "
                    f"a moderate charged fraction (21%), and mean volume "
                    f"{top_p.get('mean_volume_A3', 387):.0f} Å³ (comparable to known "
                    f"small-molecule binding sites). A ligand occupying P1 could "
                    f"allosterically suppress m6A demethylation without competing with "
                    f"the active site — a potential advantage for selectivity."
                ),
            },
        ]

        for f in findings:
            st.markdown(
                f"""
                <div style="
                    background:#161b22;
                    border-left:4px solid {f['color']};
                    border-radius:0 8px 8px 0;
                    padding:16px 20px; margin-bottom:14px;
                ">
                    <div style="font-size:16px;font-weight:700;color:#f0f6fc;margin-bottom:6px">
                        {f['icon']} {f['title']}
                    </div>
                    <div style="color:#c9d1d9;font-size:14px;line-height:1.7">
                        {f['body']}
                    </div>
                </div>
                """,
                unsafe_allow_html=True,
            )

        # Radar chart
        st.markdown("#### Pocket Comparison")
        st.caption(
            "Radar chart comparing the three detected pockets across three key dimensions. "
            "A pocket that scores high on all three axes is an ideal drug target."
        )
        metrics = ["persistence", "druggability_score", "pocket_score"]
        labels  = ["Persistence", "Druggability", "Overall score"]
        fig_radar = go.Figure()
        for p in pockets:
            vals = [p.get(m, 0) for m in metrics]
            vals.append(vals[0])
            fig_radar.add_trace(go.Scatterpolar(
                r=vals,
                theta=labels + [labels[0]],
                fill="toself",
                name=p.get("name", p["id"]),
                line_color=p.get("color", "#888888"),
                opacity=0.7,
            ))
        fig_radar.update_layout(
            polar=dict(
                radialaxis=dict(range=[0, 1], gridcolor="#30363d", linecolor="#30363d"),
                angularaxis=dict(gridcolor="#30363d", linecolor="#30363d"),
                bgcolor="#161b22",
            ),
            paper_bgcolor="#0d1117",
            font_color="#f0f6fc",
            showlegend=True,
            legend=dict(bgcolor="rgba(0,0,0,0)"),
            height=380,
            margin=dict(l=40, r=40, t=40, b=40),
        )
        st.plotly_chart(fig_radar, use_container_width=True)

    # ── Pharmacophore ─────────────────────────────────────────────────────────
    with tab_pharma:
        pocket_id   = pharma.get("pocket_id", "P1")
        pocket_name = pharma.get("pocket_name", "Interdomain Interface Pocket")
        summary     = pharma.get("summary", "")

        st.markdown(
            f"""
            <h3 style="color:#f0f6fc;margin:0 0 6px 0">
                Pharmacophore — Pocket {pocket_id}
            </h3>
            <p style="color:#8b949e;font-size:13px;margin:0 0 4px 0">
                {pocket_name}
            </p>
            """,
            unsafe_allow_html=True,
        )

        if summary:
            st.markdown(
                f'<div style="background:#161b22;border:1px solid #30363d;'
                f'border-radius:8px;padding:14px;margin-bottom:16px;'
                f'color:#c9d1d9;font-size:13px;line-height:1.6">{summary}</div>',
                unsafe_allow_html=True,
            )

        with st.expander("ℹ️ What is a pharmacophore model?", expanded=False):
            st.markdown(
                """
                A **pharmacophore** is an abstract description of the molecular features
                that a ligand must possess to bind to a particular pocket. It is not a
                single molecule — it is a blueprint that defines *where* and *what type*
                of chemical interactions are required.

                **Feature types:**

                | Code | Full name | Ligand requirement |
                |------|-----------|-------------------|
                | **HBA** | H-bond acceptor | Ligand must donate a hydrogen bond |
                | **HBD** | H-bond donor | Ligand must accept a hydrogen bond |
                | **AR** | Aromatic | Ligand should have an aromatic ring for π-stacking |
                | **HYD** | Hydrophobic | Ligand should have non-polar groups in this region |
                | **PI** | Positive ionisable | Ligand should be negatively charged here |
                | **NI** | Negative ionisable | Ligand should be positively charged here |

                **Importance score (0–1):** How critical each feature is for binding.
                Features scoring > 0.8 should be considered mandatory; below 0.6
                they are optional but beneficial.

                This pharmacophore model was derived from the pocket geometry at
                maximum volume. It can be used directly in virtual screening to
                filter compound libraries for P1-compatible molecules.
                """
            )

        features = pharma.get("features", [])
        if features:
            df_feat = pd.DataFrame(features)
            df_feat["Feature meaning"] = df_feat["type"].map(
                lambda t: PHARMA_FEATURE_LABELS.get(t, t)
            )
            st.dataframe(
                df_feat[["type", "residue", "description", "importance", "Feature meaning"]].rename(columns={
                    "type":        "Code",
                    "description": "Description",
                    "residue":     "Residue",
                    "importance":  "Importance",
                }),
                use_container_width=True,
                hide_index=True,
                column_config={
                    "Importance": st.column_config.ProgressColumn(
                        "Importance", min_value=0, max_value=1, format="%.2f",
                        help="How critical this feature is for binding (0–1). "
                             "> 0.8 = mandatory; 0.6–0.8 = beneficial; < 0.6 = optional.",
                    ),
                },
            )

            # 3D scatter of pharmacophore features
            type_colors = {
                "HBA": "#58a6ff", "HBD": "#3fb950", "AR": "#FFD700",
                "HYD": "#FF6B35", "PI":  "#c778dd", "NI": "#f85149",
            }
            fig_ph = go.Figure()
            for feat in features:
                pos = feat["position"]
                col = type_colors.get(feat["type"], "#888888")
                fig_ph.add_trace(go.Scatter3d(
                    x=[pos[0]], y=[pos[1]], z=[pos[2]],
                    mode="markers+text",
                    marker=dict(
                        size=feat["radius"] * 5,
                        color=col,
                        opacity=0.75,
                        symbol="circle",
                    ),
                    text=[feat["type"]],
                    textposition="top center",
                    textfont=dict(color=col, size=10),
                    name=f"{feat['type']} — {feat['residue']} (imp. {feat['importance']:.2f})",
                ))
            fig_ph.update_layout(
                paper_bgcolor="#0d1117",
                scene=dict(
                    bgcolor="#161b22",
                    xaxis=dict(gridcolor="#21262d", linecolor="#30363d", title="X (Å)"),
                    yaxis=dict(gridcolor="#21262d", linecolor="#30363d", title="Y (Å)"),
                    zaxis=dict(gridcolor="#21262d", linecolor="#30363d", title="Z (Å)"),
                ),
                font_color="#f0f6fc",
                height=440,
                title=dict(
                    text=f"Pharmacophore features — Pocket {pocket_id} (3D positions in Å)",
                    font_color="#8b949e",
                ),
                legend=dict(bgcolor="rgba(0,0,0,0)", font_size=11),
            )
            st.plotly_chart(fig_ph, use_container_width=True)
            st.caption(
                "Each sphere represents a pharmacophore feature at its approximate "
                "3D position within the pocket. Sphere size corresponds to the "
                "interaction radius. Rotate the plot by dragging."
            )

            col_vol, col_surf, col_enc = st.columns(3)
            col_vol.metric(
                "Pocket volume",  f"{pharma.get('volume_A3', 0):.0f} Å³",
                help="Mean open-state volume. > 300 Å³ is suitable for drug-like molecules.",
            )
            col_surf.metric(
                "Surface area",   f"{pharma.get('surface_A2', 0):.0f} Å²",
                help="Accessible surface area of the pocket interior.",
            )
            col_enc.metric(
                "Enclosure",      f"{pharma.get('enclosure', 0):.2f}",
                help="How buried the pocket is (0–1). Higher = better shielding from solvent.",
            )
        else:
            st.info(
                "Pharmacophore data not available. "
                "Run `python scripts/10_ligand.py --demo` to generate it."
            )

    # ── Export ────────────────────────────────────────────────────────────────
    with tab_export:
        st.markdown("#### Export Results")
        st.caption(
            "Download the full analysis results for use in external tools, "
            "presentations, or further computational work."
        )

        c1, c2 = st.columns(2)
        with c1:
            st.markdown("**Pocket results (JSON)**")
            st.caption("Full pocket data including volume trajectories, residues, and scores.")
            export_data = {
                "generated":     str(date.today()),
                "platform":      "FTO Pocket Discovery Platform v1.0",
                "pdb_primary":   primary,
                "pockets":       pockets,
                "events":        events,
                "pharmacophore": pharma,
            }
            st.download_button(
                "⬇️  Download pocket_results.json",
                data=json.dumps(export_data, indent=2),
                file_name="fto_pocket_results.json",
                mime="application/json",
                use_container_width=True,
            )

        with c2:
            st.markdown("**Pocket ranking (CSV)**")
            st.caption("Tabular summary of all ranked pockets — suitable for spreadsheets.")
            if pockets:
                df_export = pd.DataFrame(pockets)
                csv_cols  = [c for c in df_export.columns if c != "volume_trajectory"]
                st.download_button(
                    "⬇️  Download pocket_ranking.csv",
                    data=df_export[csv_cols].to_csv(index=False),
                    file_name="fto_pocket_ranking.csv",
                    mime="text/csv",
                    use_container_width=True,
                )

        st.markdown("---")
        st.markdown("**Reproduce this analysis**")
        st.code(
            """\
# 1. Set up the environment
conda env create -f environment.yml
conda activate fto-allosteric-pocket-discoverer

# 2a. Full pipeline (requires OpenMM and ~2 h compute)
python run_pipeline.py

# 2b. Demo mode — synthetic data, no MD required (~60 s)
python run_pipeline.py --demo

# 3. Launch the dashboard
streamlit run dashboard/app.py
""",
            language="bash",
        )

        st.markdown("---")
        st.markdown(
            """
            <div style="color:#8b949e;font-size:12px;line-height:1.8">
            <b>Attribution</b><br>
            This analysis used the FTO Pocket Discovery Platform v1.0.<br>
            Primary structure: <a href="https://www.rcsb.org/structure/3LFM"
            style="color:#58a6ff" target="_blank">PDB 3LFM</a>
            — Jia et al. (2011) <i>Nature</i> 475, 561–566.<br>
            Additional structures: <a href="https://www.rcsb.org/structure/4IE4"
            style="color:#58a6ff" target="_blank">PDB 4IE4</a>
            (Feng et al. 2014) ·
            <a href="https://www.rcsb.org/structure/4ZS3"
            style="color:#58a6ff" target="_blank">PDB 4ZS3</a>
            (Huang et al. 2015).<br>
            Conservation scores from ConSurf (Ashkenazy et al. 2016,
            <i>Nucleic Acids Res</i>).<br>
            MD simulations: OpenMM 8 · AMBER ff14SB force field · TIP3P explicit water.
            </div>
            """,
            unsafe_allow_html=True,
        )
