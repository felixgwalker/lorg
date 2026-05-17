"""
[Structure] tab – curated FTO structures, 3D viewer, metadata.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import requests
import streamlit as st

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from scripts.utils import DATA_DIR, DemoData, load_config
from dashboard.components.viewer import render_structure

STRUCTURES_DIR = DATA_DIR / "structures"


@st.cache_data(ttl=3600, show_spinner=False)
def _fetch_pdb(pdb_id: str) -> str | None:
    path = STRUCTURES_DIR / f"{pdb_id}.pdb"
    if path.exists():
        return path.read_text(encoding="utf-8", errors="ignore")
    try:
        r = requests.get(f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=20)
        if r.ok:
            return r.text
    except Exception:
        pass
    return None


def _load_metadata() -> list[dict]:
    meta_path = STRUCTURES_DIR / "metadata.json"
    if meta_path.exists():
        with open(meta_path) as f:
            d = json.load(f)
        return list(d.values())
    return DemoData().structure_metadata()


def _get_resolution(meta: dict) -> str:
    """Return resolution string, handling both 'resolution' and 'resolution_A' keys."""
    val = meta.get("resolution") or meta.get("resolution_A")
    if val is None or val == "?":
        return "?"
    try:
        return f"{float(val):.2f}"
    except (TypeError, ValueError):
        return str(val)


def render() -> None:
    config = load_config()

    st.markdown(
        """
        <div style="margin-bottom:20px">
            <h2 style="color:#f0f6fc;margin:0">Structure Set</h2>
            <p style="color:#8b949e;margin:4px 0 0 0;font-size:14px">
                Curated human FTO crystal structures from the
                <a href="https://www.rcsb.org" target="_blank"
                   style="color:#58a6ff">Protein Data Bank (RCSB)</a>
            </p>
        </div>
        """,
        unsafe_allow_html=True,
    )

    # ── Scientific context ────────────────────────────────────────────────────
    with st.expander("ℹ️ About these structures", expanded=False):
        st.markdown(
            """
            **What are these structures?**
            Each entry is an experimentally determined three-dimensional structure of
            the human FTO protein, solved by X-ray crystallography and deposited in
            the public Protein Data Bank (PDB).

            **Why multiple structures?**
            Different crystal structures capture FTO in different states: bound to
            substrates, inhibitors, or cofactors. Comparing them reveals how the
            protein moves and what sites are available for drug binding.

            **Resolution (Å):** Lower numbers mean higher precision.
            Structures below 2.5 Å are considered high-resolution and are suitable
            for detailed analysis of binding geometry.

            **Ligands shown:** Small molecules co-crystallised with the protein.
            *3DT* = 3-methylthymidine (substrate analogue), *FE2* = iron(II),
            *OGA* = N-oxalylglycine (α-KG analogue), *8XQ* / *A4F* = inhibitors.
            """
        )

    metadata = _load_metadata()

    # ── Metadata table ────────────────────────────────────────────────────────
    rows = []
    for m in metadata:
        rows.append({
            "PDB ID":      m.get("pdb_id", ""),
            "Title":       m.get("title", ""),
            "Resolution (Å)": _get_resolution(m),
            "Method":      m.get("method", ""),
            "Chains":      ", ".join(m.get("chains", [])) if m.get("chains") else "?",
            "Residues":    m.get("n_residues", "?"),
            "Ligands":     ", ".join(m.get("ligands", [])) if m.get("ligands") else "none",
        })
    df_display = pd.DataFrame(rows)

    st.markdown("#### Curated FTO Structure Set")
    st.dataframe(df_display, use_container_width=True, hide_index=True)

    # ── Structure selector ────────────────────────────────────────────────────
    pdb_ids = [m["pdb_id"] for m in metadata]
    primary = config["fto"]["primary_structure"]
    default = pdb_ids.index(primary) if primary in pdb_ids else 0

    col1, col2, col3 = st.columns([2, 2, 3])
    with col1:
        selected = st.selectbox("Select structure", pdb_ids, index=default)
    with col2:
        view_style = st.selectbox("Style", ["cartoon", "surface", "stick", "sphere"], index=0)
    with col3:
        spin = st.toggle("Auto-rotate", value=False)

    # ── Load and show structure ───────────────────────────────────────────────
    with st.spinner(f"Loading {selected} …"):
        pdb_content = _fetch_pdb(selected)

    if pdb_content is None:
        st.warning(
            f"Could not load {selected}. "
            "Run `python scripts/01_acquire.py` to download structures, "
            "or check your internet connection."
        )
        return

    cat   = config["fto"]["domains"]["catalytic"]
    cterm = config["fto"]["domains"]["cterminal"]
    active_site = config["fto"]["active_site_residues"]

    # Selected structure info card
    sel_meta = next((m for m in metadata if m["pdb_id"] == selected), {})
    res  = _get_resolution(sel_meta)
    meth = sel_meta.get("method", "?")
    ligs = sel_meta.get("ligands", [])
    ch   = sel_meta.get("chains", [])

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Resolution", f"{res} Å" if res != "?" else "?")
    m2.metric("Method", meth)
    m3.metric("Chains", ", ".join(ch) if ch else "?")
    m4.metric("Ligands", ", ".join(ligs) if ligs else "none")

    render_structure(
        pdb_content,
        height=540,
        spin=spin,
        style=view_style,
        highlight_residues=active_site,
        pocket_residues=[],
        cat_range=(cat["start"], cat["end"]),
        cterm_range=(cterm["start"], cterm["end"]),
        label=f"{selected} · {res} Å · {meth}",
    )

    st.caption(
        "**Blue** = catalytic domain (AlkB-like fold) · "
        "**Green** = C-terminal domain · "
        "**Yellow spheres** = active-site residues"
    )

    # ── Domain info ───────────────────────────────────────────────────────────
    st.markdown("#### Domain Architecture")
    c1, c2 = st.columns(2)
    with c1:
        st.markdown(
            f"""
            <div style="
                background:#161b22;border:1px solid #4A90D9;
                border-radius:8px;padding:16px;
            ">
                <div style="font-size:11px;text-transform:uppercase;
                            letter-spacing:.06em;color:#4A90D9">Catalytic Domain</div>
                <div style="font-size:20px;font-weight:700;color:#f0f6fc;margin:6px 0">
                    Residues {cat['start']}–{cat['end']}
                </div>
                <div style="color:#8b949e;font-size:13px">
                    AlkB-like double-stranded β-helix (DSBH) fold.<br>
                    Houses Fe²⁺ ion and α-ketoglutarate cofactor.<br>
                    Iron-binding triad: <b>H231 · D233 · H307</b><br><br>
                    This domain performs the demethylation chemistry.
                    It is the site of known inhibitors and the primary
                    drug-discovery target.
                </div>
            </div>
            """,
            unsafe_allow_html=True,
        )
    with c2:
        st.markdown(
            f"""
            <div style="
                background:#161b22;border:1px solid #50C878;
                border-radius:8px;padding:16px;
            ">
                <div style="font-size:11px;text-transform:uppercase;
                            letter-spacing:.06em;color:#50C878">C-terminal Domain</div>
                <div style="font-size:20px;font-weight:700;color:#f0f6fc;margin:6px 0">
                    Residues {cterm['start']}–{cterm['end']}
                </div>
                <div style="color:#8b949e;font-size:13px">
                    Unique to FTO; absent in other AlkB-family members.<br>
                    Packs tightly against the catalytic domain.<br>
                    Contains the proposed allosteric interface.<br><br>
                    The interface between these two domains is where
                    the transient <b>P1 pocket</b> forms during simulation.
                </div>
            </div>
            """,
            unsafe_allow_html=True,
        )

    # ── Active site residues ──────────────────────────────────────────────────
    with st.expander("Active site and key residues"):
        st.markdown(
            """
            The table below lists residues critical to FTO catalytic function.
            **Iron-binding residues** coordinate the Fe²⁺ cofactor that is
            essential for oxidative demethylation. **Substrate-binding residues**
            contact the N6-methyladenosine (m6A) nucleotide during catalysis.
            These residues are highlighted yellow in the 3D viewer.
            """
        )
        iron  = config["fto"]["iron_binding_residues"]
        subs  = config["fto"]["substrate_binding_residues"]
        rows  = []
        for r in sorted(set(iron + subs + active_site)):
            role = []
            if r in iron:  role.append("Fe²⁺ binding")
            if r in subs:  role.append("Substrate recognition")
            rows.append({"Residue": r, "Role": ", ".join(role) or "Active site"})
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
