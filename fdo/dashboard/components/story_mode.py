"""
Story Mode – guided narrative walkthrough of FTO findings.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import streamlit as st

sys.path.insert(0, str(Path(__file__).parent.parent.parent))


STEPS = [
    {
        "title": "Step 1 · Starting conformation",
        "icon": "🏗️",
        "text": (
            "Human FTO (fat mass and obesity-associated protein) is a 498-residue "
            "AlkB-family RNA demethylase that removes the N6-methyladenosine (m6A) "
            "modification from messenger RNA — a key epigenetic mark that regulates "
            "gene expression. The crystal structure (<b>PDB: 3LFM</b>, 2.5 Å) reveals "
            "two tightly packed domains: a <b>catalytic AlkB-like domain</b> (residues "
            "31–326, blue) housing the Fe²⁺/α-ketoglutarate active site, and a "
            "<b>C-terminal domain</b> (residues 327–498, green) unique to FTO within "
            "the AlkB family."
        ),
        "metric_label": None,
        "focus": "full_structure",
    },
    {
        "title": "Step 2 · Domain interface",
        "icon": "🔗",
        "text": (
            "The two domains share an approximately <b>1,200 Å²</b> interface stabilised "
            "by a key salt bridge <b>E244–R365</b>, two hydrogen bonds "
            "(<b>Q108–N338</b>, <b>Y205–T332</b>), and several hydrophobic contacts. "
            "This interface is the proposed allosteric hotspot: disrupting it could "
            "modulate FTO catalytic activity without directly competing with the "
            "highly conserved active site — a significant pharmacological advantage, "
            "since active-site inhibitors must outcompete natural m6A substrates."
        ),
        "metric_label": "Interface area",
        "metric_value": "~1,200 Å²",
        "focus": "interface",
    },
    {
        "title": "Step 3 · Interdomain motion",
        "icon": "🔄",
        "text": (
            "Across 20 ns of unbiased NPT molecular dynamics simulation, the catalytic "
            "domain undergoes <b>transient rotations</b> of up to <b>12°</b> relative to "
            "the C-terminal domain. The first major event occurs at <b>~7.5 ns</b>, where "
            "the two domains separate by an additional <b>3.2 Å</b> and the opening angle "
            "increases by <b>9.3°</b>. This motion is concomitant with rupture of the "
            "E244–R365 salt bridge. Because the simulation is unbiased, this represents "
            "a thermally accessible motion — not an artificially forced event."
        ),
        "metric_label": "Peak rotation",
        "metric_value": "12°",
        "focus": "motion",
    },
    {
        "title": "Step 4 · Transient pocket P1",
        "icon": "🕳️",
        "text": (
            "Within <b>200 ps</b> of the domain rotation, a new cavity — <b>Pocket P1</b> — "
            "emerges at the catalytic–CTD interface. P1 reaches a peak volume of "
            "<b>487 Å³</b> and persists in <b>38% of all trajectory frames</b>. "
            "It re-opens in two further events (at ~12 ns and ~16 ns), confirming it "
            "is a genuine metastable state rather than a random thermal fluctuation. "
            "For context, known drug-binding pockets typically range from 300–1000 Å³."
        ),
        "metric_label": "P1 persistence",
        "metric_value": "38%",
        "focus": "pocket",
    },
    {
        "title": "Step 5 · Key residues",
        "icon": "🔬",
        "text": (
            "Pocket P1 is lined by residues from both domains: "
            "<b>Y106, L109, F130, Y205</b> (catalytic domain) and "
            "<b>L335, T332, D340, R362</b> (CTD). "
            "Several are highly conserved (ConSurf ≥ 7/9) across 50 vertebrate FTO "
            "orthologues, indicating the pocket geometry is evolutionarily maintained — "
            "a hallmark of functional relevance. The adjacent active-site residues "
            "<b>H231, D233, H307</b> are invariant across all species examined (score 9/9) "
            "and show low RMSF, consistent with their essential catalytic role."
        ),
        "metric_label": "Highly conserved lining residues",
        "metric_value": "9 / 14",
        "focus": "residues",
    },
    {
        "title": "Step 6 · Allosteric potential",
        "icon": "💊",
        "text": (
            "Pocket P1 scores <b>0.82 / 1.0</b> on the composite druggability index "
            "(persistence 38%, mean volume 387 Å³, druggability 0.74). "
            "The pharmacophore model identifies <b>5 key binding features</b>: an "
            "H-bond acceptor (His231), an H-bond donor (Asp233), an aromatic centre "
            "(Tyr106), a hydrophobic patch (Leu109), and a positive-ionisable region "
            "(Arg96). A small molecule complementary to these features could "
            "<b>allosterically modulate m6A demethylation</b> without directly "
            "competing with the RNA substrate — a potential selectivity advantage "
            "over active-site inhibitors."
        ),
        "metric_label": "Composite druggability score",
        "metric_value": "0.82 / 1.0",
        "focus": "pharmacophore",
    },
]


def story_mode_ui(
    pockets: list[dict] | None = None,
    pdb_content: str | None = None,
) -> None:
    """Render the full story mode walkthrough."""
    st.markdown(
        """
        <div style="
            background: linear-gradient(135deg, #161b22 0%, #0d1117 100%);
            border: 1px solid #30363d; border-radius: 10px;
            padding: 20px 24px; margin-bottom: 20px;
        ">
            <h2 style="color:#58a6ff;margin:0 0 6px 0">📖 Story Mode</h2>
            <p style="color:#8b949e;margin:0;font-size:14px">
                A guided, plain-language walkthrough of the key FTO pocket-discovery
                findings. Use the arrows below to step through the narrative, or
                enable Autoplay to advance automatically.
            </p>
        </div>
        """,
        unsafe_allow_html=True,
    )

    # ── Session state init ────────────────────────────────────────────────────
    if "story_step" not in st.session_state:
        st.session_state.story_step = 0
    if "story_autoplay" not in st.session_state:
        st.session_state.story_autoplay = False
    if "story_last_advance" not in st.session_state:
        st.session_state.story_last_advance = 0.0

    n       = len(STEPS)
    step_idx = st.session_state.story_step

    # ── Navigation ────────────────────────────────────────────────────────────
    nav_l, nav_mid, nav_r = st.columns([1, 6, 1])
    with nav_l:
        if st.button("◀", key="story_prev", disabled=step_idx == 0):
            st.session_state.story_step -= 1
            st.session_state.story_autoplay = False
            st.rerun()
    with nav_mid:
        progress_pct = int((step_idx + 1) / n * 100)
        st.markdown(
            f"""
            <div style="text-align:center;color:#8b949e;font-size:13px;margin-bottom:4px">
                Step {step_idx + 1} of {n}
            </div>
            <div style="background:#161b22;border-radius:4px;height:4px">
                <div style="background:#58a6ff;border-radius:4px;height:4px;
                            width:{progress_pct}%;transition:width 0.3s ease"></div>
            </div>
            """,
            unsafe_allow_html=True,
        )
    with nav_r:
        if st.button("▶", key="story_next", disabled=step_idx == n - 1):
            st.session_state.story_step += 1
            st.session_state.story_autoplay = False
            st.rerun()

    # ── Step card ─────────────────────────────────────────────────────────────
    step = STEPS[step_idx]
    metric_html = ""
    if step.get("metric_label"):
        metric_html = f"""
        <div style="
            display:inline-block; background:#161b22;
            border:1px solid #30363d; border-radius:8px;
            padding:10px 18px; margin-top:14px;
        ">
            <div style="font-size:11px;color:#8b949e;text-transform:uppercase;
                        letter-spacing:0.06em">{step['metric_label']}</div>
            <div style="font-size:28px;font-weight:700;color:#58a6ff;margin-top:2px">
                {step.get('metric_value', '')}
            </div>
        </div>
        """

    st.markdown(
        f"""
        <div style="
            background:#161b22; border:1px solid #30363d;
            border-radius:12px; padding:28px; margin-top:16px;
            min-height:200px;
        ">
            <h3 style="color:#f0f6fc;margin:0 0 12px 0">
                {step['icon']} {step['title']}
            </h3>
            <p style="color:#c9d1d9;line-height:1.75;font-size:15px;margin:0">
                {step['text']}
            </p>
            {metric_html}
        </div>
        """,
        unsafe_allow_html=True,
    )

    # ── Progress dots ─────────────────────────────────────────────────────────
    dots = "".join(
        f'<div style="width:8px;height:8px;border-radius:50%;'
        f'background:{"#58a6ff" if i == step_idx else "#30363d"}"></div>'
        for i in range(n)
    )
    st.markdown(
        f'<div style="display:flex;gap:8px;justify-content:center;margin-top:14px">'
        f'{dots}</div>',
        unsafe_allow_html=True,
    )

    # ── Autoplay ──────────────────────────────────────────────────────────────
    st.markdown("<div style='height:12px'></div>", unsafe_allow_html=True)
    autoplay_on = st.toggle(
        "▶ Autoplay (4 s per step)",
        key="story_autoplay_toggle",
        value=st.session_state.story_autoplay,
    )
    st.session_state.story_autoplay = autoplay_on

    if autoplay_on:
        now = time.time()
        last = st.session_state.story_last_advance
        delay = 4.0

        remaining = delay - (now - last)
        if remaining > 0:
            time.sleep(remaining)

        if step_idx < n - 1:
            st.session_state.story_step += 1
        else:
            # Reached end — stop autoplay and reset
            st.session_state.story_step = 0
            st.session_state.story_autoplay = False

        st.session_state.story_last_advance = time.time()
        st.rerun()
