"""
[Events] tab – conformational event timeline, contact network, heatmap.
"""
from __future__ import annotations

import json
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

EVENT_COLORS = {
    "pocket_open":        "#FF6B35",
    "domain_motion":      "#58a6ff",
    "contact_break":      "#f85149",
    "contact_form":       "#3fb950",
    "cluster_transition": "#FFD700",
}

EVENT_ICONS = {
    "pocket_open":        "🕳️",
    "domain_motion":      "🔄",
    "contact_break":      "💥",
    "contact_form":       "🔗",
    "cluster_transition": "🔀",
}

EVENT_DESCRIPTIONS = {
    "pocket_open": (
        "A transient pocket opens to a volume large enough to accommodate a small molecule. "
        "These moments represent windows of opportunity for allosteric ligand binding."
    ),
    "domain_motion": (
        "A significant change in the relative orientation of the catalytic and C-terminal domains. "
        "Domain rotation is the primary driver of P1 pocket opening."
    ),
    "contact_break": (
        "A stable inter-residue interaction (salt bridge or H-bond) that was present in the "
        "closed state breaks apart. This typically gates the opening of the pocket."
    ),
    "contact_form": (
        "A new inter-residue interaction forms that was absent in the closed state. "
        "New contacts can stabilise the open conformation and prolong pocket accessibility."
    ),
    "cluster_transition": (
        "The trajectory crosses from one conformational cluster to another, indicating a "
        "transition between distinct protein states in the energy landscape."
    ),
}


def _plotly_theme() -> dict:
    return dict(
        paper_bgcolor="#0d1117", plot_bgcolor="#161b22",
        font=dict(color="#f0f6fc", size=12),
        xaxis=dict(gridcolor="#21262d", linecolor="#30363d"),
        yaxis=dict(gridcolor="#21262d", linecolor="#30363d"),
        margin=dict(l=60, r=20, t=45, b=50),
    )


def _load_events() -> list[dict]:
    path = RESULTS_DIR / "events.json"
    if path.exists():
        with open(path) as f:
            return json.load(f)
    return DemoData().events()


def _load_network() -> dict:
    path = RESULTS_DIR / "residue_network.json"
    if path.exists():
        with open(path) as f:
            return json.load(f)
    return DemoData().residue_network()


def render() -> None:
    config   = load_config()
    primary  = config["fto"]["primary_structure"]
    events   = _load_events()
    network  = _load_network()
    sim_ns   = config["simulation"]["production_ns"]

    st.markdown(
        """
        <h2 style="color:#f0f6fc;margin:0 0 4px 0">Events</h2>
        <p style="color:#8b949e;font-size:14px;margin:0 0 20px 0">
        Conformational event timeline · contact changes · allosteric mechanism
        </p>
        """,
        unsafe_allow_html=True,
    )

    with st.expander("ℹ️ What are conformational events and why do they matter?", expanded=False):
        st.markdown(
            """
            Proteins are not static — they move continuously. **Conformational events**
            are moments during the simulation where something structurally significant
            happens: a pocket opens, two domains rotate relative to each other, or a
            key salt bridge breaks.

            Understanding the *sequence* of events is crucial for mechanism. For
            example, if a domain rotation consistently precedes pocket opening, it
            suggests that the rotation is the *cause* and the pocket is the *effect*.
            This causal chain can be exploited in drug design: a molecule that locks
            the domain in the rotated state could hold the pocket open.

            **Event types explained:**

            | Type | Icon | What it means |
            |------|------|---------------|
            | Pocket open | 🕳️ | A transient cavity exceeds the volume threshold |
            | Domain motion | 🔄 | Significant change in inter-domain angle |
            | Contact break | 💥 | A stabilising salt bridge or H-bond ruptures |
            | Contact form | 🔗 | A new stabilising interaction appears |
            | Cluster transition | 🔀 | Protein moves between conformational states |

            The **timeline** below shows when each event occurs relative to the total
            simulation length. Vertical proximity on the chart indicates events that
            may be mechanistically linked.
            """
        )

    # ── Event type filter ─────────────────────────────────────────────────────
    all_types = sorted(set(e["type"] for e in events))
    sel_types = st.multiselect(
        "Filter by event type",
        options=all_types,
        default=all_types,
        format_func=lambda t: f"{EVENT_ICONS.get(t, '•')} {t.replace('_', ' ').title()}",
        key="event_type_filter",
    )
    filtered = [e for e in events if e["type"] in sel_types]

    # ── Summary counts ────────────────────────────────────────────────────────
    if events:
        type_counts = {t: sum(1 for e in events if e["type"] == t) for t in all_types}
        cols = st.columns(len(type_counts))
        for col, (etype, count) in zip(cols, type_counts.items()):
            col.metric(
                f"{EVENT_ICONS.get(etype, '')} {etype.replace('_', ' ').title()}",
                f"{count}",
            )

    st.markdown("#### Event Timeline")
    if not filtered:
        st.info("No events match the current filter. Try selecting more event types above.")
    else:
        fig_tl = go.Figure()
        for ev in sorted(filtered, key=lambda e: e["time_ps"]):
            t_start = ev["time_ps"] / 1_000
            dur_ns  = ev.get("duration_ps", 200) / 1_000
            t_end   = t_start + dur_ns
            color   = EVENT_COLORS.get(ev["type"], "#888888")
            y_label = f"{ev.get('id', '?')}  {EVENT_ICONS.get(ev['type'], '')} {ev['type'].replace('_',' ')}"

            fig_tl.add_trace(go.Scatter(
                x=[t_start, t_end],
                y=[y_label, y_label],
                mode="lines",
                line=dict(color=color, width=18),
                name=ev.get("id", ""),
                hovertemplate=(
                    f"<b>{ev.get('id', '')} — {ev['type'].replace('_',' ').title()}</b><br>"
                    f"Time: {t_start:.2f}–{t_end:.2f} ns<br>"
                    f"Duration: {dur_ns:.2f} ns<br>"
                    f"{ev['description']}"
                    "<extra></extra>"
                ),
            ))
            fig_tl.add_trace(go.Scatter(
                x=[t_start], y=[y_label],
                mode="markers",
                marker=dict(color=color, size=10, symbol="circle"),
                showlegend=False,
                hoverinfo="skip",
            ))

        _theme = {k: v for k, v in _plotly_theme().items() if k not in ("xaxis", "yaxis")}
        fig_tl.update_layout(
            **_theme,
            height=max(260, 50 * len(filtered) + 80),
            xaxis=dict(
                title="Simulation time (ns)",
                range=[0, sim_ns],
                gridcolor="#21262d",
            ),
            yaxis=dict(autorange="reversed"),
            showlegend=False,
        )
        st.plotly_chart(fig_tl, use_container_width=True)
        st.caption(
            "Each bar spans the duration of one event. Hover for details. "
            "The sequence of events — rotation → contact break → pocket open — "
            "suggests a mechanistic cascade, not random fluctuation."
        )

    # ── Event detail card ─────────────────────────────────────────────────────
    st.markdown("#### Event Details")

    ev_ids = [e.get("id", f"E{i}") for i, e in enumerate(filtered)]
    if not ev_ids:
        st.info("No events to display. Adjust the filter above.")
        return

    sel_ev_id = st.selectbox("Select event", ev_ids)
    sel_ev    = next((e for e in filtered if e.get("id") == sel_ev_id), filtered[0])
    ev_color  = EVENT_COLORS.get(sel_ev["type"], "#888888")
    ev_icon   = EVENT_ICONS.get(sel_ev["type"], "📌")
    ev_explain = EVENT_DESCRIPTIONS.get(sel_ev["type"], "")

    metrics    = sel_ev.get("metrics", {})
    metric_html = "".join([
        f'<span style="display:inline-block;background:#0d1117;border:1px solid #30363d;'
        f'border-radius:6px;padding:8px 14px;margin:4px">'
        f'<span style="color:#8b949e;font-size:11px">{k.replace("_", " ").title()}</span><br>'
        f'<span style="color:#58a6ff;font-weight:700;font-size:18px">{v}</span></span>'
        for k, v in metrics.items()
    ])

    st.markdown(
        f"""
        <div style="
            background:#161b22;border:1px solid {ev_color};
            border-radius:10px;padding:20px;margin:12px 0;
        ">
            <div style="font-size:11px;text-transform:uppercase;letter-spacing:.06em;
                        color:{ev_color}">{ev_icon} {sel_ev['type'].replace('_', ' ')}</div>
            <div style="font-size:18px;font-weight:700;color:#f0f6fc;margin:6px 0">
                {sel_ev.get('id', '')} &nbsp;·&nbsp; {sel_ev['description']}
            </div>
            <div style="color:#8b949e;font-size:13px;margin-bottom:12px">
                Time: <b>{sel_ev['time_ps'] / 1000:.2f} ns</b>
                &nbsp;·&nbsp; Duration: <b>{sel_ev.get('duration_ps', 0) / 1000:.2f} ns</b>
                &nbsp;·&nbsp; Frame: <b>{sel_ev.get('frame', '?')}</b>
            </div>
            <div style="color:#c9d1d9;font-size:13px;margin-bottom:12px;
                        background:#0d1117;border-radius:6px;padding:10px">
                {ev_explain}
            </div>
            <div>{metric_html}</div>
        </div>
        """,
        unsafe_allow_html=True,
    )

    # ── Residue contact network ───────────────────────────────────────────────
    st.markdown("#### Residue Contact Network — Changes at Pocket Opening")
    st.caption(
        "These tables show inter-residue contacts that change when the P1 pocket opens. "
        "Breaking contacts (red) lose their interaction; forming contacts (green) are new. "
        "Δ probability indicates the change in contact frequency between the open and closed states."
    )

    breaking = network.get("breaking_contacts", [])
    forming  = network.get("forming_contacts", [])
    c_break, c_form = st.columns(2)

    def _contact_table(contacts: list[dict], title: str, color: str) -> None:
        st.markdown(
            f'<div style="color:{color};font-weight:700;font-size:14px;margin-bottom:8px">'
            f'{title} ({len(contacts)})</div>',
            unsafe_allow_html=True,
        )
        if not contacts:
            st.caption("None detected")
            return
        df = pd.DataFrame(contacts)
        display_cols = ["res1", "res2", "label"]
        if "delta" in df.columns:
            df["Δ probability"] = df["delta"].apply(lambda x: f"{x:+.2f}")
            display_cols.append("Δ probability")
        rename = {"res1": "Residue 1", "res2": "Residue 2", "label": "Interaction"}
        st.dataframe(
            df[display_cols].rename(columns=rename),
            use_container_width=True, hide_index=True,
        )

    with c_break:
        _contact_table(breaking, "💥 Breaking contacts", "#f85149")
    with c_form:
        _contact_table(forming, "🔗 Forming contacts", "#3fb950")

    # ── Contact delta heatmap ─────────────────────────────────────────────────
    with st.expander("Contact change map (heatmap)"):
        st.caption(
            "Δ contact probability for residue pairs between open and closed states. "
            "Red = lost contacts (destabilised), green = gained contacts (stabilised)."
        )
        all_contacts = breaking + forming
        if all_contacts:
            all_res = sorted(set(r for c in all_contacts for r in [c["res1"], c["res2"]]))
            n       = len(all_res)
            mat     = np.zeros((n, n))
            res_idx = {r: i for i, r in enumerate(all_res)}
            for c in all_contacts:
                i = res_idx[c["res1"]]
                j = res_idx[c["res2"]]
                mat[i, j] = c["delta"]
                mat[j, i] = c["delta"]

            fig_hm = go.Figure(go.Heatmap(
                z=mat,
                x=[str(r) for r in all_res],
                y=[str(r) for r in all_res],
                colorscale=[
                    [0,   "#f85149"],
                    [0.5, "#161b22"],
                    [1.0, "#3fb950"],
                ],
                zmid=0, zmin=-1, zmax=1,
                colorbar=dict(title="Δ contact<br>probability", thickness=12),
            ))
            fig_hm.update_layout(
                **_plotly_theme(), height=420,
                xaxis_title="Residue",
                yaxis_title="Residue",
                title=dict(
                    text="Δ Contact Probability: Open State vs Closed State",
                    font_color="#8b949e",
                ),
            )
            st.plotly_chart(fig_hm, use_container_width=True)
        else:
            st.caption("No contact data available.")
