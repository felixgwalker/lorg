"""
HTML and JSON report generation.

Produces a single self-contained HTML file with inline CSS — no external
assets, no JavaScript, suitable for archiving and regulatory submission.
"""

import json
from pathlib import Path

from src.models import (
    FeatureResult,
    LayerResult,
    PipelineResult,
    RiskLevel,
    ScoreBand,
    SignalResult,
    ThreeLayerResult,
)

# ---------------------------------------------------------------------------
# Colour scheme
# ---------------------------------------------------------------------------

_RISK_COLOUR = {
    RiskLevel.LOW:      "#2d7d2d",
    RiskLevel.MEDIUM:   "#8a6000",
    RiskLevel.HIGH:     "#b84000",
    RiskLevel.CRITICAL: "#a00000",
}
_RISK_BG = {
    RiskLevel.LOW:      "#e8f5e8",
    RiskLevel.MEDIUM:   "#fff8e0",
    RiskLevel.HIGH:     "#fff0e8",
    RiskLevel.CRITICAL: "#fce8e8",
}

_SIGNAL_LABELS = {
    "is_proximity":  "IS Element Proximity",
    "gc_content":    "GC Content Deviation",
    "integron":      "Integron Association",
    "conjugative":   "Conjugative Element Homology",
    "prophage":      "Prophage Context",
}

_CSS = """
* { box-sizing: border-box; margin: 0; padding: 0; }
body {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
  font-size: 13px; color: #1a1a1a; background: #f4f5f7;
}
.page {
  max-width: 980px; margin: 0 auto; background: #fff;
  box-shadow: 0 2px 24px rgba(0,0,0,.10);
}
header {
  background: #1c2b4a; color: #fff; padding: 28px 40px;
}
header h1 { font-size: 20px; font-weight: 600; letter-spacing: .3px; }
header .meta { font-size: 11px; opacity: .65; margin-top: 6px; }
main { padding: 28px 40px; }
footer {
  background: #f0f2f5; padding: 14px 40px;
  font-size: 11px; color: #888; border-top: 1px solid #dde0e6;
}
section { margin-bottom: 32px; }
h2 {
  font-size: 14px; font-weight: 700; color: #1c2b4a;
  border-bottom: 2px solid #e0e4ea; padding-bottom: 6px; margin-bottom: 14px;
}
h3 { font-size: 13px; font-weight: 600; color: #2a2a2a; margin-bottom: 8px; }
h4 { font-size: 12px; font-weight: 600; color: #444; margin: 10px 0 5px; }
/* Risk badge */
.badge {
  display: inline-block; padding: 7px 18px; border-radius: 4px;
  font-size: 16px; font-weight: 700; letter-spacing: .4px;
}
.idx {
  font-size: 26px; font-weight: 700;
  display: inline-block; margin-left: 18px; vertical-align: middle;
}
/* Summary cards */
.cards { display: flex; gap: 16px; flex-wrap: wrap; margin-top: 16px; }
.card {
  background: #f7f8fa; border: 1px solid #e2e5ea;
  border-radius: 6px; padding: 14px 18px; flex: 1; min-width: 150px;
}
.card .lbl { font-size: 10px; text-transform: uppercase; letter-spacing: .5px;
             color: #888; margin-bottom: 4px; }
.card .val { font-size: 15px; font-weight: 600; }
/* Alerts / warnings */
.alert {
  background: #fffbe6; border-left: 4px solid #f0c040;
  padding: 10px 14px; border-radius: 0 4px 4px 0; margin-bottom: 16px;
}
.warn-inline {
  background: #fff4ec; border-left: 3px solid #e07020;
  padding: 7px 10px; border-radius: 0 3px 3px 0; font-size: 12px;
  color: #6b3000; margin: 6px 0 10px;
}
.note { font-size: 11px; color: #777; font-style: italic; margin: 4px 0 8px; }
/* Tables */
table { width: 100%; border-collapse: collapse; font-size: 12px; }
th {
  background: #f0f2f6; text-align: left; padding: 7px 9px;
  font-weight: 600; border-bottom: 2px solid #d8dce4;
}
td { padding: 6px 9px; border-bottom: 1px solid #ececec; vertical-align: top; }
tr:last-child td { border-bottom: none; }
tr:hover td { background: #fafbfc; }
.hc { color: #b84000; font-weight: 700; }
/* Bar chart */
.bar-wrap {
  display: inline-block; vertical-align: middle;
  background: #e2e5ea; border-radius: 3px; height: 14px; width: 160px;
}
.bar-fill { height: 100%; border-radius: 3px; }
.bar-val { display: inline-block; vertical-align: middle; margin-left: 7px;
           font-weight: 600; }
/* Signal detail blocks */
.sig-block {
  border: 1px solid #e2e5ea; border-radius: 6px;
  padding: 16px 18px; margin-bottom: 16px;
}
code {
  font-family: 'Courier New', monospace; font-size: 11px;
  background: #f4f4f4; padding: 1px 4px; border-radius: 2px;
}
ul { margin-left: 18px; line-height: 1.9; }
.total-row td { font-weight: 700; background: #f0f2f6; }
"""


# ---------------------------------------------------------------------------
# HTML fragments
# ---------------------------------------------------------------------------

def _bar(score: float, colour: str) -> str:
    pct = int(score * 100)
    return (
        f'<span class="bar-wrap">'
        f'<span class="bar-fill" style="width:{pct}%;background:{colour};"></span>'
        f'</span>'
        f'<span class="bar-val">{score:.3f}</span>'
    )


def _skipped_cell() -> str:
    return '<span style="color:#aaa;font-style:italic;">Skipped</span>'


def _signal_rows(signal_results: list[SignalResult],
                 colour: str,
                 active_weight_sum: float) -> str:
    rows = ""
    for s in signal_results:
        label = _SIGNAL_LABELS.get(s.signal_name, s.signal_name)
        if s.skipped or s.score is None:
            score_cell = _skipped_cell()
            contrib = "—"
            dim = ' style="opacity:.45;"'
        else:
            score_cell = _bar(s.score, colour)
            c = s.score * (s.weight / active_weight_sum) if active_weight_sum else 0.0
            contrib = f"{c:.3f}"
            dim = ""
        rows += (
            f"<tr{dim}>"
            f"<td>{label}</td>"
            f"<td>{score_cell}</td>"
            f"<td>{s.weight * 100:.0f}%</td>"
            f"<td>{contrib}</td>"
            f"</tr>\n"
        )
    return rows


# ---------------------------------------------------------------------------
# Per-signal detail sections
# ---------------------------------------------------------------------------

def _detail_gc(s: SignalResult) -> str:
    ev = s.evidence
    return f"""
<div class="sig-block">
  <h3>{_SIGNAL_LABELS['gc_content']}</h3>
  <table style="width:auto;">
    <tr><th>Metric</th><th>Value</th></tr>
    <tr><td>Query GC%</td><td>{ev.get('query_gc', '?')}%</td></tr>
    <tr><td>Host GC%</td><td>{ev.get('host_gc', '?')}%</td></tr>
    <tr><td>Absolute deviation</td><td>{ev.get('deviation_pct', '?')}%</td></tr>
    <tr><td>Score = 1.0 at deviation ≥</td><td>{ev.get('max_deviation_pct', 25)}%</td></tr>
  </table>
  {('<div class="warn-inline">' + s.warning + '</div>') if s.warning else ''}
</div>"""


def _hit_table(hits: list[dict], has_hc: bool = False) -> str:
    if not hits:
        cols = 6 if has_hc else 5
        return f'<tr><td colspan="{cols}" style="color:#aaa;">No hits found</td></tr>'
    rows = ""
    for h in hits:
        hc_cell = (
            f'<td class="hc">&#x26A0;</td>'
            if has_hc and h.get("high_concern") else
            ("<td></td>" if has_hc else "")
        )
        rows += (
            f"<tr>"
            f"<td><code>{h['subject']}</code></td>"
            f"<td>{h.get('description', '')}</td>"
            f"<td>{h['pct_identity']:.1f}%</td>"
            f"<td>{h['query_coverage']:.1f}%</td>"
            f"<td>{h['evalue']:.2e}</td>"
            f"{hc_cell}"
            f"</tr>\n"
        )
    return rows


def _detail_blast_simple(s: SignalResult, label: str) -> str:
    ev = s.evidence
    hits = ev.get("top_hits", [])
    header = (
        "<tr><th>Subject</th><th>Description</th>"
        "<th>Identity</th><th>Coverage</th><th>E-value</th></tr>"
    )
    return f"""
<div class="sig-block">
  <h3>{label}</h3>
  <p>{ev.get('hit_count', 0)} total {label.lower()} hits (top 5 shown).</p>
  {('<div class="warn-inline">' + s.warning + '</div>') if s.warning else ''}
  <table>{header}{_hit_table(hits)}</table>
</div>"""


def _detail_conjugative(s: SignalResult) -> str:
    ev = s.evidence
    hits = ev.get("top_hits", [])
    hc_count = ev.get("high_concern_hit_count", 0)
    markers = ", ".join(ev.get("high_concern_markers_found", [])) or "none"

    header = (
        "<tr><th>Subject</th><th>Description</th>"
        "<th>Identity</th><th>Coverage</th><th>E-value</th>"
        "<th title='High-concern marker'>&#x26A0;</th></tr>"
    )
    return f"""
<div class="sig-block">
  <h3>{_SIGNAL_LABELS['conjugative']}</h3>
  <p>{ev.get('hit_count', 0)} total hits.
     <strong>{hc_count}</strong> high-concern hits
     (relaxase / MOB / T4SS markers). Markers found: {markers}.</p>
  {('<div class="warn-inline">' + s.warning + '</div>') if s.warning else ''}
  <table>{header}{_hit_table(hits, has_hc=True)}</table>
</div>"""


def _detail_integron(s: SignalResult) -> str:
    ev = s.evidence
    attc = ev.get("attc_sites", [])
    blast = ev.get("top_blast_hits", [])

    attc_rows = ""
    for a in attc:
        attc_rows += (
            f"<tr>"
            f"<td>{a['start']}–{a['end']}</td>"
            f"<td>{a['strand']}</td>"
            f"<td>{a['spacer_length']} bp</td>"
            f"<td><code>{a['sequence_prefix']}</code></td>"
            f"</tr>\n"
        )
    if not attc_rows:
        attc_rows = '<tr><td colspan="4" style="color:#aaa;">None found</td></tr>'

    blast_header = (
        "<tr><th>Subject</th><th>Description</th>"
        "<th>Identity</th><th>Coverage</th><th>E-value</th></tr>"
    )

    return f"""
<div class="sig-block">
  <h3>{_SIGNAL_LABELS['integron']}</h3>
  {('<div class="warn-inline">' + s.warning + '</div>') if s.warning else ''}
  <h4>Putative attC sites ({ev.get('attc_sites_found', 0)} found)</h4>
  <p class="note">
    Detected by regex scan of both strands using the attC bottom-strand consensus
    GTT[Y][R][R]…[R][R]AAC (spacer 50–200 bp).  All hits are putative;
    INTEGRALL BLAST confirmation is the authoritative check.
  </p>
  <table>
    <tr><th>Position</th><th>Strand</th><th>Spacer</th><th>Sequence (prefix)</th></tr>
    {attc_rows}
  </table>
  <h4>INTEGRALL BLAST hits ({ev.get('integrall_blast_hits', 0)} total, top 5 shown)</h4>
  <table>{blast_header}{_hit_table(blast)}</table>
</div>"""


def _detail_prophage(s: SignalResult) -> str:
    ev = s.evidence
    source = ev.get("source", "unknown")

    if source == "phaster_api":
        regions = ev.get("regions", [])
        reg_rows = ""
        for r in regions:
            reg_rows += (
                f"<tr>"
                f"<td>{r['id']}</td>"
                f"<td><strong>{r['completeness'].title()}</strong></td>"
                f"<td>{r['start']:,}–{r['end']:,}</td>"
                f"<td>{r['gc_pct']}%</td>"
                f"<td>{r['num_cds']}</td>"
                f"</tr>\n"
            )
        if not reg_rows:
            reg_rows = '<tr><td colspan="5" style="color:#aaa;">No prophage regions detected</td></tr>'
        return f"""
<div class="sig-block">
  <h3>{_SIGNAL_LABELS['prophage']}</h3>
  <p>Source: PHASTER API. {ev.get('region_count', 0)} prophage region(s) detected.</p>
  {('<div class="warn-inline">' + s.warning + '</div>') if s.warning else ''}
  <table>
    <tr><th>Region</th><th>Completeness</th><th>Position (bp)</th>
        <th>GC%</th><th>CDS count</th></tr>
    {reg_rows}
  </table>
</div>"""

    # Local BLAST fallback or skipped
    hits = ev.get("top_hits", [])
    header = (
        "<tr><th>Subject</th><th>Description</th>"
        "<th>Identity</th><th>Coverage</th><th>E-value</th></tr>"
    )
    return f"""
<div class="sig-block">
  <h3>{_SIGNAL_LABELS['prophage']}</h3>
  <p>Source: local BLAST. {ev.get('hit_count', 0)} phage gene hits.</p>
  {('<div class="warn-inline">' + s.warning + '</div>') if s.warning else ''}
  <table>{header}{_hit_table(hits)}</table>
</div>"""


def _detail_skipped(s: SignalResult) -> str:
    label = _SIGNAL_LABELS.get(s.signal_name, s.signal_name)
    return f"""
<div class="sig-block" style="opacity:.5;">
  <h3>{label}</h3>
  <div class="warn-inline">{s.warning or 'Signal skipped.'}</div>
</div>"""


def _signal_detail(s: SignalResult) -> str:
    if s.skipped or s.score is None:
        return _detail_skipped(s)
    dispatch = {
        "gc_content":  _detail_gc,
        "is_proximity": lambda x: _detail_blast_simple(
            x, _SIGNAL_LABELS["is_proximity"]),
        "integron":    _detail_integron,
        "conjugative": _detail_conjugative,
        "prophage":    _detail_prophage,
    }
    fn = dispatch.get(s.signal_name)
    if fn:
        return fn(s)
    return f'<div class="sig-block"><h3>{s.signal_name}</h3><pre>{json.dumps(s.evidence, indent=2)}</pre></div>'


# ---------------------------------------------------------------------------
# Three-layer report section
# ---------------------------------------------------------------------------

_BAND_COLOUR = {
    "low":       "#2d7d2d",
    "moderate":  "#8a6000",
    "high":      "#b84000",
    "very_high": "#a00000",
}
_BAND_BG = {
    "low":       "#e8f5e8",
    "moderate":  "#fff8e0",
    "high":      "#fff0e8",
    "very_high": "#fce8e8",
}

_LAYER_LABELS = {
    "transfer_opportunity": "Transfer Opportunity",
    "establishment":        "Establishment",
    "consequence":          "Functional Consequence",
}

_FEATURE_DISPLAY_NAMES = {
    "is_element_match":      "IS / Mobile-element match",
    "integron_association":  "Integron / attC association",
    "conjugative_element":   "Conjugative element homology",
    "plasmid_context":       "Plasmid-context probability",
    "transposase_proximity": "Transposase proximity",
    "repeat_density":        "Flanking repeat density",
    "gc_deviation":          "GC content deviation",
    "codon_usage_distance":  "Codon usage distance",
    "taxonomic_distance":    "Taxonomic distance (donor–recipient)",
    "promoter_plausibility": "Promoter plausibility (proxy)",
    "sequence_complexity":   "Sequence complexity / length",
    "prophage_context":      "Prophage context",
    "amr_content":           "AMR gene content",
    "virulence_flags":       "Virulence / toxin flags",
    "gene_completeness":     "ORF completeness",
    "payload_count":         "Functional gene count",
}

_SOURCE_BADGE = {
    "signal_reuse": '<span style="font-size:10px;background:#e0e8f5;padding:1px 5px;border-radius:3px;">reused signal</span>',
    "computed":     '<span style="font-size:10px;background:#e8f5e0;padding:1px 5px;border-radius:3px;">computed</span>',
    "placeholder":  '<span style="font-size:10px;background:#f5f0e0;padding:1px 5px;border-radius:3px;">placeholder</span>',
}


def _feature_rows(features: list[FeatureResult], colour: str) -> str:
    rows = ""
    for f in features:
        label = _FEATURE_DISPLAY_NAMES.get(f.feature_name, f.feature_name)
        badge = _SOURCE_BADGE.get(f.source, "")
        if not f.available or f.score is None:
            score_cell = '<span style="color:#aaa;font-style:italic;">Unavailable</span>'
            dim = ' style="opacity:.45;"'
        else:
            score_cell = _bar(f.score, colour)
            dim = ""
        interp = f.interpretation[:100] + "…" if len(f.interpretation) > 100 else f.interpretation
        rows += (
            f"<tr{dim}>"
            f"<td>{label} {badge}</td>"
            f"<td>{score_cell}</td>"
            f"<td>{f.weight * 100:.0f}%</td>"
            f"<td style='font-size:11px;color:#666;'>{interp}</td>"
            f"</tr>\n"
        )
    return rows


def _layer_block(lr: LayerResult, profile_weight: float, colour: str) -> str:
    label = _LAYER_LABELS.get(lr.layer_name, lr.layer_name)
    completeness_pct = f"{lr.completeness:.0%}"
    rows = _feature_rows(lr.feature_results, colour)
    return f"""
<div class="sig-block">
  <h3 style="display:flex;justify-content:space-between;">
    <span>{label}</span>
    <span>
      {_bar(lr.layer_score, colour)}
      <span style="font-size:11px;color:#888;margin-left:12px;">
        profile weight {profile_weight * 100:.0f}% &nbsp;|&nbsp;
        {completeness_pct} features available
      </span>
    </span>
  </h3>
  <table style="margin-top:8px;">
    <tr>
      <th>Feature</th>
      <th>Score</th>
      <th>Layer weight</th>
      <th>Interpretation</th>
    </tr>
    {rows}
  </table>
</div>"""


def _three_layer_section(tl: ThreeLayerResult) -> str:
    band = tl.score_band.value
    colour = _BAND_COLOUR.get(band, "#666")
    bg     = _BAND_BG.get(band, "#f0f0f0")
    profile = tl.weight_profile

    layer_blocks = (
        _layer_block(tl.transfer_layer,     profile.get("transfer_opportunity", 0.4), colour)
        + _layer_block(tl.establishment_layer, profile.get("establishment", 0.35), colour)
        + _layer_block(tl.consequence_layer,   profile.get("consequence", 0.25), colour)
    )

    top_list = "".join(f"<li>{c}</li>" for c in tl.top_contributors) or "<li>(none above threshold)</li>"
    red_list = "".join(f"<li>{r}</li>" for r in tl.risk_reducers) or "<li>(none below threshold)</li>"
    miss_list = "".join(f"<li>{m}</li>" for m in tl.missing_important_features) or "<li>(none)</li>"

    completeness_colour = "#2d7d2d" if tl.overall_completeness >= 0.7 else "#b84000"

    return f"""
<!-- ===== Three-layer model ===== -->
<section id="three-layer">
  <h2>Three-Layer HGT Risk Assessment (v0.2)</h2>

  <div style="margin-bottom:16px;">
    <span class="badge"
          style="background:{bg};color:{colour};border:2px solid {colour};">
      {band.replace('_', ' ').upper()}
    </span>
    <span class="idx" style="color:{colour};">{tl.hgt_risk_index:.3f}</span>
    <span style="color:#888;font-size:12px;margin-left:8px;">HGT Risk Index / 1.000</span>
    <span style="margin-left:20px;font-size:12px;color:{completeness_colour};">
      &#9679; {tl.overall_completeness:.0%} features assessed
      &nbsp;|&nbsp; profile: <strong>{tl.weight_profile_name}</strong>
    </span>
  </div>

  <div class="sig-block" style="background:#fafbfc;">
    <h3 style="margin-bottom:8px;">Explanation</h3>
    <p style="line-height:1.7;">{tl.explanation}</p>
    <div style="display:flex;gap:24px;margin-top:14px;flex-wrap:wrap;">
      <div>
        <div style="font-size:11px;font-weight:700;text-transform:uppercase;letter-spacing:.4px;color:#888;margin-bottom:4px;">Top contributors</div>
        <ul style="margin-left:16px;font-size:12px;">{top_list}</ul>
      </div>
      <div>
        <div style="font-size:11px;font-weight:700;text-transform:uppercase;letter-spacing:.4px;color:#888;margin-bottom:4px;">Risk-reducing factors</div>
        <ul style="margin-left:16px;font-size:12px;">{red_list}</ul>
      </div>
      <div>
        <div style="font-size:11px;font-weight:700;text-transform:uppercase;letter-spacing:.4px;color:#888;margin-bottom:4px;">Important missing data</div>
        <ul style="margin-left:16px;font-size:12px;">{miss_list}</ul>
      </div>
    </div>
  </div>

  <h3 style="margin:16px 0 10px;">Layer Scores</h3>
  {layer_blocks}

  <p class="note" style="margin-top:8px;">
    Layer weights: Transfer {profile.get('transfer_opportunity', 0)*100:.0f}% /
    Establishment {profile.get('establishment', 0)*100:.0f}% /
    Consequence {profile.get('consequence', 0)*100:.0f}%.
    Unavailable features are excluded; remaining weights are re-normalised.
    Scores reflect sequence-level indicators only.
  </p>
</section>"""


# ---------------------------------------------------------------------------
# Main generators
# ---------------------------------------------------------------------------

def generate_html(result: PipelineResult, output_path: Path) -> None:
    agg = result.aggregation
    q = result.query
    h = result.host
    level = agg.risk_level
    colour = _RISK_COLOUR[level]
    bg = _RISK_BG[level]

    skipped_n = len(agg.skipped_signals)
    skipped_banner = ""
    if skipped_n:
        skipped_banner = f"""
<div class="alert">
  <strong>Warning:</strong> {skipped_n} of 5 signal(s) were skipped
  ({", ".join(agg.skipped_signals)}).
  The risk index is computed from {5 - skipped_n} signals
  (re-normalised weight sum = {agg.active_weight_sum:.2f}) and may
  <strong>underestimate true HGT risk</strong>.
  Install the missing BLAST databases or enable network access for a complete assessment.
</div>"""

    score_rows = _signal_rows(agg.signal_results, colour, agg.active_weight_sum)
    detail_sections = "\n".join(_signal_detail(s) for s in agg.signal_results)
    three_layer_section = (
        _three_layer_section(result.three_layer)
        if result.three_layer is not None
        else ""
    )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>HGT Risk Assessment — {q.identifier}</title>
<style>{_CSS}</style>
</head>
<body>
<div class="page">

<header>
  <h1>HGT Risk Assessment Report</h1>
  <div class="meta">
    Generated: {result.run_timestamp} &nbsp;|&nbsp;
    hgt-risk-assessor v{result.pipeline_version} &nbsp;|&nbsp;
    Query: {q.identifier}
  </div>
</header>

<main>

<!-- ===== Summary ===== -->
<section id="summary">
  <h2>Overall Risk Assessment</h2>
  <div>
    <span class="badge" style="background:{bg};color:{colour};border:2px solid {colour};">
      {level.value.upper()}
    </span>
    <span class="idx" style="color:{colour};">{agg.risk_index:.3f}</span>
    <span style="color:#888;font-size:12px;margin-left:8px;">/ 1.000</span>
  </div>
  {skipped_banner}
  <div class="cards">
    <div class="card">
      <div class="lbl">Active Signals</div>
      <div class="val">{5 - skipped_n} / 5</div>
    </div>
    <div class="card">
      <div class="lbl">Active Weight Sum</div>
      <div class="val">{agg.active_weight_sum:.2f}</div>
    </div>
    <div class="card">
      <div class="lbl">Query Length</div>
      <div class="val">{q.length:,} bp</div>
    </div>
    <div class="card">
      <div class="lbl">Query GC%</div>
      <div class="val">{q.gc_content:.1%}</div>
    </div>
    <div class="card">
      <div class="lbl">Host GC%</div>
      <div class="val">{h.gc_content:.1%}</div>
    </div>
  </div>
</section>

<!-- ===== Organism info ===== -->
<section id="organisms">
  <h2>Sequence and Host Information</h2>
  <table style="width:auto;min-width:500px;">
    <tr><th>Field</th><th>Query</th><th>Host</th></tr>
    <tr><td>Identifier</td>
        <td>{q.identifier}</td>
        <td>{h.identifier}</td></tr>
    <tr><td>Description</td>
        <td>{q.source_format.value.upper()} input</td>
        <td>{h.organism_name or '—'}</td></tr>
    <tr><td>GC Content</td>
        <td>{q.gc_content:.2%}</td>
        <td>{h.gc_content:.2%}</td></tr>
    <tr><td>Length</td>
        <td>{q.length:,} bp</td>
        <td>—</td></tr>
    <tr><td>GC Source</td>
        <td>computed</td>
        <td>{h.source}</td></tr>
  </table>
</section>

<!-- ===== Signal scores ===== -->
<section id="scores">
  <h2>Signal Scores</h2>
  <table>
    <tr>
      <th>Signal</th>
      <th>Score</th>
      <th>Weight</th>
      <th>Weighted Contribution</th>
    </tr>
    {score_rows}
    <tr class="total-row">
      <td>Overall Risk Index</td>
      <td></td>
      <td></td>
      <td>{agg.risk_index:.3f}</td>
    </tr>
  </table>
</section>

<!-- ===== Signal detail ===== -->
<section id="detail">
  <h2>Signal Detail</h2>
  {detail_sections}
</section>

{three_layer_section}

<!-- ===== Methodology ===== -->
<section id="methodology">
  <h2>Methodology and Thresholds</h2>
  <p>
    The risk index is a weighted sum of five normalised signal scores (each 0.0–1.0).
    When signals are skipped, remaining weights are re-normalised to sum to 1.0 so the
    index remains in range.  Skipped signals are flagged prominently; the result should
    be interpreted conservatively when the active weight sum is substantially below 1.0.
  </p>

  <h3 style="margin-top:14px;">Signal Weights</h3>
  <table style="width:auto;">
    <tr><th>Signal</th><th>Default Weight</th><th>Rationale</th></tr>
    <tr><td>IS Element Proximity</td><td>0.25</td>
        <td>IS elements are the primary drivers of intra- and inter-genomic mobilisation.</td></tr>
    <tr><td>Conjugative Element Homology</td><td>0.25</td>
        <td>Self-transmissibility is the most direct route to inter-organism HGT.</td></tr>
    <tr><td>Integron Association</td><td>0.20</td>
        <td>Integrons capture and disseminate gene cassettes; attC sites are direct evidence.</td></tr>
    <tr><td>GC Content Deviation</td><td>0.15</td>
        <td>Sequence composition anomalies are a classical marker of recent horizontal acquisition.</td></tr>
    <tr><td>Prophage Context</td><td>0.15</td>
        <td>Prophage induction can transduce flanking sequences into new host cells.</td></tr>
  </table>

  <h3 style="margin-top:14px;">Risk Classification</h3>
  <table style="width:auto;">
    <tr><th>Level</th><th>Index Range</th><th>Recommended Action</th></tr>
    <tr><td style="color:#2d7d2d;font-weight:700;">Low</td>
        <td>0.000 – 0.249</td>
        <td>Minimal sequence-level indicators. Standard laboratory biosafety procedures apply.</td></tr>
    <tr><td style="color:#8a6000;font-weight:700;">Medium</td>
        <td>0.250 – 0.499</td>
        <td>Some indicators present. Expert review recommended before scale-up or release.</td></tr>
    <tr><td style="color:#b84000;font-weight:700;">High</td>
        <td>0.500 – 0.749</td>
        <td>Multiple significant indicators. Formal contained use risk assessment required (UK SACGM).</td></tr>
    <tr><td style="color:#a00000;font-weight:700;">Critical</td>
        <td>0.750 – 1.000</td>
        <td>Strong HGT risk signals. Do not proceed without institutional biosafety officer review.</td></tr>
  </table>

  <h3 style="margin-top:14px;">Data Sources</h3>
  <ul>
    <li><strong>ISfinder</strong> — IS element reference database (isfinder.biotoul.fr)</li>
    <li><strong>INTEGRALL</strong> — integron and gene cassette reference (integrall.bio.ua.pt)</li>
    <li><strong>NCBI Protein</strong> — conjugative element protein sequences (relaxases, T4SS)</li>
    <li><strong>PHASTER</strong> — prophage identification API (phaster.ca)</li>
    <li><strong>NCBI RefSeq Viral</strong> — prophage BLAST fallback</li>
  </ul>

  <p style="margin-top:14px;font-size:11px;color:#888;">
    This tool provides sequence-level risk indicators only.  Output must be interpreted
    alongside experimental data and reviewed by a qualified biosafety officer before any
    regulatory submission.  This report does not constitute a formal risk assessment under
    UK GMO Contained Use Regulations or Scottish Environment Protection Agency requirements.
  </p>
</section>

</main>

<footer>
  hgt-risk-assessor v{result.pipeline_version} &nbsp;|&nbsp;
  {result.run_timestamp} &nbsp;|&nbsp;
  For research use only — not a substitute for expert biosafety review.
</footer>

</div>
</body>
</html>"""

    output_path.write_text(html, encoding="utf-8")


def generate_json(result: PipelineResult, output_path: Path) -> None:
    """Write a machine-readable JSON summary."""

    def _serial(obj):
        if hasattr(obj, "value"):           # Enum
            return obj.value
        if hasattr(obj, "__dataclass_fields__"):
            return {k: _serial(v) for k, v in vars(obj).items()}
        if isinstance(obj, list):
            return [_serial(i) for i in obj]
        if isinstance(obj, dict):
            return {k: _serial(v) for k, v in obj.items()}
        return obj

    output_path.write_text(
        json.dumps(_serial(result), indent=2),
        encoding="utf-8",
    )
