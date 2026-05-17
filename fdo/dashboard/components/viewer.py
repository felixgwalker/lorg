"""
3D molecular viewer component using 3Dmol.js embedded via Streamlit.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import streamlit.components.v1 as components

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

# Domain colour constants kept in sync with config
CAT_COLOR  = "#4A90D9"   # catalytic domain – blue
CTD_COLOR  = "#50C878"   # C-terminal domain – green
POCK_COLOR = "#FF6B35"   # pocket highlight – orange-red
BG_COLOR   = "0x0d1117"

# CDN for 3Dmol.js (pinned version)
_3DMOL_CDN = "https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.4/3Dmol-min.js"


def _escape(s: str) -> str:
    return s.replace("\\", "\\\\").replace("`", "\\`").replace("$", "\\$")


def render_structure(
    pdb_content: str,
    height: int = 520,
    spin: bool = False,
    highlight_residues: list[int] | None = None,
    pocket_residues: list[int] | None = None,
    cat_range: tuple[int, int] = (31, 326),
    cterm_range: tuple[int, int] = (327, 498),
    style: str = "cartoon",
    label: str = "",
) -> None:
    """
    Render a PDB structure in a dark-themed 3Dmol.js viewer.

    Parameters
    ----------
    pdb_content : str
        Full PDB file text.
    height : int
        Viewer height in pixels.
    spin : bool
        Auto-rotate the structure.
    highlight_residues : list[int]
        Residue numbers to highlight in yellow (e.g. active site).
    pocket_residues : list[int]
        Residue numbers to colour orange-red (pocket-lining residues).
    cat_range : tuple
        (start, end) residue numbers for the catalytic domain.
    cterm_range : tuple
        (start, end) residue numbers for the C-terminal domain.
    style : str
        "cartoon" | "surface" | "stick" | "sphere"
    label : str
        Optional label to overlay.
    """
    highlight_residues = highlight_residues or []
    pocket_residues    = pocket_residues    or []
    pdb_escaped        = _escape(pdb_content)

    spin_js    = "viewer.spin(true);" if spin else ""
    label_html = (
        f'<div style="position:absolute;top:6px;left:8px;color:#8b949e;'
        f'font-size:12px;font-family:monospace;pointer-events:none">{label}</div>'
        if label else ""
    )

    cat_sel   = f"{{resi: [{cat_range[0]}, {cat_range[1]}], chain: 'A'}}"
    cterm_sel = f"{{resi: [{cterm_range[0]}, {cterm_range[1]}], chain: 'A'}}"
    hl_sel    = f"{{resi: {highlight_residues}}}" if highlight_residues else "{}"
    pock_sel  = f"{{resi: {pocket_residues}}}"    if pocket_residues    else "{}"

    style_map = {
        "cartoon": (
            f"viewer.setStyle({cat_sel},   {{cartoon: {{color: '{CAT_COLOR}',  opacity: 0.92}}}});\n"
            f"    viewer.setStyle({cterm_sel}, {{cartoon: {{color: '{CTD_COLOR}', opacity: 0.92}}}});\n"
        ),
        "surface": (
            f"viewer.setStyle({cat_sel},   {{surface: {{color: '{CAT_COLOR}',  opacity: 0.5}}, "
            f"cartoon: {{color: '{CAT_COLOR}'}}}});\n"
            f"    viewer.setStyle({cterm_sel}, {{surface: {{color: '{CTD_COLOR}', opacity: 0.5}}, "
            f"cartoon: {{color: '{CTD_COLOR}'}}}});\n"
        ),
        "stick":   "viewer.setStyle({}, {stick: {colorscheme: 'chain'}});\n",
        "sphere":  "viewer.setStyle({}, {sphere: {radius: 0.4, colorscheme: 'chain'}});\n",
    }

    base_style = style_map.get(style, style_map["cartoon"])

    hl_js = (
        f"viewer.setStyle({hl_sel}, {{sphere: {{color: '#FFD700', radius: 0.9}}, "
        f"cartoon: {{color: '#FFD700'}}}});\n"
        if highlight_residues else ""
    )
    pock_js = (
        f"viewer.setStyle({pock_sel}, {{sphere: {{color: '{POCK_COLOR}', radius: 0.8}}, "
        f"cartoon: {{color: '{POCK_COLOR}'}}}});\n"
        if pocket_residues else ""
    )

    html = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<script src="{_3DMOL_CDN}"></script>
<style>
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  html, body {{ background: #0d1117; width: 100%; height: {height}px; overflow: hidden; }}
  #viewer-container {{
    width: 100%; height: {height}px; position: relative;
    border: 1px solid #30363d; border-radius: 8px; overflow: hidden;
  }}
  .legend {{
    position: absolute; bottom: 8px; left: 8px;
    display: flex; gap: 10px; font-size: 11px;
    font-family: -apple-system, sans-serif; color: #f0f6fc;
    background: rgba(13,17,23,0.75); padding: 4px 8px; border-radius: 4px;
    pointer-events: none;
  }}
  .legend span {{ display: flex; align-items: center; gap: 4px; }}
  .dot {{ width: 8px; height: 8px; border-radius: 50%; flex-shrink: 0; }}
</style>
</head>
<body>
<div id="viewer-container">
  {label_html}
  <div id="mol-viewer" style="width:100%;height:{height}px;"></div>
  <div class="legend">
    <span><div class="dot" style="background:{CAT_COLOR}"></div>Catalytic</span>
    <span><div class="dot" style="background:{CTD_COLOR}"></div>C-terminal</span>
    {'<span><div class="dot" style="background:' + POCK_COLOR + '"></div>Pocket</span>' if pocket_residues else ''}
    {'<span><div class="dot" style="background:#FFD700"></div>Active site</span>' if highlight_residues else ''}
  </div>
</div>

<script>
(function() {{
  var element = document.getElementById("mol-viewer");
  var config  = {{ backgroundColor: "{BG_COLOR}", antialias: true }};
  var viewer  = $3Dmol.createViewer(element, config);
  var pdbData = `{pdb_escaped}`;

  viewer.addModel(pdbData, "pdb");
  viewer.setStyle({{}}, {{cartoon: {{color: "#555555"}}}});

  {base_style}
  {hl_js}
  {pock_js}

  viewer.zoomTo();
  {spin_js}
  viewer.render();
}})();
</script>
</body>
</html>"""

    components.html(html, height=height + 4, scrolling=False)


def render_pocket_sphere(
    pdb_content: str,
    pocket: dict,
    height: int = 480,
    cat_range: tuple[int, int] = (31, 326),
    cterm_range: tuple[int, int] = (327, 498),
) -> None:
    """
    Render a structure with the pocket shown as a translucent sphere
    centred on the pocket centroid, plus highlighted lining residues.
    """
    center     = pocket.get("center", [0, 0, 0])
    vol        = pocket.get("mean_volume_A3", 300)
    pock_color = pocket.get("color", POCK_COLOR)
    pdb_escaped = _escape(pdb_content)

    # Sphere radius estimated from pocket volume (V = 4/3 π r³)
    radius = max(4.0, (3 * vol / (4 * math.pi)) ** (1 / 3))

    cx, cy, cz = float(center[0]), float(center[1]), float(center[2])

    pocket_res = pocket.get("residues", [])
    pock_sel   = f"{{resi: {pocket_res}}}" if pocket_res else "{}"
    pock_js    = (
        f"viewer.setStyle({pock_sel}, "
        f"{{sphere: {{color: '{pock_color}', radius: 0.8}}, "
        f"cartoon: {{color: '{pock_color}'}}}});\n"
        if pocket_res else ""
    )

    label_text = f"Pocket {pocket.get('id', '')} · vol {vol:.0f} Å³ · radius {radius:.1f} Å"

    cat_sel   = f"{{resi: [{cat_range[0]}, {cat_range[1]}], chain: 'A'}}"
    cterm_sel = f"{{resi: [{cterm_range[0]}, {cterm_range[1]}], chain: 'A'}}"

    html = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<script src="{_3DMOL_CDN}"></script>
<style>
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  html, body {{ background: #0d1117; width: 100%; height: {height}px; overflow: hidden; }}
  #viewer-container {{
    width: 100%; height: {height}px; position: relative;
    border: 1px solid #30363d; border-radius: 8px; overflow: hidden;
  }}
  .label {{
    position: absolute; top: 6px; left: 8px; color: #8b949e;
    font-size: 12px; font-family: monospace; pointer-events: none;
  }}
  .legend {{
    position: absolute; bottom: 8px; left: 8px;
    display: flex; gap: 10px; font-size: 11px;
    font-family: -apple-system, sans-serif; color: #f0f6fc;
    background: rgba(13,17,23,0.75); padding: 4px 8px; border-radius: 4px;
    pointer-events: none;
  }}
  .legend span {{ display: flex; align-items: center; gap: 4px; }}
  .dot {{ width: 8px; height: 8px; border-radius: 50%; flex-shrink: 0; }}
</style>
</head>
<body>
<div id="viewer-container">
  <div class="label">{label_text}</div>
  <div id="mol-viewer" style="width:100%;height:{height}px;"></div>
  <div class="legend">
    <span><div class="dot" style="background:{CAT_COLOR}"></div>Catalytic</span>
    <span><div class="dot" style="background:{CTD_COLOR}"></div>C-terminal</span>
    <span><div class="dot" style="background:{pock_color}"></div>Pocket lining</span>
    <span><div class="dot" style="background:{pock_color};opacity:0.35;border:1px solid {pock_color}"></div>Pocket volume</span>
  </div>
</div>

<script>
(function() {{
  var element = document.getElementById("mol-viewer");
  var viewer  = $3Dmol.createViewer(element,
    {{ backgroundColor: "{BG_COLOR}", antialias: true }});
  var pdbData = `{pdb_escaped}`;

  viewer.addModel(pdbData, "pdb");

  // Domain colouring
  viewer.setStyle({{}}, {{cartoon: {{color: "#555555"}}}});
  viewer.setStyle({cat_sel},   {{cartoon: {{color: "{CAT_COLOR}",  opacity: 0.90}}}});
  viewer.setStyle({cterm_sel}, {{cartoon: {{color: "{CTD_COLOR}",  opacity: 0.90}}}});

  // Highlight pocket-lining residues
  {pock_js}

  // Translucent sphere at pocket centroid
  viewer.addSphere({{
    center: {{ x: {cx}, y: {cy}, z: {cz} }},
    radius: {radius:.2f},
    color:  "{pock_color}",
    opacity: 0.22,
    wireframe: false,
  }});

  viewer.zoomTo();
  viewer.render();
}})();
</script>
</body>
</html>"""

    components.html(html, height=height + 4, scrolling=False)


def render_ghost_trail(
    frames_pdb: list[str],
    height: int = 520,
    cat_range: tuple[int, int] = (31, 326),
    cterm_range: tuple[int, int] = (327, 498),
) -> None:
    """
    Render multiple frames as superimposed semi-transparent cartoons
    ('ghost trail') to visualise domain motion.
    """
    if not frames_pdb:
        return

    n = len(frames_pdb)
    models_js = ""
    for i, pdb in enumerate(frames_pdb):
        opacity = 0.15 + 0.7 * (i / max(n - 1, 1))
        esc = _escape(pdb)
        models_js += f"""
    viewer.addModel(`{esc}`, "pdb");
    viewer.setStyle({{model: {i}}}, {{cartoon: {{
        color: '{CAT_COLOR}', opacity: {opacity:.2f}
    }}}});
"""

    html = f"""<!DOCTYPE html>
<html>
<head>
<script src="{_3DMOL_CDN}"></script>
<style>
  html, body {{ background: #0d1117; margin: 0; width: 100%; height: {height}px; }}
  #viewer {{ width: 100%; height: {height}px; border: 1px solid #30363d; border-radius: 8px; }}
</style>
</head>
<body>
<div id="viewer"></div>
<script>
(function() {{
  var viewer = $3Dmol.createViewer(document.getElementById("viewer"),
    {{ backgroundColor: "{BG_COLOR}", antialias: true }});
  {models_js}
  viewer.zoomTo();
  viewer.render();
}})();
</script>
</body>
</html>"""
    components.html(html, height=height + 4, scrolling=False)
