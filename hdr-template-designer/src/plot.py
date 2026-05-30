"""Annotated HDR template schematic diagram."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from pathlib import Path


def _gc_sliding(seq: str, window: int = 10) -> list[float]:
    """Return per-position GC fraction using a sliding window."""
    if not seq:
        return []
    result = []
    half = window // 2
    for i in range(len(seq)):
        lo = max(0, i - half)
        hi = min(len(seq), i + half + 1)
        chunk = seq[lo:hi]
        result.append((chunk.count("G") + chunk.count("C")) / len(chunk))
    return result


def plot_template_diagram(
    template: dict,
    qc_checks: list[dict],
    output_dir: Path,
    fmt: str = "png",
    pam_mutations: list[dict] | None = None,
) -> Path:
    """Render an annotated schematic of the HDR donor template.

    The figure contains two stacked panels:

    1. **Block diagram** – coloured boxes for left arm, edit, and right arm,
       with PAM-mutation sites marked by vertical tick marks.
    2. **GC content track** – sliding-window GC % plotted along the full
       template length, with a shaded 40–60 % optimal band.

    Parameters
    ----------
    template : dict
        Template dict from ``build_template``.
    qc_checks : list[dict]
        Per-check results from ``run_all_checks``.
    output_dir : Path
        Directory to write the diagram into.
    fmt : str
        Image format, ``"png"`` or ``"svg"``.
    pam_mutations : list[dict] | None
        Optional PAM-disruption mutations from ``suggest_pam_disruption``
        (keys: ``pos``, ``ref``, ``alt``).

    Returns
    -------
    Path
        Path to the saved figure.
    """
    out_path = output_dir / f"template_diagram.{fmt}"

    left_arm = template["left_arm"]
    edit_seq = template["edit_seq"]
    right_arm = template["right_arm"]
    left_len = len(left_arm)
    edit_len = max(len(edit_seq), 1)
    right_len = len(right_arm)
    total = left_len + edit_len + right_len

    left_gc = (left_arm.count("G") + left_arm.count("C")) / max(len(left_arm), 1)
    right_gc = (right_arm.count("G") + right_arm.count("C")) / max(len(right_arm), 1)

    left_flag = any(c["flag"] for c in qc_checks if "left" in c.get("check", ""))
    right_flag = any(c["flag"] for c in qc_checks if "right" in c.get("check", ""))

    left_color = "#FF7043" if left_flag else "#42A5F5"
    right_color = "#FF7043" if right_flag else "#66BB6A"
    edit_color = "#FDD835"

    # ------------------------------------------------------------------ #
    #  Figure layout: two panels (block diagram + GC track)              #
    # ------------------------------------------------------------------ #
    fig, (ax_blocks, ax_gc) = plt.subplots(
        2, 1,
        figsize=(13, 5),
        gridspec_kw={"height_ratios": [3, 2], "hspace": 0.45},
    )

    # ------------------------------------------------------------------ #
    #  Panel 1: Block diagram                                             #
    # ------------------------------------------------------------------ #
    ax_blocks.set_xlim(0, total)
    ax_blocks.set_ylim(0, 2.4)
    ax_blocks.axis("off")

    # Left arm block
    rect_left = mpatches.FancyBboxPatch(
        (0, 0.6), left_len, 1.0,
        boxstyle="round,pad=0.02",
        facecolor=left_color, edgecolor="white",
        linewidth=1.5, alpha=0.85,
    )
    ax_blocks.add_patch(rect_left)
    ax_blocks.text(
        left_len / 2, 1.1,
        f"Left Arm\n{left_len} bp\nGC={left_gc:.0%}",
        ha="center", va="center", fontsize=9, fontweight="bold", color="white",
    )

    # Edit block
    edit_x = left_len
    # Scale the edit block width for visibility — at least 1% of total width.
    display_edit_len = max(edit_len, max(1, int(total * 0.015)))
    rect_edit = mpatches.FancyBboxPatch(
        (edit_x, 0.5), display_edit_len, 1.2,
        boxstyle="round,pad=0.02",
        facecolor=edit_color, edgecolor="#F57F17",
        linewidth=2.0, alpha=0.95,
    )
    ax_blocks.add_patch(rect_edit)
    label_edit = edit_seq if len(edit_seq) <= 6 else f"{edit_seq[:6]}…"
    ax_blocks.text(
        edit_x + display_edit_len / 2, 1.1,
        f"Edit\n{label_edit}",
        ha="center", va="center", fontsize=8, fontweight="bold", color="#4E342E",
    )

    # Right arm block (offset by display_edit_len so it doesn't overlap)
    right_x = edit_x + display_edit_len
    display_right_len = total - left_len - display_edit_len
    if display_right_len > 0:
        rect_right = mpatches.FancyBboxPatch(
            (right_x, 0.6), display_right_len, 1.0,
            boxstyle="round,pad=0.02",
            facecolor=right_color, edgecolor="white",
            linewidth=1.5, alpha=0.85,
        )
        ax_blocks.add_patch(rect_right)
        ax_blocks.text(
            right_x + display_right_len / 2, 1.1,
            f"Right Arm\n{right_len} bp\nGC={right_gc:.0%}",
            ha="center", va="center", fontsize=9, fontweight="bold", color="white",
        )

    # Boundary tick marks and labels
    for bx, label in [(left_len, f"{left_len}"), (left_len + edit_len, f"{left_len + edit_len}")]:
        ax_blocks.axvline(x=bx, ymin=0.25, ymax=0.85, color="#5D4037", linewidth=1.2,
                          linestyle="--", alpha=0.7)

    # PAM mutation tick marks
    if pam_mutations:
        for mut in pam_mutations:
            px = mut["pos"]
            if 0 <= px < total:
                ax_blocks.annotate(
                    f"{mut['ref']}->{mut['alt']}",
                    xy=(px, 1.6),
                    xytext=(px, 2.1),
                    ha="center", va="bottom", fontsize=7, color="#B71C1C",
                    arrowprops=dict(arrowstyle="-|>", color="#B71C1C", lw=1.0),
                )

    # Ruler ticks along x axis (block panel)
    tick_interval = max(1, round(total / 8 / 10) * 10) if total > 20 else 5
    tick_positions = list(range(0, total + 1, tick_interval))
    for tp in tick_positions:
        ax_blocks.plot([tp, tp], [0.5, 0.6], color="#757575", lw=0.8)
        ax_blocks.text(tp, 0.42, str(tp), ha="center", va="top", fontsize=7, color="#757575")

    # QC flag banner
    flags = [c for c in qc_checks if c["flag"]]
    if flags:
        flag_text = "QC WARNINGS: " + " | ".join(c["check"] for c in flags)
        ax_blocks.text(
            total / 2, 0.15, flag_text,
            ha="center", va="center", fontsize=8, color="#B71C1C",
            bbox=dict(boxstyle="round", facecolor="#FFCDD2", edgecolor="#B71C1C", alpha=0.85),
        )

    legend_handles = [
        mpatches.Patch(facecolor="#42A5F5", label="Left arm (pass)"),
        mpatches.Patch(facecolor="#66BB6A", label="Right arm (pass)"),
        mpatches.Patch(facecolor="#FF7043", label="Arm (QC flag)"),
        mpatches.Patch(facecolor="#FDD835", edgecolor="#F57F17", label="Edit sequence"),
    ]
    if pam_mutations:
        legend_handles.append(
            mpatches.Patch(facecolor="#FFCDD2", edgecolor="#B71C1C", label="PAM mutation")
        )
    ax_blocks.legend(
        handles=legend_handles, loc="upper right", fontsize=8,
        bbox_to_anchor=(1.0, 2.4), framealpha=0.9,
    )
    ax_blocks.set_title("HDR Donor Template Schematic", fontsize=13, fontweight="bold", pad=8)

    # ------------------------------------------------------------------ #
    #  Panel 2: GC content sliding window track                          #
    # ------------------------------------------------------------------ #
    full_seq = template["full_template"]
    window = max(10, len(full_seq) // 30)
    gc_vals = _gc_sliding(full_seq, window=window)

    x_pos = list(range(len(gc_vals)))
    ax_gc.fill_between(x_pos, 0.40, 0.60, color="#E8F5E9", alpha=0.6, label="Optimal GC (40–60%)")
    ax_gc.plot(x_pos, gc_vals, color="#1565C0", linewidth=1.2, label=f"GC% (w={window}bp)")
    ax_gc.axhline(0.40, color="#388E3C", linewidth=0.8, linestyle=":")
    ax_gc.axhline(0.60, color="#388E3C", linewidth=0.8, linestyle=":")

    # Shade arm/edit regions differently.
    ax_gc.axvspan(0, left_len, alpha=0.07, color="#42A5F5")
    ax_gc.axvspan(left_len, left_len + edit_len, alpha=0.18, color="#FDD835")
    ax_gc.axvspan(left_len + edit_len, len(full_seq), alpha=0.07, color="#66BB6A")

    # PAM mutation markers on GC track
    if pam_mutations:
        for mut in pam_mutations:
            px = mut["pos"]
            if 0 <= px < len(gc_vals):
                ax_gc.axvline(x=px, color="#B71C1C", linewidth=1.0, linestyle="--", alpha=0.7)

    ax_gc.set_xlim(0, len(full_seq))
    ax_gc.set_ylim(0, 1.05)
    ax_gc.set_xlabel("Template position (bp)", fontsize=9)
    ax_gc.set_ylabel("GC fraction", fontsize=9)
    ax_gc.set_title("GC Content Along Template", fontsize=10, fontweight="bold")
    ax_gc.legend(loc="upper right", fontsize=8, framealpha=0.9)
    ax_gc.tick_params(labelsize=8)
    ax_gc.spines["top"].set_visible(False)
    ax_gc.spines["right"].set_visible(False)

    # ------------------------------------------------------------------ #
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path
