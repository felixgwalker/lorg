"""Annotated HDR template diagram."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path


def plot_template_diagram(template: dict, qc_checks: list[dict], output_dir: Path, fmt: str = "png") -> Path:
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

    fig, ax = plt.subplots(figsize=(12, 3))
    ax.set_xlim(0, total)
    ax.set_ylim(0, 2)
    ax.axis("off")

    # Draw left arm
    rect_left = mpatches.FancyBboxPatch((0, 0.5), left_len, 1.0,
                                         boxstyle="round,pad=0.02",
                                         facecolor=left_color, edgecolor="white",
                                         linewidth=1.5, alpha=0.85)
    ax.add_patch(rect_left)
    ax.text(left_len / 2, 1.0, f"Left Arm\n{left_len} bp\nGC={left_gc:.0%}",
            ha="center", va="center", fontsize=9, fontweight="bold", color="white")

    # Draw edit
    edit_x = left_len
    rect_edit = mpatches.FancyBboxPatch((edit_x, 0.4), edit_len, 1.2,
                                         boxstyle="round,pad=0.02",
                                         facecolor=edit_color, edgecolor="#F57F17",
                                         linewidth=2.0, alpha=0.95)
    ax.add_patch(rect_edit)
    ax.text(edit_x + edit_len / 2, 1.0, f"Edit\n{edit_seq}",
            ha="center", va="center", fontsize=8, fontweight="bold", color="#4E342E")

    # Draw right arm
    right_x = edit_x + edit_len
    rect_right = mpatches.FancyBboxPatch((right_x, 0.5), right_len, 1.0,
                                          boxstyle="round,pad=0.02",
                                          facecolor=right_color, edgecolor="white",
                                          linewidth=1.5, alpha=0.85)
    ax.add_patch(rect_right)
    ax.text(right_x + right_len / 2, 1.0, f"Right Arm\n{right_len} bp\nGC={right_gc:.0%}",
            ha="center", va="center", fontsize=9, fontweight="bold", color="white")

    # QC flag indicators
    flags = [c for c in qc_checks if c["flag"]]
    if flags:
        flag_text = "QC WARNINGS: " + " | ".join(c["check"] for c in flags)
        ax.text(total / 2, 0.2, flag_text, ha="center", va="center",
                fontsize=8, color="#B71C1C",
                bbox=dict(boxstyle="round", facecolor="#FFCDD2", edgecolor="#B71C1C", alpha=0.8))

    legend_handles = [
        mpatches.Patch(facecolor="#42A5F5", label="Left arm (pass)"),
        mpatches.Patch(facecolor="#66BB6A", label="Right arm (pass)"),
        mpatches.Patch(facecolor="#FF7043", label="Arm (QC flag)"),
        mpatches.Patch(facecolor="#FDD835", edgecolor="#F57F17", label="Edit sequence"),
    ]
    ax.legend(handles=legend_handles, loc="upper right", fontsize=8,
              bbox_to_anchor=(1.0, 2.0), framealpha=0.9)

    ax.set_title("HDR Donor Template Diagram", fontsize=13, fontweight="bold", pad=10)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path
