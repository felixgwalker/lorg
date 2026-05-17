import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import os


def plot_conservation_heatmap(scored_elements, species_list, output_dir):
    if not scored_elements:
        return []

    element_ids = [e["element_id"] for e in scored_elements]
    data = np.zeros((len(element_ids), len(species_list)))
    for i, elem in enumerate(scored_elements):
        for j, sp in enumerate(species_list):
            data[i, j] = elem["per_species_identity"].get(sp, 0.0)

    fig, ax = plt.subplots(figsize=(max(6, len(species_list) * 1.8), max(5, len(element_ids) * 0.55)))
    im = ax.imshow(data, aspect="auto", cmap="RdYlGn", vmin=0, vmax=1)

    ax.set_xticks(range(len(species_list)))
    ax.set_xticklabels(species_list, rotation=30, ha="right", fontsize=10)
    ax.set_yticks(range(len(element_ids)))
    ax.set_yticklabels(element_ids, fontsize=9)
    ax.set_title("Regulatory Element Conservation\n(Sequence Identity per Species)", fontsize=12)

    for i in range(len(element_ids)):
        for j in range(len(species_list)):
            val = data[i, j]
            color = "black" if 0.3 < val < 0.8 else "white"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7, color=color)

    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.04)
    cbar.set_label("Sequence Identity")

    classification_colors = {
        "conserved": "#2ca02c",
        "partially_conserved": "#ff7f0e",
        "diverged": "#d62728",
    }
    for i, elem in enumerate(scored_elements):
        cl = elem["classification"]
        ax.get_yticklabels()[i].set_color(classification_colors.get(cl, "black"))

    from matplotlib.patches import Patch
    legend_handles = [
        Patch(facecolor=v, label=k.replace("_", " ").title())
        for k, v in classification_colors.items()
    ]
    ax.legend(handles=legend_handles, loc="upper left", bbox_to_anchor=(1.15, 1),
              title="Classification", fontsize=8)

    plt.tight_layout()

    paths = []
    for ext in ("png", "svg"):
        out = os.path.join(output_dir, f"conservation_heatmap.{ext}")
        fig.savefig(out, dpi=150, bbox_inches="tight")
        paths.append(out)
    plt.close(fig)
    return paths
