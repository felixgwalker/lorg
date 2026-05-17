"""Write pairwise distance matrix CSV, Newick tree, and ranked comparison table."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def write_distance_matrix(names: list[str], matrix: np.ndarray, output_dir: Path) -> Path:
    """Write pairwise distance matrix as a labeled CSV."""
    df = pd.DataFrame(matrix, index=names, columns=names)
    out = output_dir / "distance_matrix.csv"
    df.to_csv(out)
    return out


def write_newick(newick: str, output_dir: Path) -> Path:
    """Write Newick tree string to a .nwk file."""
    out = output_dir / "phylogenetic_tree.nwk"
    out.write_text(newick + "\n")
    return out


def write_ranked_comparisons(
    names: list[str],
    matrix: np.ndarray,
    normalized_matrix: np.ndarray,
    output_dir: Path,
) -> Path:
    """Write ranked pairwise species comparison table CSV."""
    rows = []
    n = len(names)
    for i in range(n):
        for j in range(i + 1, n):
            rows.append({
                "species_1": names[i],
                "species_2": names[j],
                "distance": round(float(matrix[i, j]), 6),
                "normalized_distance": round(float(normalized_matrix[i, j]), 6),
            })
    df = pd.DataFrame(rows)
    df.sort_values("distance", inplace=True)
    df.reset_index(drop=True, inplace=True)
    df.insert(0, "rank", df.index + 1)
    out = output_dir / "ranked_comparisons.csv"
    df.to_csv(out, index=False)
    return out
