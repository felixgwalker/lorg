"""Neighbor-joining tree construction using Biopython DistanceTreeConstructor."""

from __future__ import annotations

import numpy as np


def build_nj_tree(names: list[str], matrix: np.ndarray) -> str:
    """Build a Neighbor-Joining tree and return Newick string.

    Uses Biopython's DistanceTreeConstructor with the NJ algorithm.
    Falls back to a simple star topology Newick if Biopython NJ fails.
    """
    try:
        from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
        lower_tri: list[list[float]] = []
        for i in range(len(names)):
            row: list[float] = []
            for j in range(i + 1):
                row.append(float(matrix[i, j]))
            lower_tri.append(row)
        dm = DistanceMatrix(names, lower_tri)
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(dm)
        from io import StringIO
        from Bio import Phylo
        sio = StringIO()
        Phylo.write(tree, sio, "newick")
        return sio.getvalue().strip()
    except Exception:
        return _star_newick(names, matrix)


def _star_newick(names: list[str], matrix: np.ndarray) -> str:
    """Fallback: produce a star-topology Newick from mean distances."""
    n = len(names)
    if n == 0:
        return "();"
    parts = []
    for i, name in enumerate(names):
        mean_dist = matrix[i, :].sum() / max(n - 1, 1)
        safe_name = name.replace("(", "_").replace(")", "_").replace(",", "_").replace(":", "_").replace(";", "_")
        parts.append(f"{safe_name}:{mean_dist:.6f}")
    return "(" + ",".join(parts) + ");"
