import numpy as np


def fdr_correct(p_values, method="bh"):
    """Apply multiple-testing correction to a list of p-values.

    Parameters
    ----------
    p_values : list of float
        Raw p-values to correct.
    method : str
        Correction method.  Currently only Benjamini-Hochberg ('bh') is
        supported.  Case-insensitive.

    Returns
    -------
    list of float
        FDR-adjusted p-values in the same order as the input.

    Raises
    ------
    ValueError
        If an unrecognised method is requested.
    """
    if method.lower() in ("bh", "benjamini-hochberg", "benjamini_hochberg"):
        return benjamini_hochberg(p_values)
    raise ValueError(
        f"Unknown FDR method '{method}'. Supported methods: 'bh' (Benjamini-Hochberg)."
    )


def benjamini_hochberg(pvalues):
    n = len(pvalues)
    if n == 0:
        return []
    indexed = sorted(enumerate(pvalues), key=lambda x: x[1])
    adjusted = [0.0] * n
    prev_adj = 1.0
    for rank, (orig_idx, pval) in enumerate(reversed(indexed), 1):
        adj = pval * n / (n - rank + 1)
        adj = min(adj, prev_adj)
        prev_adj = adj
        adjusted[orig_idx] = adj
    return adjusted


def apply_fdr(tested_genes, alpha=0.05):
    pvals = [g["lrt_pval"] for g in tested_genes]
    adj_pvals = benjamini_hochberg(pvals)
    for i, gene in enumerate(tested_genes):
        gene["adjusted_pval"] = adj_pvals[i]
        gene["significant"] = adj_pvals[i] < alpha and gene["selection_signal"]
    return tested_genes
