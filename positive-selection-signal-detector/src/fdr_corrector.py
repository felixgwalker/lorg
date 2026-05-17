import numpy as np


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
