# Gene Regulatory Network Builder

Infers gene regulatory network edges between transcription factors and target genes from expression data using GENIE3, ARACNE, or Pearson correlation.

## Overview

Given a gene expression matrix and an optional TF list, this tool infers directed regulatory edges from regulators to target genes and reports a ranked edge list with hub gene scores.

## Approach

**Inputs:** Gene expression matrix (genes × samples, TPM or log-normalised counts); optional text file of TF gene IDs.

**Core method:** **GENIE3** trains a random forest to predict each target gene's expression from all regulators (or all genes if no TF list), extracting feature importances as edge weights. **ARACNE** computes pairwise mutual information between all regulator-target pairs and applies the data processing inequality to prune indirect edges. **Correlation** uses absolute Pearson r between regulator-target pairs, with FDR correction. For all methods, the top `n` edges per target gene are retained (default 10). Hub scores are computed as the normalised out-degree (for TFs) or in-degree (for targets).

**Outputs:** TSV of ranked edges (`regulatory_edges.tsv`); TSV of node hub scores (`network_nodes.tsv`); optional network graph (for small networks).

**How it ships:** `python run_builder.py --expression expression.tsv --tf-list tfs.txt`; `main.py` delegates to `src.pipeline.main()` which loads `run_builder.py` via `importlib`.

## Usage

```bash
# Build a regulatory network with GENIE3
python run_builder.py --expression expression.tsv --tf-list tfs.txt -o results/

# Synthetic demo (no real input required)
python run_builder.py --demo -o results/

# Use ARACNE
python run_builder.py --expression expression.tsv --tf-list tfs.txt --method ARACNE -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- scipy>=1.10.0
- matplotlib>=3.7.0

## Status

Planned
