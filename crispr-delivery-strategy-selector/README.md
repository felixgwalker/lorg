# CRISPR Delivery Strategy Selector

Ranks CRISPR delivery modalities for a given cell type, payload, and experimental context.

## Overview

Selecting the right delivery method (plasmid, RNP, LNP, AAV, lentivirus) is as important as guide
design. Each modality has different efficiency, immunogenicity, integration risk, and payload size
constraints across cell types. This tool applies a rule-based compatibility scoring table to rank
delivery strategies and flag contraindications.

## Approach

**Inputs:** Delivery specification JSON with keys: `cell_type` (e.g. HEK293, primary_T, neuron,
hepatocyte), `organism` (human/mouse), `payload_type` (RNP/plasmid/mRNA), `payload_size_bp` (int),
`target_tissue` (optional), `allow_integration` (bool).

**Core method:** Five delivery modalities are scored:
- **Plasmid transfection**: high efficiency in dividing cells, low in primary/non-dividing cells.
- **RNP electroporation**: transient, low immunogenicity, effective in primary cells including T cells and HSCs.
- **LNP (mRNA+sgRNA)**: liver/lung tropism, no integration risk, limited payload.
- **AAV**: high efficiency in post-mitotic cells (neurons, hepatocytes), limited to 4.7 kb payload.
- **Lentivirus**: stable integration (flagged when `allow_integration` is False), broad cell tropism.
Scores are normalised to 0–1 and include size constraint checks.

**Outputs:** Strategy ranking TSV; notes JSON; optional bar chart.

**How it ships:** `python run_selector.py delivery_spec.json`; delegated from
`main.py → src.pipeline.main() → run_selector.py`.

## Usage

```bash
python run_selector.py delivery_spec.json -o results/
python run_selector.py --demo -o results/
```

## Dependencies

- pandas>=2.0.0
- numpy>=1.24.0
- matplotlib>=3.7.0

## Status

Planned
