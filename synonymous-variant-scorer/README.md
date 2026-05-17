# Synonymous Variant Scorer

## Project Overview

A functional impact scoring tool for synonymous (silent) variants. The tool challenges the embedded assumption in clinical and research genomics that synonymous mutations — variants that do not change the amino acid sequence — are functionally inert. By integrating evidence across four distinct biological mechanisms, it produces a per-variant score reflecting the real probability of functional consequence.

---

## Mechanisms Assessed

The pipeline evaluates four mechanisms through which a synonymous variant can have functional impact:

1. **Splicing Disruption** — synonymous variants within or near splice sites can alter splice site recognition, create cryptic splice sites, or disrupt exonic splicing enhancers/silencers. Scored using SpliceAI predictions.
2. **Codon Usage Bias** — substitution to a rare codon for the same amino acid can slow ribosomal elongation, affecting protein folding and yield. Scored against organism-specific codon usage frequency tables.
3. **mRNA Stability** — synonymous changes alter local mRNA secondary structure, affecting transcript half-life and translational efficiency. Scored using mRNA stability predictors.
4. **Cotranslational Folding** — ribosomal pausing at rare codons influences the folding trajectory of the nascent peptide. Variants that disrupt established pausing patterns can produce misfolded protein despite an identical final sequence. Scored using pause-site conservation and codon ramp models.

Each mechanism contributes a component score; scores are aggregated into a composite functional impact index.

---

## Data Sources

| Source | Use |
|---|---|
| gnomAD | Population frequency filtering and variant context |
| SpliceAI | Splicing disruption prediction |
| Codon usage tables (Kazusa / CoCoPUTs) | Codon frequency scoring per organism |
| mRNA stability predictors | Secondary structure and stability delta scoring |

---

## Technical Stack

| Component | Role |
|---|---|
| Python pipeline | VCF parsing, mechanism scoring, report generation |
| VCF input | Standard variant input format |
| Next.js report interface | Per-variant score breakdown, mechanism visualisation |
| VCF platform integration | Connects to existing variant interpretation platform work |

---

## Novel Contribution

Existing tools address individual mechanisms in isolation: SpliceAI for splicing, codon usage tools for elongation effects. No accessible, integrated tool currently scores synonymous variants across **all four functional mechanisms simultaneously** in a single workflow with a unified output format. This is the first such tool, and it is designed to be usable by clinical bioinformaticians without specialist knowledge of each underlying mechanism.

---

## Wider Relevance

- **Clinical genomics**: synonymous variants are systematically under-investigated in diagnostic pipelines; this tool provides a rapid triage score to prioritise follow-up.
- **Research variant interpretation**: functional genomics studies routinely discard synonymous variants; this tool enables their inclusion in mechanistic hypotheses.
- **Scientific narrative**: the project has a clear intellectual argument — that a widely held assumption (synonymous = silent) is empirically incomplete — which makes it a strong candidate for publication and public communication.
- Designed to integrate with existing VCF platform infrastructure.
