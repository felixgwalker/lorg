# Genomic Resurrection Scorer

## Project Overview

A feasibility assessment tool for de-extinction projects. Given an ancient DNA dataset and a proxy living relative genome, the tool produces a structured, quantitative report on how tractable resurrection of the extinct species would be. The aim is to move de-extinction feasibility discussions from narrative speculation to reproducible, evidence-based scoring.

---

## Assessment Layers

The pipeline evaluates five layers in sequence:

1. **Ancient DNA Quality** — damage patterns, read depth, fragment length distribution, contamination estimates (MapDamage, GATK).
2. **Genomic Completeness** — proportion of the extinct genome recoverable above confidence thresholds relative to the proxy reference.
3. **Divergence Quantification** — genome-wide and functional-region divergence between extinct species and proxy; coding, regulatory, and non-coding partitions scored separately.
4. **Edit Burden Estimation** — number and complexity of edits required to convert the proxy genome toward the extinct genome; accounts for SNPs, indels, and structural variants.
5. **Ethical and Ecological Flagging** — structured flag generation covering habitat availability, ecological role, welfare considerations, and existing regulatory or conservation conflicts.

Each layer produces a sub-score; sub-scores are combined into an overall feasibility index with documented weighting logic.

---

## Initial Case Study

**Thylacine (*Thylacinus cynocephalus*) vs Tasmanian Devil (*Sarcophilus harrisii*)**

The initial case study uses the University of Melbourne 2022 thylacine genome publication as the ancient DNA source. This provides a well-characterised dataset with published quality metrics against which the pipeline output can be validated. The Tasmanian devil serves as the proxy living relative for divergence and edit burden calculations.

---

## Technical Stack

| Component | Role |
|---|---|
| Python pipeline | Orchestration, scoring logic, report generation |
| MapDamage | Ancient DNA damage assessment |
| GATK | Variant calling and genotype refinement |
| Ensembl VEP | Functional annotation of divergent variants |
| Custom scoring logic | Layer aggregation, weighting, feasibility index |
| Next.js front end | Interactive report interface, score visualisation |

---

## Novel Contribution

The individual components (MapDamage, GATK, VEP) are established tools. The novel contribution is the **integration and scoring framework**: a single pipeline that takes ancient DNA and a proxy genome as inputs and produces a standardised, multi-layer feasibility report. No comparable integrated tool currently exists. The framework makes feasibility assessment reproducible and comparable across different de-extinction proposals.

---

## Wider Relevance

- Methods paper candidate: the scoring framework is novel and generalisable.
- Directly relevant to the Colossal Biosciences thylacine programme, which requires exactly this kind of feasibility quantification.
- Generalises to any extinct/proxy species pair with available ancient DNA (mammoth/elephant, passenger pigeon/band-tailed pigeon, etc.).
- Supports evidence-based policy discussion around de-extinction regulation and resource allocation.
