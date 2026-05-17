# HGT Risk Assessor

## Project Overview

A sequence-level horizontal gene transfer (HGT) risk scoring tool for engineered organisms. Given an annotated genome or sequence of a synthetic construct and its intended host, the tool produces a standardised biosafety risk report before any wet-lab work begins. The aim is to provide the synthetic biology community with a reproducible, evidence-based HGT risk assessment that can accompany regulatory submissions and internal safety reviews.

---

## Scoring Framework

The tool exposes **two complementary models** that run in parallel on every invocation.

### 1 — Flat signal model (v0.1, preserved)

Five BLAST-derived signals are individually scored and combined into a single aggregate risk index using equal weighting. The index maps to four risk levels:

| Index range | Level |
|---|---|
| < 0.25 | Low |
| 0.25 – 0.50 | Medium |
| 0.50 – 0.75 | High |
| ≥ 0.75 | Critical |

### 2 — Three-layer HGT Risk Index (v0.2)

A richer, interpretable model that groups 14 named features into three mechanistic layers. Each layer gets a sub-score; the sub-scores are combined using a configurable weight profile to produce the **HGT Risk Index** (0.0–1.0).

```
HGT Risk Index = w_T × Transfer_score
               + w_E × Establishment_score
               + w_C × Consequence_score
```

#### Layers and features

| Layer | Feature | Source |
|---|---|---|
| **Transfer Opportunity** | IS element match | ISfinder BLAST |
| | Integron association | INTEGRALL BLAST |
| | Conjugative element homology | NCBI conjugative DB BLASTX |
| | Plasmid context | *(pending integration)* |
| | Transposase proximity | *(pending integration)* |
| | Repeat density | Computed from sequence |
| **Establishment** | GC deviation | Computed vs host GC |
| | Codon usage distance | Total variation distance, hardcoded E. coli / B. subtilis tables |
| | Taxonomic distance | Curated offline lineage lookup |
| | Promoter plausibility | σ-70 −35/−10 box scan |
| | Sequence complexity | Computed from sequence |
| **Functional Consequence** | Prophage context | PHASTER API |
| | AMR content | *(pending CARD integration)* |
| | Virulence flags | *(pending VFDB integration)* |
| | Gene completeness | ORF scan (6 frames) |
| | Payload count | ORF scan (6 frames) |

Features marked *pending* always return `available=False` and are excluded from the index until the relevant database is integrated. The weights of available features are re-normalised to sum to 1.0 within each layer so unavailable features never inflate or deflate the result.

#### Score bands

| Index | Band |
|---|---|
| < 0.25 | **low** — Minimal sequence-level HGT indicators. Standard biosafety procedures apply. |
| 0.25 – 0.50 | **moderate** — Some indicators present. Expert review recommended before scale-up or release. |
| 0.50 – 0.75 | **high** — Multiple significant indicators. Formal contained use risk assessment required. |
| ≥ 0.75 | **very_high** — Strong HGT risk signals. Do not proceed without biosafety officer review. |

#### Weight profiles

The relative importance of the three layers can be tuned for the biosafety scenario in question.

| Profile | Transfer | Establishment | Consequence | Best for |
|---|---|---|---|---|
| `default` | 0.40 | 0.35 | 0.25 | General biosafety review |
| `environmental` | 0.45 | 0.40 | 0.15 | Open-environment or field release |
| `clinical_amr` | 0.30 | 0.25 | 0.45 | Clinical AMR / antibiotic resistance payload |

Select a profile with `--weight-profile` on the command line.

#### Missing data behaviour

- Any feature whose signal module is skipped or whose dependencies are unavailable reports `available=False`.
- Its configured weight is removed from the denominator, and the remaining feature weights are re-normalised.
- A **completeness percentage** is reported alongside the index so that users understand how much of the model contributed to the score.
- Features that are both important and unavailable (AMR content, virulence flags, plasmid context, transposase proximity, conjugative element, IS element) are listed explicitly in the report as *missing important features*.

---

## Signals Integrated

The pipeline integrates five sequence-level signals associated with elevated HGT potential (flat model; all five also contribute to the three-layer model):

1. **Insertion Sequence (IS) Proximity** — proximity of the construct or insertion site to known IS elements, which are primary vectors of mobilisation. Scored against the ISfinder database.
2. **GC Content Deviation** — significant deviation from host genomic GC content is a hallmark of recent horizontal acquisition and may indicate ongoing mobility or instability in the host background.
3. **Integron Association** — presence of integron-associated gene cassettes or attI/attC recombination sites, which facilitate capture and dissemination of gene cassettes across species. Scored against the INTEGRALL database.
4. **Conjugative Element Homology** — homology to known conjugative plasmid or integrative conjugative element (ICE) components, indicating potential for self-transmissibility or co-mobilisation.
5. **Prophage Context** — proximity to or integration within prophage regions, which can package and transduce flanking genomic material during induction events.

---

## Usage

```
python -m src.cli \
  --input sequence.fasta \
  --host "Escherichia coli K-12" \
  --host-gc 0.509 \
  --weight-profile clinical_amr \
  --donor-taxon "Klebsiella pneumoniae" \
  --output-dir results/
```

Key options:

| Flag | Description |
|---|---|
| `--input FILE` | Input sequence file (FASTA or VCF) |
| `--host HOST_ID` | Host organism (NCBI accession or common name) |
| `--host-gc FLOAT` | Host GC content 0.0–1.0 (skips Entrez lookup) |
| `--input-format {fasta,vcf}` | Input file format (default: fasta) |
| `--ref FASTA` | Reference FASTA (required with `--input-format vcf`) |
| `--weight-profile PROFILE` | `default` / `environmental` / `clinical_amr` |
| `--donor-taxon TAXON` | Donor organism for taxonomic distance feature |
| `--recipient-taxon TAXON` | Recipient organism (defaults to `--host`) |
| `--no-network` | Disable PHASTER API and Entrez lookups (offline mode) |
| `--entrez-email EMAIL` | NCBI Entrez email (required unless `--host-gc` or `--no-network`) |
| `--skip-signal SIGNAL` | Skip one or more signals: `gc_content is_proximity integron conjugative prophage` |
| `--data-dir PATH` | Directory containing BLAST databases (default: `data/`) |
| `--output-dir PATH` | Directory for output reports (default: `results/`) |
| `--threads N` | BLAST thread count (default: 4) |

### Example console output

```
  ── Flat signal model ──────────────────────
  Risk Index : 0.412  (Medium)

  ── Three-layer model (clinical_amr) ───────
  HGT Risk Index : 0.538  (high)
  Transfer       : 0.320
  Establishment  : 0.610
  Consequence    : 0.650
  Completeness   : 64%

  Moderate-to-high establishment signals detected (GC deviation: 0.61).
  Codon usage distance is elevated, suggesting the insert was recently
  acquired from a distant donor. No AMR or virulence content could be
  confirmed (database integration pending); the consequence score may be
  underestimated. A formal contained use risk assessment is recommended.
```

---

## Data Sources

| Source | Use |
|---|---|
| ISfinder | Insertion sequence element reference database |
| INTEGRALL | Integron and gene cassette reference database |
| NCBI Prokaryotic Genome Database | Host genome context, comparative GC content, conjugative element references |
| PHASTER API | Prophage region identification |
| CARD *(pending)* | Antimicrobial resistance gene database |
| VFDB *(pending)* | Virulence factor database |

Database download helper: `python data/download_databases.py --data-dir data/`

---

## Technical Stack

| Component | Role |
|---|---|
| Python 3.10+ pipeline | Sequence analysis, signal scoring, report generation |
| BLAST+ subprocess wrapper | blastn / blastx against local databases |
| VCF or FASTA input | Flexible input for variant-level or full-sequence analysis |
| Three-layer scoring engine | `src/scoring/features.py`, `src/scoring/layers.py` |
| Rule-based explanation engine | `src/scoring/explanation.py` |
| HTML + JSON report output | Self-contained HTML (inline CSS), machine-readable JSON sidecar |
| Central configuration | `src/config.py` — all thresholds and weight profiles in one place |

---

## Project Structure

```
hgt-risk-assessor/
├── src/
│   ├── config.py           # Score bands, weight profiles, feature weights
│   ├── models.py           # Dataclasses: PipelineResult, ThreeLayerResult, …
│   ├── pipeline.py         # Orchestrates signals → aggregation → three-layer
│   ├── aggregator.py       # Flat signal aggregation (v0.1)
│   ├── blast.py            # BLAST+ subprocess wrapper
│   ├── input_parser.py     # FASTA / VCF ingestion
│   ├── report.py           # HTML + JSON report rendering
│   ├── cli.py              # argparse entry point
│   ├── signals/
│   │   ├── gc_content.py
│   │   ├── is_proximity.py
│   │   ├── integron.py
│   │   ├── conjugative.py
│   │   └── prophage.py
│   └── scoring/
│       ├── features.py     # 14 feature extractors (three-layer model)
│       ├── layers.py       # Layer aggregation, classify_band()
│       └── explanation.py  # Rule-based NL explanation generator
├── tests/
│   └── test_scoring.py     # 50-test suite for three-layer model
├── data/
│   └── download_databases.py
└── README.md
```

---

## Novel Contribution

HGT risk considerations exist in biosafety literature, but there is no standardised, accessible, sequence-level tool that consolidates them into a single automated assessment. Biosafety evaluations currently rely on expert judgement applied inconsistently across institutions. This is the **first standardised accessible HGT risk framework for synthetic biology**, designed to make risk assessment reproducible, auditable, and available to teams without specialist biosafety expertise.

---

## Wider Relevance

- **UK and Scottish GMO regulatory context**: the tool is designed to produce output aligned with APHA/HSE contained use risk assessment requirements, supporting regulatory submissions.
- **iGEM responsible science emphasis**: directly addresses iGEM's biosafety and responsible science requirements; strong project narrative for the competition.
- **Academic publication candidate**: the framework and validation methodology constitute a novel methods contribution to the biosafety literature.
- **Institutional development**: best developed with support from biosafety officers and regulatory specialists at the University of Glasgow; institutional involvement strengthens both the tool and the publication case.
- Generalises beyond synthetic biology to natural HGT surveillance in clinical and environmental microbiology.
