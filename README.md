# lorg

A focused de-extinction and ancient-DNA toolkit for proxy-genome engineering.

*lorg* is Scottish Gaelic for "trace" or "track" — following the genomic trail
from ancient DNA to a living proxy, and from there to a viable restored population.

## Organising principle

Every tool must answer: *name the established package it would be cited instead of,
and the de-extinction-specific reason yours exists anyway.* If the second answer is
"none," it is a wrapper — see PORTFOLIO_TRIAGE.md for the full rationale.

The defensible niche is **degraded/ancient DNA + proxy-genome engineering** — the
exact setting where scikit-allel, PLINK, CRISPOR, BEAST, etc. silently assume clean
modern data and a single reference individual.

---

## Core (novel) — 14 tools

Problems no existing package addresses. The contribution is the *problem definition
and decision output*, even where the internal math is standard.

| Tool | Why it exists |
|------|---------------|
| [`proxy-species-edit-burden-calculator`](proxy-species-edit-burden-calculator/) | **Flagship.** Computes "edits required to move a proxy genome toward a reconstructed target." The metric itself is the publishable contribution. |
| [`large-insert-tolerance-predictor`](large-insert-tolerance-predictor/) | Position-effect risk for multi-gene cassettes — novel problem, nobody ships this. |
| [`erv-risk-mapper`](erv-risk-mapper/) | Novel biosafety framing for genome-engineering campaigns involving ERV activation risk. |
| [`regulatory-element-conservation-scorer`](regulatory-element-conservation-scorer/) | Non-coding component of edit burden; de-extinction-specific decision layer. |
| [`regulatory-rewiring-analyser`](regulatory-rewiring-analyser/) | Gained/lost/relocated enhancers between proxy and target — novel comparative output. |
| [`genetic-rescue-viability-estimator`](genetic-rescue-viability-estimator/) | Novel intervention-outcome model for managed genetic rescue. |
| [`genetic-rescue-candidate-selector`](genetic-rescue-candidate-selector/) | Donor ranking balancing outbreeding risk — novel composite score. |
| [`population-viability-genomics-estimator`](population-viability-genomics-estimator/) | Genomic MVP / extinction-probability integration. |
| [`conservation-priority-ranker`](conservation-priority-ranker/) | Multi-population composite urgency score. |
| [`inbreeding-risk-forecaster`](inbreeding-risk-forecaster/) | Forward projection of inbreeding trajectory, not just current measurement. |
| [`ancient-sample-authenticator`](ancient-sample-authenticator/) | Composite authenticity verdict tuned for de-extinction reference targets. |
| [`phylogenetic-distance-estimator`](phylogenetic-distance-estimator/) | Framed as proxy/donor ranking; wraps standard distance methods but owns the decision layer. |
| [`genome-edit-feasibility-scorer`](genome-edit-feasibility-scorer/) | Project-level decision layer integrating edit-burden inputs across loci and delivery constraints. |
| [`multiplex-edit-planner`](multiplex-edit-planner/) | Germline-scale, many-edit planning is the de-extinction reality. Schedules batched edits across embryo/cell programmes. |

---

## aDNA/proxy-aware — 14 tools

Each wraps an established package and contributes **one novel layer**: a
damage-aware null model, a single-proxy-individual mode, ancient-tip handling,
or an interpretation/threshold layer. The baseline each differs from is stated
in each tool's `README.md → Approach` section.

| Tool | Baseline | Novel layer |
|------|----------|-------------|
| [`ancestral-state-reconstructor`](ancestral-state-reconstructor/) | PastML / FastML | Target-genome reconstructor: reconstructs extinct sequence from proxy + aDNA, propagating damage uncertainty. **Second flagship.** |
| [`ancient-dna-damage-classifier`](ancient-dna-damage-classifier/) | mapDamage2 / PMDtools | Wraps damage profiling tools; adds Bayesian per-read classifier and de-extinction authentication output. Absorbs `dna-fragmentation-profiler` scope. |
| [`contamination-estimator`](contamination-estimator/) | ANGSD / schmutzi / ContamMix | Wraps all three; adds composite verdict layer. Does not re-derive likelihoods. |
| [`palaeogenomic-coverage-assessor`](palaeogenomic-coverage-assessor/) | samtools depth | aDNA-specific metrics (endogenous fraction, duplication rate, breadth at low coverage) in one report. |
| [`ancestral-gene-content-reconstructor`](ancestral-gene-content-reconstructor/) | Count / BadiRate | Input to target reconstruction; wraps gene-family evolution models. |
| [`roh-interpreter`](roh-interpreter/) | PLINK / bcftools ROH | Interpretation layer on top of ROH calls; aDNA-aware (low-coverage, pseudo-haploid). |
| [`positive-selection-signal-detector`](positive-selection-signal-detector/) | selscan / iHS | Reframed as "traits the proxy must recapitulate," with a damage-aware null. Not generic sweep-scanning. |
| [`introgression-detector`](introgression-detector/) | Dsuite / ADMIXTOOLS | Justified by *ancient* introgression (degraded D-stat null). Otherwise it is Dsuite. |
| [`effective-population-size-estimator`](effective-population-size-estimator/) | NeEstimator / PSMC | Keeps only the temporal/aDNA Ne angle; drops parts that duplicate NeEstimator. |
| [`lineage-divergence-dater`](lineage-divergence-dater/) | BEAST / OxCal | Keeps the aDNA tip-dating mode; drops strict-clock parts BEAST owns. |
| [`adaptive-diversity-scorer`](adaptive-diversity-scorer/) | scikit-allel | Composite is fine; depends on scikit-allel for Fst, owns the adaptive/neutral classification layer. |
| [`gene-loss-detector`](gene-loss-detector/) | CESAR / Ortholog Finder | Lineage-specific loss mapped to proxy; kept only if degradation-robust. |
| [`kinship-coefficient-calculator`](kinship-coefficient-calculator/) | KING / READ / lcMLkin | Novel only as aDNA kinship (low-coverage). Otherwise it is KING. |
| [`synteny-block-visualiser`](synteny-block-visualiser/) | MCScanX / JCVI | Consolidated proxy-vs-target synteny tool. Absorbs `conserved-synteny-detector` and `genome-rearrangement-mapper` (both stubs, deleted). |

---

## Shared infra — 3 items

| Item | Role |
|------|------|
| [`proxy-edit-designer/`](proxy-edit-designer/) | Consolidated CRISPR replacement for 20 deleted tools. Wraps CRISPOR/PrimeDesign; adds ancient-target scoring, multiplex ranking, and edit-burden emission. `src/salvaged/` contains production kernels (CFD scorer, off-target finder, base-editor analyser, HDR template builder) rescued from the deleted cluster. |
| [`deextinct_core/`](deextinct_core/) | Shared Python package (`pip install -e deextinct_core/`). Defines `ProxyGenome`, `TargetReconstruction`, `DamageProfile` — the object graph passed between tools. All dataclasses support JSON round-trips for provenance emission. |
| [`validation/`](validation/) | Simulation-based test harness. msprime + SLiM for popgen tools; `post-mortem-damage-simulator` (Briggs model) for aDNA tools. Per-tool accuracy benchmarks against known ground truth — non-negotiable for publication. |

---

## What was removed (and why)

See [PORTFOLIO_TRIAGE.md](PORTFOLIO_TRIAGE.md) for the full rationale.
In brief: **~70 tools deleted** across six categories:

| Category | Count | Reason |
|----------|-------|--------|
| Generic population genetics | 9 | scikit-allel / PLINK / ADMIXTURE / dadi / BEAST own this |
| Clinical / human-variant | 11 | Wrong domain — VEP / CADD / SpliceAI / InterVar |
| Generic assembly / annotation / RNA-seq | 16 | QUAST / BUSCO / AGAT / RepeatMasker / WGCNA |
| Generic comparative genomics | 9 | OrthoFinder / CAFE / MCScanX / liftOver |
| Synthetic biology | 6 | Different product — COBRApy / Benchling |
| CRISPR cluster (near-duplicates) | 20 | CRISPOR / PrimeDesign / inDelphi / BE-Hive; consolidated into `proxy-edit-designer` |

---

## Directory count

| Category | Count |
|----------|-------|
| Core (novel) tools | 14 |
| aDNA/proxy-aware (PIVOT) tools | 14 |
| Shared infra additions | 3 |
| **Total shipped** | **31** |
| Internal benchmark utility (`validation/post-mortem-damage-simulator`) | 1 |
