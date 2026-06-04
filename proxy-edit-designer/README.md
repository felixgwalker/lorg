# Proxy Edit Designer

Consolidated CRISPR edit design tool for de-extinction genome-engineering campaigns.
Wraps CRISPOR and PrimeDesign, then adds the layer those tools omit: scoring every
candidate edit against a *reconstructed ancient target sequence* (produced by
`ancestral-state-reconstructor`) rather than a modern reference.

## Why this tool exists

**Baseline:** CRISPOR (Concordet & Haeussler 2018) and PrimeDesign (Hsu et al. 2021)
are best-in-class guide design tools. This tool does not re-implement their scoring
logic — it calls them as dependencies.

**The gap they leave:** Both tools assume (a) a clean modern reference genome and
(b) a single target sequence. De-extinction projects face neither: the target is a
*reconstructed* ancient sequence with per-site damage uncertainty, and the engineering
campaign requires *dozens to hundreds* of simultaneous edits across a proxy genome
(germline multiplex scale). Neither tool surfaces aggregate edit burden, off-target
enrichment near ancient indels, or bystander-edit risk in multi-guide arrays targeting
reconstructed loci.

**This tool's novel layer:**
1. Accept a `TargetReconstruction` (from `deextinct_core`) with per-site posterior
   uncertainty, and propagate that uncertainty into the guide specificity score.
2. Re-score CRISPOR off-target hits that overlap ancient indels or low-confidence
   reconstructed positions (these are systematically under-penalised by standard CFD
   scoring against a modern proxy).
3. Support multiplex array design: rank guide sets jointly by aggregate off-target
   burden, PAM availability, and bystander-edit collision across the full edit plan.
4. Emit edit specs to `proxy-species-edit-burden-calculator` for project-level
   feasibility scoring.

## Approach

```
Input: proxy genome FASTA + TargetReconstruction (from ancestral-state-reconstructor)
       + list of desired coding/regulatory changes

Step 1  CRISPOR / PrimeDesign (subprocess call) — generate candidate guides and
        pegRNA designs for each locus.

Step 2  Off-target re-scoring against proxy genome using CFD scoring
        (salvaged/cfd_scorer.py; Doench et al. 2016) with ancient-locus penalty
        for sites overlapping reconstructed positions with posterior < threshold.

Step 3  Bystander-edit analysis per guide using base-editor window mapping
        (salvaged/base_editor_analyser.py + salvaged/base_editor_config.py).

Step 4  PAM disruption check for HDR templates
        (salvaged/pam_disruptor.py + salvaged/edit_applicator.py).

Step 5  Multiplex ranking: score guide *sets* by total off-target burden,
        inter-guide spacing, and combined bystander collision rate.

Step 6  Output ranked guide set + pegRNA designs + edit specs for burden calculator.
```

## Salvaged modules

The `src/salvaged/` directory contains production-quality implementations
rescued from the 20 deleted CRISPR tools (see stage1f commit):

| Module | Origin | Purpose |
|---|---|---|
| `cfd_scorer.py` | guide-rna-off-target-scorer | Doench 2016 CFD off-target scoring |
| `off_target_finder.py` | guide-rna-off-target-scorer | Seed-and-extend genome search |
| `genome_indexer.py` | guide-rna-off-target-scorer | 20-mer k-mer genome index |
| `guide_parser.py` | guide-rna-off-target-scorer | FASTA/CSV/raw guide input parsing |
| `edit_applicator.py` | hdr-template-designer | SNP/insertion/deletion template builder |
| `pam_disruptor.py` | hdr-template-designer | Silent PAM-breaking mutation suggester |
| `reference_reader.py` | hdr-template-designer | FASTA locus extractor |
| `base_editor_analyser.py` | crispr-base-editor-window-visualiser | Per-position base-editor window analysis |
| `base_editor_config.py` | crispr-base-editor-window-visualiser | Published editor profiles (ABE7.10, BE4max, …) |

## Stub modules (to implement)

- `src/ancient_target_scorer.py` — re-score CRISPOR hits against `TargetReconstruction`
- `src/multiplex_planner.py` — joint ranking of guide sets at germline scale
- `src/crispor_runner.py` — subprocess wrapper for CRISPOR
- `src/primedesign_runner.py` — subprocess wrapper for PrimeDesign
- `src/burden_emitter.py` — format edit specs for `proxy-species-edit-burden-calculator`

## Dependencies

- CRISPOR (external, called via subprocess)
- PrimeDesign (external, called via subprocess)
- `deextinct_core` (this repo) — `ProxyGenome`, `TargetReconstruction`, `DamageProfile`
- numpy, pandas, biopython, matplotlib

## Status

Stub + salvaged kernels. Step 1 (CRISPOR/PrimeDesign subprocess wrappers) and
Steps 2–4 (salvaged scoring logic) are the starting point for implementation.
