# deextinct_core

Shared Python package defining the aDNA-aware data model passed between all
tools in this repository. The integration that no other package offers: a single
object graph that carries proxy-genome identity, target reconstruction state, and
damage uncertainty together, so each tool receives a consistent view of the project.

## Data model

### `ProxyGenome`
Represents the living proxy species genome being engineered.

Fields:
- `species` — proxy species name (e.g. `"Loxodonta africana"`)
- `assembly_id` — reference assembly identifier
- `fasta_path` — path to proxy FASTA
- `annotation_path` — optional GFF3/GTF annotation
- `metadata` — dict of arbitrary provenance fields

### `TargetReconstruction`
The reconstructed extinct-target sequence at one or more loci.

Fields:
- `target_species` — extinct species name (e.g. `"Mammuthus primigenius"`)
- `loci` — list of `ReconstructedLocus` objects (chrom, start, end, sequence,
  per-site posterior probabilities)
- `source_samples` — list of aDNA sample identifiers used in reconstruction
- `damage_profile` — associated `DamageProfile`
- `method` — reconstruction method (e.g. `"ancestral-state-reconstructor v0.1"`)

### `DamageProfile`
Encodes post-mortem damage characteristics for one aDNA library.

Fields:
- `sample_id`
- `ct_rate_5prime` — list of C→T rates per position from 5' end
- `ga_rate_3prime` — list of G→A rates per position from 3' end
- `mean_fragment_length` — float
- `authenticity_posterior` — float in [0, 1]
- `model` — decay model used (`"geometric"`, `"briggs"`, …)

## Usage pattern

```python
from deextinct_core import ProxyGenome, TargetReconstruction, DamageProfile

# Typically constructed by ancestral-state-reconstructor and ancient-dna-damage-classifier;
# consumed by proxy-species-edit-burden-calculator, proxy-edit-designer, and PIVOT tools.
proxy = ProxyGenome(species="Loxodonta africana", assembly_id="loxAfr3", ...)
damage = DamageProfile(sample_id="ULM001", ct_rate_5prime=[0.32, 0.12, ...], ...)
target = TargetReconstruction(target_species="Mammuthus primigenius", ..., damage_profile=damage)
```

## Serialisation

All three dataclasses support `.to_dict()` / `.from_dict()` for JSON round-trips,
enabling provenance/manifest emission (see PORTFOLIO_TRIAGE.md ADD item 4).
