"""Joint ranking of guide sets for germline-scale multiplex editing.

TODO: implement. Key logic:
- Accept a list of per-locus guide candidates (output of ancient_target_scorer).
- Score guide *sets* jointly: aggregate off-target burden, inter-guide spacing
  constraints (avoid NHEJ-mediated deletions between nearby cuts), combined
  bystander-edit collision rate across the full array.
- Return ranked list of guide sets with composite scores.

Depends on: deextinct_core.ProxyGenome, .salvaged.cfd_scorer
"""

from __future__ import annotations


def rank_multiplex_sets(
    locus_candidates: list[list[dict]],
    proxy_genome,  # deextinct_core.ProxyGenome
    max_guides: int = 10,
) -> list[dict]:
    """Stub — rank guide combinations for multiplex germline editing."""
    raise NotImplementedError
