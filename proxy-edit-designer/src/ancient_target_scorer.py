"""Re-score CRISPOR off-target hits against a TargetReconstruction.

TODO: implement. Key logic:
- Accept TargetReconstruction from deextinct_core with per-site posteriors.
- For each CRISPOR off-target hit overlapping a reconstructed locus, apply a
  penalty proportional to the site's damage uncertainty (low posterior ->
  higher effective off-target risk, because the reference may not reflect the
  true ancient sequence).
- Return augmented hit list with `ancient_adjusted_cfd` field.

Depends on: deextinct_core.TargetReconstruction, .salvaged.cfd_scorer
"""

from __future__ import annotations


def rescore_hits_against_reconstruction(
    hits: list[dict],
    reconstruction,  # deextinct_core.TargetReconstruction
    uncertainty_penalty_scale: float = 2.0,
) -> list[dict]:
    """Stub — apply ancient-reconstruction uncertainty penalty to CFD scores."""
    raise NotImplementedError
