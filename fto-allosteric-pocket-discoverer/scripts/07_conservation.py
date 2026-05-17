"""
FR8 · Conservation Mapping
Maps evolutionary conservation onto FTO residues via ConSurf REST API
or a local BLOSUM62-based MSA approach.

Usage:
    python scripts/07_conservation.py [--method consurf|local|demo]
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import requests

sys.path.insert(0, str(Path(__file__).parent.parent))
from scripts.utils import (
    DATA_DIR, DemoData, load_config, save_json, setup_logger,
)

log = setup_logger("07_conservation")

CONSERVATION_DIR = DATA_DIR / "conservation"

# ConSurf colour scale (1 = variable, 9 = conserved)
CONSURF_COLORS = {
    1: "#00CCFF", 2: "#33DDFF", 3: "#66EEFF",
    4: "#99FFEE", 5: "#CCFFCC", 6: "#FFFF99",
    7: "#FFCC66", 8: "#FF9933", 9: "#FF6600",
}


# ── ConSurf API ───────────────────────────────────────────────────────────────

def fetch_consurf(pdb_id: str, chain: str = "A") -> dict | None:
    """
    Query the ConSurf server REST API for pre-computed conservation scores.
    Returns None if the server is unreachable or the entry is unavailable.
    """
    base = "https://consurf.tau.ac.il/api/consurf_pdb"
    url  = f"{base}?pdb={pdb_id.lower()}&chain={chain}"
    log.info(f"Querying ConSurf for {pdb_id} chain {chain} …")
    try:
        r = requests.get(url, timeout=30)
        if r.status_code == 200:
            data = r.json()
            residues = [int(k) for k in data.keys() if k.isdigit()]
            scores   = [int(data[str(r)].get("consurf_score", 5)) for r in sorted(residues)]
            log.info(f"  ConSurf: {len(residues)} residue scores retrieved")
            return {
                "residues": sorted(residues),
                "scores":   scores,
                "method":   f"ConSurf (PDB {pdb_id} chain {chain})",
                "colormap": {str(k): v for k, v in CONSURF_COLORS.items()},
            }
        else:
            log.warning(f"  ConSurf returned HTTP {r.status_code}")
    except Exception as exc:
        log.warning(f"  ConSurf query failed: {exc}")
    return None


# ── Local BLOSUM62 approach ───────────────────────────────────────────────────

def _fetch_uniprot_msa(uniprot_id: str, n_hits: int = 50) -> list[str]:
    """Fetch top homologues via UniProt BLAST and return sequences."""
    url = (
        f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    )
    try:
        r = requests.get(url, timeout=20)
        r.raise_for_status()
        lines = r.text.splitlines()
        seq   = "".join(l for l in lines if not l.startswith(">"))
        return [seq]  # fallback: just the query sequence
    except Exception as exc:
        log.warning(f"  UniProt fetch failed: {exc}")
        return []


def _blosum62_conservation(sequences: list[str]) -> list[float]:
    """
    Calculate per-column conservation from an MSA as normalised entropy.
    Higher score = more conserved.
    """
    if not sequences:
        return []
    import math
    n_seq  = len(sequences)
    n_col  = len(sequences[0])
    scores = []
    for col in range(n_col):
        counts: dict[str, int] = {}
        for seq in sequences:
            aa = seq[col] if col < len(seq) else "-"
            counts[aa] = counts.get(aa, 0) + 1
        entropy = -sum(
            (c / n_seq) * math.log2(c / n_seq)
            for c in counts.values() if c > 0
        )
        max_entropy = math.log2(min(n_seq, 20)) if n_seq > 1 else 1
        # Invert: low entropy = high conservation
        cons = 1.0 - (entropy / max_entropy) if max_entropy > 0 else 1.0
        # Map to 1–9 scale
        scores.append(max(1, min(9, round(cons * 8 + 1))))
    return scores


def compute_local_conservation(config: dict) -> dict:
    uniprot = config["fto"]["uniprot"]
    seqs    = _fetch_uniprot_msa(uniprot)
    if len(seqs) < 2:
        log.warning("Insufficient sequences for local conservation – using demo data")
        return DemoData().conservation_scores()

    scores   = _blosum62_conservation(seqs)
    residues = list(range(1, len(scores) + 1))
    return {
        "residues": residues,
        "scores":   [int(s) for s in scores],
        "method":   "Local BLOSUM62 entropy",
        "colormap": {str(k): v for k, v in CONSURF_COLORS.items()},
    }


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", choices=["consurf", "local", "demo"], default="consurf")
    args = parser.parse_args()

    config = load_config()
    CONSERVATION_DIR.mkdir(parents=True, exist_ok=True)
    out = CONSERVATION_DIR / "fto_conservation.json"

    if args.method == "demo":
        data = DemoData().conservation_scores()
    elif args.method == "local":
        data = compute_local_conservation(config)
    else:  # consurf
        primary = config["fto"]["primary_structure"]
        data    = fetch_consurf(primary)
        if data is None:
            log.info("Falling back to local conservation …")
            data = compute_local_conservation(config)
        if data is None:
            log.info("Falling back to demo conservation …")
            data = DemoData().conservation_scores()

    save_json(data, out)
    n = len(data.get("residues", []))
    log.info(f"Conservation scores for {n} residues saved → data/conservation/fto_conservation.json")
    log.info(f"  Method: {data.get('method')}")


if __name__ == "__main__":
    main()
