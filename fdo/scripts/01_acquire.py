"""
FR1 · Structure Acquisition
Downloads and curates human FTO structures from RCSB PDB.

Usage:
    python scripts/01_acquire.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).parent.parent))
from scripts.utils import DATA_DIR, load_config, save_json, setup_logger

log = setup_logger("01_acquire")
STRUCTURES_DIR = DATA_DIR / "structures"


def download_pdb(pdb_id: str) -> Path:
    out = STRUCTURES_DIR / f"{pdb_id}.pdb"
    if out.exists() and out.stat().st_size > 1_000:
        log.info(f"  {pdb_id}: already cached")
        return out
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    log.info(f"  {pdb_id}: downloading from RCSB …")
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    out.write_text(r.text, encoding="utf-8")
    log.info(f"  {pdb_id}: saved ({out.stat().st_size:,} bytes)")
    return out


def fetch_rcsb_meta(pdb_id: str) -> dict:
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    try:
        d = requests.get(url, timeout=20).json()
        info = d.get("rcsb_entry_info", {})
        return {
            "title":       d.get("struct", {}).get("title", ""),
            "resolution":  (info.get("resolution_combined") or [None])[0],
            "method":      info.get("experimental_method", ""),
            "deposition":  d.get("rcsb_accession_info", {}).get("deposit_date", ""),
            "release":     d.get("rcsb_accession_info", {}).get("initial_release_date", ""),
            "n_atoms":     info.get("deposited_atom_count"),
        }
    except Exception as exc:
        log.warning(f"  {pdb_id}: metadata fetch failed — {exc}")
        return {}


def parse_pdb(pdb_file: Path) -> dict:
    chains, residues, ligands = set(), set(), set()
    with open(pdb_file, encoding="utf-8", errors="ignore") as f:
        for line in f:
            rec = line[:6].strip()
            if rec == "ATOM":
                chains.add(line[21])
                residues.add((line[21], line[22:26].strip()))
            elif rec == "HETATM":
                name = line[17:20].strip()
                if name not in ("HOH", "WAT", "DOD"):
                    ligands.add(name)
    return {
        "chains":    sorted(chains),
        "n_residues": len(residues),
        "ligands":   sorted(ligands),
    }


_FTO_KEYWORDS = {"fto", "fat mass", "obesity-associated", "ftm", "alkbh9"}


def _validate_fto(pdb_id: str, title: str) -> None:
    """Warn if the structure title does not look like an FTO entry."""
    title_lower = title.lower()
    if not any(kw in title_lower for kw in _FTO_KEYWORDS):
        log.warning(
            f"  {pdb_id}: POSSIBLE WRONG STRUCTURE — title does not mention FTO:\n"
            f"    '{title}'\n"
            f"  Verify this PDB ID maps to a human FTO structure before proceeding."
        )


def acquire_structures(config: dict) -> dict:
    STRUCTURES_DIR.mkdir(parents=True, exist_ok=True)
    pdb_ids  = config["fto"]["pdb_ids"]
    metadata = {}

    for pdb_id in pdb_ids:
        log.info(f"Processing {pdb_id} …")
        entry: dict = {"pdb_id": pdb_id}
        try:
            pdb_path = download_pdb(pdb_id)
            entry.update(fetch_rcsb_meta(pdb_id))
            entry.update(parse_pdb(pdb_path))
            entry["file_path"] = str(pdb_path.relative_to(Path(__file__).parent.parent))
            entry["status"] = "ok"
            _validate_fto(pdb_id, entry.get("title", ""))
            log.info(
                f"  {pdb_id}: {entry.get('resolution')} Å  "
                f"chains={entry.get('chains')}  "
                f"residues={entry.get('n_residues')}  "
                f"ligands={entry.get('ligands')}"
            )
        except Exception as exc:
            log.error(f"  {pdb_id}: FAILED — {exc}")
            entry["status"] = "failed"
            entry["error"]  = str(exc)
        metadata[pdb_id] = entry

    save_json(metadata, STRUCTURES_DIR / "metadata.json")
    ok = sum(1 for v in metadata.values() if v["status"] == "ok")
    log.info(f"Acquired {ok}/{len(metadata)} structures → data/structures/metadata.json")
    return metadata


if __name__ == "__main__":
    cfg = load_config()
    result = acquire_structures(cfg)
    failed = [k for k, v in result.items() if v["status"] != "ok"]
    if failed:
        log.warning(f"Failed: {failed}")
    sys.exit(0 if not failed else 1)
