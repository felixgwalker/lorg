"""
Download and index BLAST databases required by hgt-risk-assessor.

Usage
-----
    python data/download_databases.py --all
    python data/download_databases.py --isfinder --integrall
    python data/download_databases.py --conjugative --entrez-email you@example.com
    python data/download_databases.py --phage

Requirements
------------
    pip install biopython requests          (always)
    BLAST+ makeblastdb on PATH              (or use --makeblastdb /path/to/makeblastdb)
"""

import argparse
import gzip
import json
import logging
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import requests

logging.basicConfig(level=logging.INFO, format="%(levelname)s  %(message)s")
logger = logging.getLogger(__name__)

DATA_DIR = Path(__file__).parent

# ---------------------------------------------------------------------------
# Database registry
# ---------------------------------------------------------------------------

DATABASES = {
    "isfinder": {
        "description":  "ISfinder insertion sequence database",
        "url":          "https://isfinder.biotoul.fr/blast/isfinder.fasta",
        "filename":     "isfinder.fasta",
        "dbtype":       "nucl",
        "subdir":       "isfinder",
        "db_name":      "isfinder",
        "fallback": (
            "Visit https://isfinder.biotoul.fr/blast.php and download the FASTA "
            "manually.  Save to data/isfinder/isfinder.fasta, then:\n"
            "  makeblastdb -in data/isfinder/isfinder.fasta -dbtype nucl "
            "-out data/isfinder/isfinder -title ISfinder"
        ),
    },
    "integrall": {
        "description":  "INTEGRALL integron sequence database",
        "url":          "http://integrall.bio.ua.pt/download/integrall_sequences.fasta",
        "filename":     "integrall.fasta",
        "dbtype":       "nucl",
        "subdir":       "integrall",
        "db_name":      "integrall",
        "fallback": (
            "Visit http://integrall.bio.ua.pt/download and download the sequences "
            "FASTA manually.  Save to data/integrall/integrall.fasta, then:\n"
            "  makeblastdb -in data/integrall/integrall.fasta -dbtype nucl "
            "-out data/integrall/integrall -title INTEGRALL"
        ),
    },
    "conjugative": {
        "description":  "Conjugative element protein database (MOB relaxases, T4SS)",
        "url":          None,   # built via NCBI Entrez
        "filename":     "conjugative_proteins.fasta",
        "dbtype":       "prot",
        "subdir":       "conjugative",
        "db_name":      "conjugative_proteins",
        "fallback": (
            "Manually collect MOB-family relaxase and T4SS protein sequences in "
            "FASTA format, save to data/conjugative/conjugative_proteins.fasta, then:\n"
            "  makeblastdb -in data/conjugative/conjugative_proteins.fasta "
            "-dbtype prot -out data/conjugative/conjugative_proteins "
            "-title ConjugativeProteins"
        ),
    },
    "phage": {
        "description":  "NCBI RefSeq viral genomes (prophage fallback database)",
        "url": (
            "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/"
            "viral.1.1.genomic.fna.gz"
        ),
        "filename":     "phage_genes.fasta.gz",
        "dbtype":       "nucl",
        "subdir":       "phage",
        "db_name":      "phage_genes",
        "fallback": (
            "Download from https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/ "
            "(viral.1.1.genomic.fna.gz), decompress, save to data/phage/phage_genes.fasta, then:\n"
            "  makeblastdb -in data/phage/phage_genes.fasta -dbtype nucl "
            "-out data/phage/phage_genes -title PhageGenes"
        ),
    },
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _download(url: str, dest: Path, chunk: int = 1024 * 1024) -> None:
    logger.info(f"Downloading {url}")
    with requests.get(url, stream=True, timeout=60) as resp:
        resp.raise_for_status()
        total = int(resp.headers.get("content-length", 0))
        done = 0
        with open(dest, "wb") as fh:
            for buf in resp.iter_content(chunk):
                fh.write(buf)
                done += len(buf)
                if total:
                    print(f"\r  {done/1e6:.1f}/{total/1e6:.1f} MB  "
                          f"({done/total*100:.0f}%)", end="", flush=True)
    print()
    logger.info(f"Saved to {dest}")


def _decompress_gz(gz_path: Path) -> Path:
    out = gz_path.with_suffix("")
    logger.info(f"Decompressing {gz_path.name}")
    with gzip.open(gz_path, "rb") as fi, open(out, "wb") as fo:
        while chunk := fi.read(1024 * 1024):
            fo.write(chunk)
    gz_path.unlink()
    logger.info(f"Decompressed → {out.name}")
    return out


def _makeblastdb(fasta: Path, db: Path, dbtype: str, title: str, binary: str) -> bool:
    cmd = [binary, "-in", fasta.as_posix(), "-dbtype", dbtype,
           "-out", db.as_posix(), "-title", title, "-parse_seqids"]
    logger.info(f"Running makeblastdb for {title}")
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        logger.error(f"makeblastdb failed:\n{r.stderr}")
        return False
    logger.info("makeblastdb completed.")
    return True


def _write_manifest(db_dir: Path, name: str, url: Optional[str],
                    fasta: Path, db: Path) -> None:
    try:
        n = sum(1 for ln in open(fasta) if ln.startswith(">"))
    except Exception:
        n = -1
    manifest = {
        "name": name,
        "source": url or "ncbi_entrez",
        "built": datetime.now(timezone.utc).isoformat(),
        "fasta": str(fasta),
        "blast_db": str(db),
        "record_count": n,
    }
    path = db_dir / f"{name}.manifest.json"
    path.write_text(json.dumps(manifest, indent=2))
    logger.info(f"Manifest: {path}")


def _db_ready(db_dir: Path, info: dict) -> bool:
    db = db_dir / info["db_name"]
    ext = ".nhr" if info["dbtype"] == "nucl" else ".phr"
    return db.with_suffix(ext).exists()


# ---------------------------------------------------------------------------
# Per-database build functions
# ---------------------------------------------------------------------------

def _build_standard(key: str, data_dir: Path, makeblastdb: str, force: bool) -> bool:
    info = DATABASES[key]
    db_dir = data_dir / info["subdir"]
    db_dir.mkdir(parents=True, exist_ok=True)

    if _db_ready(db_dir, info) and not force:
        logger.info(f"{key}: already built. Use --force to rebuild.")
        return True

    raw = db_dir / info["filename"]
    try:
        _download(info["url"], raw)
    except Exception as exc:
        logger.error(f"Download failed: {exc}")
        logger.info(f"\nFallback:\n{info['fallback']}\n")
        return False

    fasta = _decompress_gz(raw) if info["filename"].endswith(".gz") else raw
    db = db_dir / info["db_name"]
    ok = _makeblastdb(fasta, db, info["dbtype"], key.title(), makeblastdb)
    if ok:
        _write_manifest(db_dir, key, info["url"], fasta, db)
    return ok


def _build_conjugative(data_dir: Path, entrez_email: str,
                        makeblastdb: str, force: bool) -> bool:
    info = DATABASES["conjugative"]
    db_dir = data_dir / info["subdir"]
    db_dir.mkdir(parents=True, exist_ok=True)

    if _db_ready(db_dir, info) and not force:
        logger.info("conjugative: already built. Use --force to rebuild.")
        return True

    if not entrez_email:
        logger.error("--entrez-email is required to build the conjugative database.")
        logger.info(f"\nFallback:\n{info['fallback']}\n")
        return False

    try:
        from Bio import Entrez, SeqIO
    except ImportError:
        logger.error("biopython is required: pip install biopython")
        return False

    Entrez.email = entrez_email

    queries = [
        "relaxase[Title] AND (plasmid[Title] OR conjugat*[Title]) AND bacteria[Organism]",
        "VirB4[Protein Name] AND bacteria[Organism]",
        "VirD4[Protein Name] AND bacteria[Organism]",
        "MOB[Title] AND relaxase AND bacteria[Organism]",
        "TraI[Title] AND conjugat* AND bacteria[Organism]",
    ]

    fasta_path = db_dir / info["filename"]
    all_records = []

    for query in queries:
        try:
            sh = Entrez.esearch(db="protein", term=query, retmax=500)
            ids = Entrez.read(sh)["IdList"]
            sh.close()
            if not ids:
                continue
            logger.info(f"  {len(ids)} hits: {query[:60]}…")
            fh = Entrez.efetch(
                db="protein", id=",".join(ids[:500]),
                rettype="fasta", retmode="text"
            )
            all_records.extend(SeqIO.parse(fh, "fasta"))
            fh.close()
        except Exception as exc:
            logger.warning(f"  Query failed: {exc}")

    if not all_records:
        logger.error("No sequences retrieved — check NCBI connectivity.")
        logger.info(f"\nFallback:\n{info['fallback']}\n")
        return False

    seen: set[str] = set()
    unique = [r for r in all_records if r.id not in seen and not seen.add(r.id)]  # type: ignore[func-returns-value]
    logger.info(f"Writing {len(unique)} unique protein sequences to {fasta_path.name}")
    with open(fasta_path, "w") as fh:
        SeqIO.write(unique, fh, "fasta")

    db = db_dir / info["db_name"]
    ok = _makeblastdb(fasta_path, db, "prot", "ConjugativeProteins", makeblastdb)
    if ok:
        _write_manifest(db_dir, "conjugative", None, fasta_path, db)
    return ok


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    p = argparse.ArgumentParser(
        description="Download and index BLAST databases for hgt-risk-assessor"
    )
    p.add_argument("--all",         action="store_true", help="Build all databases")
    p.add_argument("--isfinder",    action="store_true")
    p.add_argument("--integrall",   action="store_true")
    p.add_argument("--conjugative", action="store_true",
                   help="Build conjugative protein DB (requires --entrez-email)")
    p.add_argument("--phage",       action="store_true",
                   help="Download NCBI RefSeq viral genomes (phage fallback)")
    p.add_argument("--data-dir",    type=Path, default=DATA_DIR)
    p.add_argument("--makeblastdb", default="makeblastdb",
                   help="Path to makeblastdb binary")
    p.add_argument("--entrez-email", metavar="EMAIL",
                   help="NCBI Entrez email (required for --conjugative)")
    p.add_argument("--force", action="store_true",
                   help="Rebuild even if databases already exist")
    args = p.parse_args()

    targets = []
    if args.all:
        targets = list(DATABASES)
    else:
        for key in ("isfinder", "integrall", "conjugative", "phage"):
            if getattr(args, key):
                targets.append(key)

    if not targets:
        p.print_help()
        sys.exit(0)

    results: dict[str, bool] = {}
    for key in targets:
        logger.info(f"{'='*60}")
        logger.info(f"Building: {DATABASES[key]['description']}")
        if key == "conjugative":
            results[key] = _build_conjugative(
                args.data_dir, args.entrez_email or "",
                args.makeblastdb, args.force
            )
        else:
            results[key] = _build_standard(key, args.data_dir, args.makeblastdb, args.force)

    print("\n--- Summary ---")
    for key, ok in results.items():
        print(f"  {key:<20}  {'OK' if ok else 'FAILED'}")

    if not all(results.values()):
        sys.exit(1)


if __name__ == "__main__":
    main()
