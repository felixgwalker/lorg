"""
FTO Pocket Discovery Pipeline Orchestrator
Runs all 10 analysis stages in sequence.

Usage:
    python run_pipeline.py                # full pipeline
    python run_pipeline.py --demo         # demo data only (no MD required)
    python run_pipeline.py --from 04      # resume from a stage
    python run_pipeline.py --only 01,02   # run specific stages
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

from rich.console import Console
from rich.panel   import Panel
from rich.table   import Table
from rich         import box

ROOT    = Path(__file__).parent
console = Console()

STAGES = {
    "01": ("01_acquire.py",       "Structure Acquisition",     False),
    "02": ("02_prepare.py",       "Structure Preparation",     False),
    "03": ("03_simulate.py",      "Molecular Dynamics",        False),
    "04": ("04_motion_analysis.py","Motion Analysis",          True),
    "05": ("05_pocket_detection.py","Pocket Detection",        True),
    "06": ("06_pocket_ranking.py", "Pocket Ranking",           True),
    "07": ("07_conservation.py",   "Conservation Mapping",     True),
    "08": ("08_event_detection.py","Event Detection",          True),
    "09": ("09_residue_network.py","Residue Network",          True),
    "10": ("10_ligand.py",         "Ligand / Pharmacophore",   True),
}


def run_stage(key: str, demo: bool = False) -> tuple[bool, float]:
    script, name, supports_demo = STAGES[key]
    script_path = ROOT / "scripts" / script

    cmd = [sys.executable, str(script_path)]
    if demo and supports_demo:
        cmd.append("--demo")

    console.rule(f"[bold cyan]Stage {key}: {name}")
    t0 = time.perf_counter()
    result = subprocess.run(cmd, text=True)
    elapsed = time.perf_counter() - t0

    ok = result.returncode == 0
    status = "[green]✓ OK" if ok else "[red]✗ FAILED"
    console.print(f"  {status}  ({elapsed:.1f}s)")
    return ok, elapsed


def print_summary(results: dict[str, tuple[bool, float]]) -> None:
    table = Table(
        title="Pipeline Summary", box=box.ROUNDED,
        style="dim", header_style="bold cyan",
    )
    table.add_column("Stage", style="bold")
    table.add_column("Name")
    table.add_column("Status")
    table.add_column("Time (s)", justify="right")

    for key, (ok, elapsed) in results.items():
        _, name, _ = STAGES[key]
        status = "[green]✓ OK" if ok else "[red]✗ FAILED"
        table.add_row(key, name, status, f"{elapsed:.1f}")

    console.print()
    console.print(table)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--demo",  action="store_true",
                        help="Use synthetic demo data for analysis stages (no MD required)")
    parser.add_argument("--from",  dest="from_stage", default=None,
                        help="Resume from stage (e.g. --from 04)")
    parser.add_argument("--only",  default=None,
                        help="Comma-separated list of stages to run (e.g. --only 01,04,06)")
    parser.add_argument("--skip",  default=None,
                        help="Comma-separated list of stages to skip")
    args = parser.parse_args()

    console.print(Panel.fit(
        "[bold white]FTO Interdomain Pocket Discovery Platform[/bold white]\n"
        "[dim]Pipeline Orchestrator v1.0[/dim]",
        border_style="cyan",
    ))

    # Determine which stages to run
    all_keys = sorted(STAGES.keys())
    if args.only:
        keys = [k.zfill(2) for k in args.only.split(",")]
    elif args.from_stage:
        start = args.from_stage.zfill(2)
        keys  = [k for k in all_keys if k >= start]
    else:
        keys = all_keys

    if args.skip:
        skip = {k.zfill(2) for k in args.skip.split(",")}
        keys = [k for k in keys if k not in skip]

    if args.demo:
        # In demo mode, skip the heavy compute stages 02 and 03
        keys = [k for k in keys if k not in ("02", "03")]
        console.print("[yellow]Demo mode: skipping structure preparation and MD simulation.[/yellow]\n")

    console.print(f"Running stages: [bold cyan]{', '.join(keys)}[/bold cyan]\n")

    results: dict[str, tuple[bool, float]] = {}
    total_start = time.perf_counter()

    for key in keys:
        ok, elapsed = run_stage(key, demo=args.demo)
        results[key] = (ok, elapsed)
        if not ok and key in ("01", "02", "03") and not args.demo:
            console.print(f"\n[red]Critical stage {key} failed – aborting pipeline.[/red]")
            console.print("Tip: run with [bold]--demo[/bold] to generate synthetic data.")
            break

    total = time.perf_counter() - total_start
    print_summary(results)

    n_ok = sum(1 for ok, _ in results.values() if ok)
    console.print(
        f"\n[bold]{'All' if n_ok == len(results) else n_ok}/{len(results)} stages completed "
        f"in {total:.1f}s[/bold]"
    )
    console.print("\n[bold cyan]Launch dashboard:[/bold cyan]  streamlit run dashboard/app.py")


if __name__ == "__main__":
    main()
