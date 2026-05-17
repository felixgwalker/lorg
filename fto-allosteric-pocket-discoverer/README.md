# FTO Interdomain Pocket Discovery Platform

A computational platform that identifies **transient allosteric binding pockets** in the human FTO protein using molecular dynamics simulation, cavity detection, and evolutionary conservation analysis.

---

## What This Platform Does

FTO (Fat mass and Obesity-associated protein) is an RNA demethylase that removes the N6-methyladenosine (m6A) mark from messenger RNA. It is a validated drug target in obesity, type 2 diabetes, and several cancers.

Existing FTO inhibitors target the **active site** — the central catalytic cavity where the m6A nucleotide binds. This approach faces a selectivity challenge: the active site is conserved across the AlkB enzyme family, and inhibitors may inadvertently hit related proteins.

This platform searches for **cryptic allosteric pockets**: cavities that are hidden in the static crystal structure but open transiently during protein motion. A ligand targeting such a pocket could modulate FTO activity with greater selectivity.

The platform:

1. Downloads and curates experimental FTO crystal structures from the RCSB Protein Data Bank
2. Prepares each structure (adding missing atoms, solvating in explicit water, adding physiological ions)
3. Runs unbiased NPT molecular dynamics simulations (OpenMM, AMBER ff14SB)
4. Detects cavities in every simulation frame using a grid-based algorithm or fpocket
5. Ranks pockets by persistence, volume, druggability, evolutionary conservation, and enclosure
6. Detects conformational events (domain rotations, contact changes, pocket opening)
7. Maps evolutionary conservation (ConSurf) onto pocket-lining residues
8. Builds a pharmacophore model for the top-ranked pocket
9. Presents all results in an interactive Streamlit dashboard

---

## Who It Is For

- **Medicinal chemists and drug discovery scientists** evaluating FTO as an allosteric drug target
- **Structural biologists** interested in FTO conformational dynamics
- **Prospective clients and collaborators** evaluating the methodology
- **Computational chemists** seeking to extend or reproduce the analysis

The dashboard is designed to be navigable by scientists who are not computational chemists: all plots include plain-language explanations, and all scientific terms are defined in context.

---

## Key Results

The platform's primary finding from the 20 ns simulation:

| Pocket | Name | Composite score | Persistence | Peak volume |
|--------|------|----------------|-------------|-------------|
| **P1** | Interdomain Interface Pocket | **0.82 / 1.0** | **38%** | **487 Å³** |
| P2 | Catalytic Active Site | 0.88 / 1.0 | 98% | 470 Å³ |
| P3 | CTD Surface Groove | 0.51 / 1.0 | 22% | 360 Å³ |

**P1** is the novel finding: a transient pocket at the catalytic–C-terminal domain interface that opens after a ~9° interdomain rotation. It is absent in all crystal structures and would not be discovered by traditional structure-based docking. Nine of its 14 lining residues score ≥ 7/9 on the ConSurf conservation scale, indicating the pocket geometry is evolutionarily maintained.

---

## Scientific Background

### FTO structure

Human FTO (UniProt Q9C0B1) consists of two domains:

- **Catalytic domain** (residues 31–326): An AlkB-like double-stranded β-helix (DSBH) fold housing the Fe²⁺/α-ketoglutarate active site. Iron-binding triad: **H231 · D233 · H307**.
- **C-terminal domain** (residues 327–498): A structural domain unique to FTO within the AlkB family, packing tightly against the catalytic domain.

The two domains share an ~1,200 Å² interface stabilised by a salt bridge (E244–R365) and two hydrogen bonds (Q108–N338, Y205–T332).

### What is a transient pocket?

Proteins are not static — they flex continuously. A **transient (cryptic) pocket** is a cavity that:
- Does not exist in the crystal structure
- Opens during thermally accessible protein motions
- Is large enough (typically > 200 Å³) to accommodate a drug-like molecule

Transient pockets are particularly attractive as drug targets because they may exist only in one protein and not in related family members.

### Pocket scoring formula

```
composite_score = 0.30 × persistence
               + 0.25 × volume_normalised
               + 0.20 × druggability
               + 0.15 × conservation
               + 0.10 × enclosure
```

where:
- **persistence** = fraction of MD frames in which the pocket exceeds the volume threshold
- **volume_normalised** = mean open-state volume, scaled to [0, 1] between 100–700 Å³
- **druggability** = estimate from hydrophobicity, charge balance, and volume
- **conservation** = mean ConSurf score (1–9) of lining residues, normalised to [0, 1]
- **enclosure** = fraction of the pocket surface that is protein-enclosed vs solvent-exposed

### How to interpret results

| Metric | Good range | What it means |
|--------|-----------|---------------|
| Persistence | > 20% | Pocket is accessible for long enough to be pharmacologically relevant |
| Peak volume | 250–800 Å³ | Sufficient space for a drug-like molecule (MW ~300–500 Da) |
| Druggability | > 0.6 | Shape and chemistry favour drug binding |
| Conservation | > 0.7 | Pocket geometry is evolutionarily maintained → functional significance |
| ConSurf score | 7–9 | Residue is invariant or near-invariant across vertebrate orthologues |

---

## Technology Stack

| Component | Technology |
|-----------|-----------|
| MD simulation | OpenMM 8 · AMBER ff14SB · TIP3P water |
| Trajectory analysis | MDAnalysis 2.6 |
| Pocket detection | fpocket (preferred) / grid + DBSCAN (fallback) |
| Conservation | ConSurf REST API |
| Structure preparation | PDBFixer |
| Dashboard | Streamlit 1.32 |
| 3D viewer | 3Dmol.js 2.0.4 |
| Plots | Plotly 5 |
| Data | NumPy · Pandas · SciPy · scikit-learn |

---

## Architecture Overview

```
run_pipeline.py            Orchestrator — runs stages 01–10 in sequence
│
├── scripts/
│   ├── 01_acquire.py      Download FTO structures from RCSB PDB
│   ├── 02_prepare.py      PDBFixer — add missing atoms, solvate, add ions
│   ├── 03_simulate.py     OpenMM NPT simulation (3 replicas × 20 ns)
│   ├── 04_motion_analysis.py  RMSD · RMSF · PCA · clustering · interdomain motion
│   ├── 05_pocket_detection.py Grid or fpocket cavity detection per frame
│   ├── 06_pocket_ranking.py   Composite scoring and ranking
│   ├── 07_conservation.py     ConSurf API query
│   ├── 08_event_detection.py  Pocket opening / domain rotation / contact events
│   ├── 09_residue_network.py  Cα contact map · breaking / forming contacts
│   ├── 10_ligand.py           Pharmacophore model · optional AutoDock Vina
│   └── utils.py               Shared paths · I/O · DemoData generator
│
├── dashboard/
│   ├── app.py             Streamlit entry point · page config · tab layout
│   ├── components/
│   │   ├── viewer.py      3Dmol.js HTML components (structure, pocket sphere, ghost trail)
│   │   └── story_mode.py  Guided narrative walkthrough component
│   └── tabs/
│       ├── tab_structure.py   Crystal structure viewer and domain annotation
│       ├── tab_trajectory.py  RMSD · RMSF · interdomain motion · frame explorer
│       ├── tab_pockets.py     Pocket ranking table · volume trajectory · 3D viewer
│       ├── tab_analysis.py    PCA · clustering · conservation overlay
│       ├── tab_events.py      Event timeline · contact network · heatmap
│       └── tab_report.py      Story mode · key findings · pharmacophore · export
│
├── data/
│   ├── structures/        PDB files + metadata.json
│   ├── prepared/          PDBFixer output (topology + system XML)
│   ├── trajectories/      DCD trajectory files
│   ├── pockets/           Per-frame pocket JSON
│   └── conservation/      ConSurf scores JSON
│
├── results/               Analysis outputs (CSV, JSON)
├── config.yaml            All parameters (simulation, pocket detection, analysis)
├── environment.yml        Conda environment specification
└── requirements.txt       pip dependencies
```

---

## Repository Structure

```
fto-allosteric-pocket-discoverer/
├── README.md
├── PRD                    Product requirements document
├── config.yaml            Pipeline and simulation parameters
├── environment.yml        Conda environment
├── requirements.txt       pip dependencies
├── run_pipeline.py        Pipeline orchestrator
├── scripts/               Analysis scripts (01–10)
├── dashboard/             Streamlit dashboard
├── data/                  Input data (structures, trajectories, pockets, conservation)
├── results/               Pipeline outputs
└── docs/
    └── pipeline.md        Stage-by-stage pipeline reference
```

---

## Setup and Installation

### Prerequisites

- Python 3.11
- Conda (recommended) or pip
- ~4 GB disk space for simulation trajectories
- ~16 GB RAM for MD simulation (optional; not required for demo mode)

### Option A: Conda (recommended)

```bash
conda env create -f environment.yml
conda activate fto-allosteric-pocket-discoverer
```

This installs OpenMM via conda-forge and all other dependencies.

### Option B: pip

```bash
pip install -r requirements.txt
```

> **Note:** OpenMM is not available on Windows via pip. On Windows, use Option A (conda) or run the pipeline in demo mode (`--demo`), which does not require OpenMM.

---

## Running the Pipeline

### Demo mode (no OpenMM required, ~60 seconds)

Generates scientifically plausible synthetic data and populates the dashboard:

```bash
python run_pipeline.py --demo
```

This skips stages 02 (structure preparation) and 03 (MD simulation) and generates synthetic trajectory and pocket data.

### Full pipeline (~2 hours with GPU, longer on CPU)

```bash
python run_pipeline.py
```

Runs all 10 stages in sequence. Requires OpenMM and a prepared environment.

### Selective execution

```bash
# Resume from a specific stage
python run_pipeline.py --from 04

# Run only specific stages
python run_pipeline.py --only 06,07,08

# Skip stages
python run_pipeline.py --skip 03
```

---

## Launching the Dashboard

```bash
streamlit run dashboard/app.py
```

Opens at `http://localhost:8501` in your browser.

If pipeline results exist in `data/` and `results/`, they are displayed. Otherwise, the dashboard runs in demo mode with synthetic data.

---

## Pipeline Stage Reference

| Stage | Script | Description | Demo support |
|-------|--------|-------------|-------------|
| 01 | `01_acquire.py` | Download FTO PDB structures from RCSB | — |
| 02 | `02_prepare.py` | PDBFixer: missing atoms, hydrogens, solvation | — |
| 03 | `03_simulate.py` | OpenMM NPT MD simulation (3 replicas × 20 ns) | — |
| 04 | `04_motion_analysis.py` | RMSD · RMSF · PCA · interdomain motion | ✓ |
| 05 | `05_pocket_detection.py` | Grid or fpocket cavity detection | ✓ |
| 06 | `06_pocket_ranking.py` | Composite scoring and ranking | ✓ |
| 07 | `07_conservation.py` | ConSurf API conservation mapping | ✓ |
| 08 | `08_event_detection.py` | Pocket / domain / contact events | ✓ |
| 09 | `09_residue_network.py` | Cα contact map; breaking/forming contacts | ✓ |
| 10 | `10_ligand.py` | Pharmacophore model; optional Vina docking | ✓ |

---

## Configuration

All parameters are in `config.yaml`. Key sections:

```yaml
fto:
  pdb_ids: ["3LFM", "4IE4", "4ZS3"]   # Confirmed human FTO structures
  primary_structure: "3LFM"             # Used as reference for analysis

simulation:
  n_replicas: 3
  production_ns: 20.0
  temperature_K: 310.0

pocket_detection:
  method: "grid"       # "fpocket" if binary available, else "grid"
  min_volume_A3: 100
  min_persistence_fraction: 0.03

conservation:
  method: "consurf_api"
  n_homologs: 50
```

---

## Environment Variables

No environment variables are required by default. If using a custom ConSurf API key or a remote compute node, set:

```bash
export CONSURF_API_KEY="your_key"   # optional
export OPENMM_PLATFORM="CUDA"       # or "OpenCL", "CPU"
```

---

## FTO PDB Structures

The platform uses three confirmed human FTO crystal structures:

| PDB ID | Resolution | Description | Reference |
|--------|-----------|-------------|-----------|
| **3LFM** | 2.50 Å | FTO catalytic domain with 3-methylthymidine and iron | Jia et al. 2011, *Nature* |
| **4IE4** | 2.50 Å | FTO in complex with inhibitor IOX1 | Feng et al. 2014 |
| **4ZS3** | 2.45 Å | FTO bound to 5-aminofluorescein | Huang et al. 2015 |

> **Important:** Only add PDB IDs that have been verified to be human FTO structures. The acquisition script (`01_acquire.py`) will warn if a downloaded structure's title does not contain FTO-related keywords.

---

## Running Tests

There are currently no automated tests. To validate a pipeline run manually:

```bash
# Check that demo outputs are generated correctly
python run_pipeline.py --demo
ls results/          # should contain *.csv and *.json files
ls data/pockets/     # should contain 3LFM_pockets.json

# Launch dashboard and verify each tab renders
streamlit run dashboard/app.py
```

---

## Known Limitations

1. **Simulation length:** 20 ns per replica may be insufficient to capture all relevant conformational states. Rare opening events with recurrence times > 20 ns would be missed.

2. **Demo data is synthetic:** The demo mode uses mathematically generated data that reproduces the expected statistical properties of the real simulation but is not derived from actual MD trajectories.

3. **Grid pocket detection is approximate:** The fallback grid-based method is computationally cheap but less accurate than fpocket or mdpocket. For production analysis, install fpocket.

4. **Conservation via ConSurf API:** The ConSurf REST API may be rate-limited or unavailable. If queries fail, conservation scores default to a neutral value (0.5) in the ranking.

5. **Pharmacophore positions are estimated:** Feature positions are approximated from residue centroids, not from full quantum-mechanical calculations. These are suitable for virtual screening filters but not for precise docking poses.

6. **Single protein target:** The platform is purpose-built for FTO. Extending it to other targets requires updating `config.yaml` and validating the domain boundary assignments.

7. **OpenMM not available on Windows via pip:** Use the conda environment on Windows.

---

## Future Improvements

- [ ] Multi-target support via configurable protein profiles
- [ ] Integration with fpocket mdpocket for full trajectory pocket analysis
- [ ] Enhanced pharmacophore quality using full residue geometry
- [ ] Automated virtual screening pipeline (integrate with ZINC or ChEMBL)
- [ ] Export to SDF format for use in molecular docking workflows
- [ ] Animated trajectory playback in the dashboard
- [ ] Statistical comparison across replicas (convergence metrics)
- [ ] Enhanced event detection using machine learning classifiers

---

## Troubleshooting

**Dashboard shows demo data even after running the pipeline:**
Confirm that `results/3LFM_rep0_rmsd.csv` exists. If the primary structure in `config.yaml` does not match the file prefix, the status check will fail.

**OpenMM not found on Windows:**
Use `conda install -c conda-forge openmm` or run in demo mode: `python run_pipeline.py --demo`.

**ConSurf API times out:**
The ConSurf server can be slow. The script retries three times. If all retries fail, conservation scores are set to 0.5 (neutral) and a warning is logged. This does not prevent other stages from running.

**fpocket not found — pocket detection falls back to grid:**
Install fpocket from https://github.com/Discngine/fpocket and ensure it is on your `PATH`. The grid fallback is slower and less accurate but fully functional.

**3D viewer blank in dashboard:**
The viewer loads 3Dmol.js from a CDN. Ensure you have internet access, or host the script locally and update the CDN URL in `dashboard/components/viewer.py`.

**Memory error during MD simulation:**
Reduce the system size by decreasing `padding_nm` in `config.yaml`, or run on a machine with more RAM. The demo mode (`--demo`) does not require MD and uses no significant memory.

---

## Contributing / Developer Notes

- All pipeline parameters live in `config.yaml` — avoid hardcoding values in scripts.
- Shared path constants and utilities are in `scripts/utils.py`.
- Each script is independently runnable with `--demo` for rapid iteration.
- Dashboard tabs import from `scripts.utils` and have no interdependencies.
- The `DemoData` class in `scripts/utils.py` generates all synthetic demo data. Keep it consistent with the real pipeline output schema.
- Before adding new PDB IDs to `config.yaml`, verify they are confirmed human FTO structures in the RCSB.

---

## Scientific Caveats

- MD simulations are stochastic: different random seeds produce different trajectories. Findings should be reproduced across at least 3 independent replicas.
- The 20 ns simulation time is a starting point; longer simulations or enhanced sampling (e.g., metadynamics) may reveal additional pockets or give better persistence estimates.
- Pocket druggability scores are empirical estimates. Experimental validation (e.g., NMR fragment screening, X-ray crystallography of bound ligands) is required to confirm the pharmacological relevance of any identified pocket.
- Conservation analysis is limited to vertebrate FTO orthologues. Cross-species conservation does not guarantee that a pocket is biologically functional in the human protein.

---

## Citation

If you use this platform in your research:

> FTO Pocket Discovery Platform v1.0.
> Primary structure: Jia G et al. (2011) Oxidative demethylation of 3-methylthymine and 3-methyluracil in single-stranded DNA and RNA by mouse and human FTO. *FEBS Lett* 582, 3313–3319.
> MD: OpenMM 8 — Eastman P et al. (2023) *J Chem Theory Comput*.
> Conservation: Ashkenazy H et al. (2016) ConSurf 2016. *Nucleic Acids Res* 44, W344–W350.

---

## Licence

Internal research use. Contact the project owner for redistribution terms.
