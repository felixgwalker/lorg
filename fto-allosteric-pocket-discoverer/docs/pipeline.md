# FTO Pocket Discovery — Pipeline Guide

## Overview

```
01_acquire       Download & curate PDB structures
02_prepare       Fix, protonate, solvate systems
03_simulate      OpenMM NPT molecular dynamics
04_motion        RMSD · RMSF · PCA · clustering
05_pockets       Grid / fpocket cavity detection
06_ranking       Score and rank candidate pockets
07_conservation  Evolutionary conservation mapping
08_events        Conformational event detection
09_network       Residue interaction network
10_ligand        Pharmacophore + optional docking
```

## Quick Start

### 1 · Set up environment
```bash
conda env create -f environment.yml
conda activate fto-allosteric-pocket-discoverer
```

### 2 · Run the full pipeline
```bash
python run_pipeline.py
```

### 3 · Launch the dashboard
```bash
streamlit run dashboard/app.py
```

---

## Stage Details

### FR1 · Structure Acquisition (`01_acquire.py`)
Downloads 3–10 human FTO structures from RCSB PDB.  
Output: `data/structures/{PDB_ID}.pdb`, `data/structures/metadata.json`

### FR2 · Structure Preparation (`02_prepare.py`)
Uses **PDBFixer** to:
- Fill missing residues and atoms
- Add hydrogens at physiological pH (7.4)
- Solvate in a truncated octahedron TIP3P box
- Add Na⁺/Cl⁻ to 150 mM ionic strength
Output: `data/prepared/{PDB_ID}_prepared.pdb`, `{PDB_ID}_system.xml`

### FR3 · Molecular Dynamics (`03_simulate.py`)
Runs **OpenMM** NPT simulations:
- AMBER ff14SB protein force field
- TIP3P explicit water
- 2 fs timestep with H-mass repartitioning
- Monte Carlo barostat (1 atm)
- Langevin thermostat (310 K)
- 3 independent replicas per structure, 20–50 ns each
Output: `data/trajectories/{PDB_ID}_rep{N}.dcd`, `*.log`

### FR4 · Motion Analysis (`04_motion_analysis.py`)
Calculates with **MDAnalysis**:
- Backbone RMSD vs time
- Per-residue RMSF
- Interdomain distance and angle
- PCA (top 10 components)
- K-means clustering
Output: `results/{PDB_ID}_rmsd.csv`, `_rmsf.csv`, `_pca.npz`, `_clusters.csv`

### FR5 · Pocket Detection (`05_pocket_detection.py`)
Two modes:
- **fpocket** (preferred): `mdpocket` for trajectory analysis
- **grid** (fallback): Voronoi-based grid burial scoring
Output: `data/pockets/{PDB_ID}_pocket_trajectory.json`

### FR6 · Pocket Ranking (`06_pocket_ranking.py`)
Composite druggability score:
```
score = 0.30×persistence + 0.25×volume_norm + 0.20×druggability
      + 0.15×conservation + 0.10×enclosure
```
Output: `results/pocket_ranking.csv`, `results/pocket_scores.json`

### FR7 · Conservation (`07_conservation.py`)
Queries **ConSurf** REST API or calculates from MSA.  
Output: `data/conservation/fto_conservation.json`

### FR8 · Event Detection (`08_event_detection.py`)
Detects:
- Pocket opening / closing (volume threshold crossing)
- Domain rotation events (angle change > 5°)
- Salt bridge / H-bond breaking / forming
Output: `results/events.json`

### FR9 · Residue Network (`09_residue_network.py`)
Builds Cα contact map (< 8 Å) for open vs closed states.  
Identifies contacts that change when P1 opens.  
Output: `results/residue_network.json`

### FR10 · Ligand Module (`10_ligand.py`)
Generates pharmacophore features for the top-ranked pocket.  
Optionally runs AutoDock Vina docking.  
Output: `results/pharmacophore.json`, `results/docking_results.json`

---

## Data Directories

```
data/
  structures/     Raw PDB downloads
  prepared/       PDBFixer output
  trajectories/   DCD trajectory files
  pockets/        Per-frame pocket JSON
  conservation/   ConSurf scores
results/
  figures/        Publication-quality PNGs/SVGs
  animations/     MP4/GIF outputs
  reports/        Final PDF report
```

## Demo Mode

If no MD data exists, the dashboard runs in **demo mode** with
scientifically plausible synthetic data. Run stage 01 (acquire)
first to get real PDB structures in the Structure tab.
