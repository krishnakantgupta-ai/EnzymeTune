# ⚡ EnzymeTuner

> **Automated Computational Pipeline for Decoupled Enzyme Thermostabilization and Kinetic Optimization ($K_m$ & $V_{\max}$)**

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![PyTorch](https://img.shields.io/badge/PyTorch-2.0+-ee4c2c.svg)](https://pytorch.org/)
[![FoldX](https://img.shields.io/badge/Engine-FoldX_5.0-success.svg)](https://foldxsuite.crg.eu/)
[![AutoDock Vina](https://img.shields.io/badge/Docking-AutoDock_Vina-orange.svg)](https://vina.scripps.edu/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

---

## 📌 Overview

A classic hurdle in rational protein engineering is the **activity–stability trade-off**:
* **Catalytic machinery & gating loops** require dynamic flexibility to bind substrates, lower transition-state activation barriers ($\Delta G^\ddagger$), and clear products.
* **Global thermostability** typically demands structural rigidification, which often penalizes catalytic turnover ($k_{\text{cat}}$ or $V_{\max}$) when applied near active pockets.

**EnzymeTuner** solves this dilemma through **geometric shell partitioning**. By mapping your enzyme into three distinct physical shells relative to the catalytic center, it decouples scaffold rigidification from pocket reshaping and pairs the best independent mutations into synergistic multi-point variants.

---

## 🔬 Spatial Architecture

```text
                       [ ENZYME SCAFFOLD ]
               > 8.5 Å from active center (FoldX ΔΔG)
            ┌─────────────────────────────────────────┐
            │       [ KINETIC TUNING SHELL ]          │
            │   4.5–8.5 Å (Vina ΔΔG_bind + ESM-2)     │
            │      ┌───────────────────────────┐      │
            │      │     [ CATALYTIC CORE ]    │      │
            │      │  < 4.5 Å (STRICTLY LOCKED)│      │
            │      │         [Active]          │      │
            │      └───────────────────────────┘      │
            └─────────────────────────────────────────┘
# Create dedicated conda environment
conda create -n enzymetuner python=3.10 -y
conda activate enzymetuner

# Install scientific dependencies, OpenBabel, and AutoDock Vina
conda install -c conda-forge openbabel vina biopython pandas matplotlib -y
pip install torch fair-esm

### 🧬 FoldX 5.0 Installation & Setup

FoldX is proprietary software with **free licenses for academic and non-profit research**. Because redistribution of the binary is restricted, you must download it directly from the official portal.

#### Step 1: Request an Academic License
1. Visit the [FoldX Suite Registration Page](https://foldxsuite.crg.eu/)[cite: 8].
2. Register using your institutional email address (`.edu`, `.ac.*`, or institutional domain).
3. Confirm your account via the verification email and sign in.

#### Step 2: Download the Binary
1. Navigate to the **Downloads** section.
2. Download the version corresponding to your operating system:
   * **Linux:** `FoldX 5.0 (Linux 64-bit)`
   * **macOS:** `FoldX 5.0 (macOS 64-bit)`
3. Extract the downloaded archive (usually a `.tar.gz` or `.zip` file):
   ```bash
   tar -xvf foldx5Linux64.tar.gz

| Zone | Distance | Engineering Objective | Strategy |
| :--- | :---: | :--- | :--- |
| **Catalytic Core** | `0.0 – 4.5 Å` | **Preserve Chemistry** | Strictly conserved[cite: 8]. Locked to avoid catalytic inactivation[cite: 8]. |
| **Kinetic Shell** | `4.5 – 8.5 Å` | **Optimize $K_m$ & $V_{\max}$** | Reshape binding cleft & transition-state electrostatics[cite: 8]. Evaluated via baseline-calibrated AutoDock Vina[cite: 8]. |
| **Scaffold Zone** | `> 8.5 Å` | **Elevate $T_m$ & Half-life** | Maximize thermodynamic folding stability ($\Delta\Delta G_{\text{folding}} < 0$) using FoldX 5.0 with zero catalytic penalty[cite: 8]. |

python enzymetuner.py \
  --pdb 1XYZ_clean.pdb \
  --chain A \
  --active_site 35 52 108 \
  --ligand substrate.sdf \
  --foldx_bin /usr/local/bin/foldx \
  --max_foldx 40 \
  --max_vina 20 \
  --out single_candidates.csv \
  --out_combos combinatorial_designs.csv \
  --plot optimization_landscape.png

| Flag | Default | Description |
| :--- | :---: | :--- |
| `--pdb` | *Required* | Path to clean, input protein PDB file[cite: 9]. |
| `--active_site` | *Required* | Space-separated residue numbers for catalytic machinery (e.g., `35 52 108`)[cite: 9]. |
| `--chain` | `A` | Target protein chain identifier[cite: 9]. |
| `--ligand` | `None` | Path to ligand/transition-state analog file (`.sdf`, `.mol2`, or `.pdb`)[cite: 9]. |
| `--foldx_bin` | `foldx` | Path to the FoldX binary executable[cite: 9]. |
| `--vina_bin` | `vina` | Path to the AutoDock Vina binary executable[cite: 3]. |
| `--obabel_bin` | `obabel` | Path to the OpenBabel binary executable[cite: 3]. |
| `--max_foldx` | `30` | Number of top ESM-filtered scaffold variants to evaluate in FoldX[cite: 9]. |
| `--max_vina` | `15` | Number of top second-shell variants to evaluate via Vina docking[cite: 9]. |
| `--out` | `single_candidates.csv` | Output file for single-point mutations[cite: 9]. |
| `--out_combos`| `combinatorial_designs.csv`| Output file for multi-point combinatorial variants[cite: 9]. |
| `--plot` | `stability_vs_affinity.png`| Destination path for the 4-quadrant optimization landscape plot[cite: 9]. |
