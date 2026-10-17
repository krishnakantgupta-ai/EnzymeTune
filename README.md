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
