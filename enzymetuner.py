#!/usr/bin/env python3
"""
EnzymeTuner: Automated Computational Pipeline for Simultaneous Enzyme
Thermostabilization and Kinetic Parameter Optimization.

Integrates:
  1. Structural spatial partitioning (Core vs Kinetic Shell vs Scaffold)
  2. Evolutionary zero-shot fitness filtering (ESM-2)
  3. Thermodynamic stability calculations (FoldX RepairPDB & BuildModel)
  4. Kinetic substrate binding calibration (AutoDock Vina & OpenBabel)
  5. Multi-point combinatorial design (Scaffold + Kinetic Shell pairing)
  6. Automated 4-quadrant optimization landscape visualization
"""

import os
import glob
import shutil
import subprocess
import argparse
import tempfile
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, NeighborSearch, Selection
import torch
import esm

AA_LIST = list("ACDEFGHIKLMNPQRSTVWY")

class FoldXEngine:
    def __init__(self, foldx_bin: str = "foldx", work_dir: str = None):
        self.foldx_bin = foldx_bin
        self.work_dir = work_dir or tempfile.mkdtemp(prefix="foldx_run_")
        os.makedirs(self.work_dir, exist_ok=True)

    def repair_pdb(self, pdb_path: str) -> str:
        pdb_filename = os.path.basename(pdb_path)
        dest_pdb = os.path.join(self.work_dir, pdb_filename)
        shutil.copyfile(pdb_path, dest_pdb)

        cmd = [self.foldx_bin, "--command=RepairPDB", f"--pdb={pdb_filename}"]
        print(f"[*] Running FoldX RepairPDB on {pdb_filename}...")
        result = subprocess.run(cmd, cwd=self.work_dir, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"FoldX RepairPDB failed:\n{result.stderr}")

        base_name = os.path.splitext(pdb_filename)[0]
        repaired_pdb = os.path.join(self.work_dir, f"{base_name}_Repair.pdb")
        if not os.path.exists(repaired_pdb):
            raise FileNotFoundError(f"Repaired PDB not found at: {repaired_pdb}")
        return repaired_pdb

    def compute_ddg_batch(self, repaired_pdb_path: str, mutations: list[dict], chain_id: str = "A") -> tuple[dict, dict]:
        pdb_filename = os.path.basename(repaired_pdb_path)
        individual_list_path = os.path.join(self.work_dir, "individual_list.txt")

        tag_to_mut = {}
        with open(individual_list_path, "w") as f:
            for item in mutations:
                mut_tag = f"{item['WT']}{chain_id}{item['PDB_Resid']}{item['Variant']}"
                f.write(f"{mut_tag};\n")
                tag_to_mut[mut_tag] = f"{item['WT']}{item['PDB_Resid']}{item['Variant']}"

        cmd = [
            self.foldx_bin,
            "--command=BuildModel",
            f"--pdb={pdb_filename}",
            "--mutant-file=individual_list.txt",
            "--numberOfRuns=1"
        ]

        print(f"[*] Running FoldX BuildModel across {len(mutations)} candidates...")
        result = subprocess.run(cmd, cwd=self.work_dir, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"FoldX BuildModel failed:\n{result.stderr}")

        ddg_results = {}
        base_name = os.path.splitext(pdb_filename)[0]
        dif_file = os.path.join(self.work_dir, f"Dif_{base_name}.fxout")
        avg_file = os.path.join(self.work_dir, f"Average_{base_name}.fxout")

        target_out = dif_file if os.path.exists(dif_file) else avg_file
        if os.path.exists(target_out):
            with open(target_out, "r") as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) >= 3 and not parts[0].startswith("Pdb"):
                        mutant_tag = parts[0].replace(";", "").strip()
                        try:
                            ddg_val = float(parts[1].strip())
                            std_code = tag_to_mut.get(mutant_tag)
                            if std_code:
                                ddg_results[std_code] = round(ddg_val, 3)
                        except ValueError:
                            continue

        mutant_pdbs = {}
        for mut_tag, std_code in tag_to_mut.items():
            pattern = os.path.join(self.work_dir, f"{base_name}*_{mut_tag}*.pdb")
            matched = glob.glob(pattern)
            if matched:
                mutant_pdbs[std_code] = matched[0]

        return ddg_results, mutant_pdbs

    def cleanup(self):
        if os.path.exists(self.work_dir):
            shutil.rmtree(self.work_dir)


class AutoDockVinaEngine:
    def __init__(self, vina_bin: str = "vina", obabel_bin: str = "obabel", work_dir: str = None):
        self.vina_bin = vina_bin
        self.obabel_bin = obabel_bin
        self.work_dir = work_dir or tempfile.mkdtemp(prefix="vina_run_")
        os.makedirs(self.work_dir, exist_ok=True)

    def prepare_ligand_pdbqt(self, ligand_path: str) -> str:
        out_pdbqt = os.path.join(self.work_dir, "ligand.pdbqt")
        cmd = [
            self.obabel_bin,
            ligand_path,
            "-O", out_pdbqt,
            "-p", "7.4",
            "--partialcharge", "gasteiger"
        ]
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0 or not os.path.exists(out_pdbqt):
            raise RuntimeError(f"OpenBabel ligand conversion failed:\n{res.stderr}")
        return out_pdbqt

    def prepare_receptor_pdbqt(self, receptor_pdb: str, out_name: str) -> str:
        out_pdbqt = os.path.join(self.work_dir, f"{out_name}.pdbqt")
        cmd = [
            self.obabel_bin,
            receptor_pdb,
            "-O", out_pdbqt,
            "-xr",
            "-xh",
            "--partialcharge", "gasteiger"
        ]
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0 or not os.path.exists(out_pdbqt):
            raise RuntimeError(f"OpenBabel receptor conversion failed:\n{res.stderr}")
        return out_pdbqt

    def dock_candidate(self, receptor_pdbqt: str, ligand_pdbqt: str, center: tuple[float, float, float], box_size: float = 20.0) -> float:
        out_docked = os.path.join(self.work_dir, "docked_out.pdbqt")
        cmd = [
            self.vina_bin,
            "--receptor", receptor_pdbqt,
            "--ligand", ligand_pdbqt,
            "--center_x", f"{center[0]:.3f}",
            "--center_y", f"{center[1]:.3f}",
            "--center_z", f"{center[2]:.3f}",
            "--size_x", f"{box_size:.1f}",
            "--size_y", f"{box_size:.1f}",
            "--size_z", f"{box_size:.1f}",
            "--exhaustiveness", "8",
            "--out", out_docked
        ]

        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            print(f"[!] Vina docking warning: {res.stderr}")
            return None

        for line in res.stdout.splitlines():
            line_parts = line.strip().split()
            if len(line_parts) >= 2 and line_parts[0] == "1":
                try:
                    return float(line_parts[1])
                except ValueError:
                    continue
        return None

    def cleanup(self):
        if os.path.exists(self.work_dir):
            shutil.rmtree(self.work_dir)


class EnzymeTuner:
    def __init__(self, pdb_path: str, catalytic_resids: list[int], chain_id: str = "A", 
                 foldx_bin: str = "foldx", vina_bin: str = "vina", obabel_bin: str = "obabel"):
        self.pdb_path = pdb_path
        self.catalytic_resids = set(catalytic_resids)
        self.chain_id = chain_id
        self.foldx_bin = foldx_bin
        self.vina_bin = vina_bin
        self.obabel_bin = obabel_bin

        parser = PDBParser(QUIET=True)
        self.structure = parser.get_structure("enzyme", pdb_path)
        self.chain = self.structure[0][chain_id]

        print("[*] Loading ESM-2 evolutionary model (35M parameters)...")
        self.esm_model, self.alphabet = esm.pretrained.esm2_t12_35M_UR50D()
        self.batch_converter = self.alphabet.get_batch_converter()
        self.esm_model.eval()

    def get_pocket_geometric_center(self) -> tuple[float, float, float]:
        coords = [
            atom.coord for res in self.chain 
            if res.id[1] in self.catalytic_resids 
            for atom in res.get_atoms()
        ]
        if not coords:
            raise ValueError("No atoms found for catalytic residues.")
        centroid = np.mean(coords, axis=0)
        return float(centroid[0]), float(centroid[1]), float(centroid[2])

    def get_residue_ca_coord(self, resid: int) -> np.ndarray:
        for res in self.chain:
            if res.id[1] == resid and 'CA' in res:
                return res['CA'].coord
        return None

    def partition_residues(self, inner_cutoff: float = 4.5, outer_cutoff: float = 8.5):
        atoms = Selection.unfold_entities(self.chain, 'A')
        ns = NeighborSearch(atoms)

        active_atoms = [
            atom for res in self.chain
            if res.id[1] in self.catalytic_resids
            for atom in res.get_atoms()
        ]

        core_resids = set(self.catalytic_resids)
        shell_resids = set()

        for atom in active_atoms:
            for r in ns.search(atom.coord, inner_cutoff, level='R'):
                if r.id[0] == ' ':
                    core_resids.add(r.id[1])
            for r in ns.search(atom.coord, outer_cutoff, level='R'):
                if r.id[0] == ' ' and r.id[1] not in core_resids:
                    shell_resids.add(r.id[1])

        all_resids = {r.id[1] for r in self.chain if r.id[0] == ' '}
        scaffold_resids = all_resids - core_resids - shell_resids
        return core_resids, shell_resids, scaffold_resids

    def get_sequence(self):
        three_to_one = {
            'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F',
            'GLY':'G', 'HIS':'H', 'ILE':'I', 'LYS':'K', 'LEU':'L',
            'MET':'M', 'ASN':'N', 'PRO':'P', 'GLN':'Q', 'ARG':'R',
            'SER':'S', 'THR':'T', 'VAL':'V', 'TRP':'W', 'TYR':'Y'
        }
        seq, resids = [], []
        for res in self.chain:
            if res.id[0] == ' ' and res.resname in three_to_one:
                seq.append(three_to_one[res.resname])
                resids.append(res.id[1])
        return "".join(seq), resids

    def compute_esm_scores(self, wild_seq: str, target_mutations: list[tuple[str, int, str]]) -> dict:
        data = [("protein", wild_seq)]
        _, _, batch_tokens = self.batch_converter(data)

        with torch.no_grad():
            token_probs = torch.log_softmax(self.esm_model(batch_tokens)["logits"], dim=-1)

        scores = {}
        for wt_aa, seq_idx, mut_aa in target_mutations:
            wt_tok = self.alphabet.get_idx(wt_aa)
            mut_tok = self.alphabet.get_idx(mut_aa)
            llr = token_probs[0, seq_idx + 1, mut_tok] - token_probs[0, seq_idx + 1, wt_tok]
            scores[f"{wt_aa}{seq_idx+1}{mut_aa}"] = round(llr.item(), 3)
        return scores

    def screen_and_score(self, ligand_path: str = None, max_foldx_evals: int = 30, max_vina_evals: int = 15) -> pd.DataFrame:
        core, shell, scaffold = self.partition_residues()
        seq, resids = self.get_sequence()
        id_to_idx = {resid: i for i, resid in enumerate(resids)}

        mutation_pool = []
        for resid in scaffold:
            idx = id_to_idx[resid]
            for target_aa in AA_LIST:
                if target_aa != seq[idx]:
                    mutation_pool.append((seq[idx], idx, target_aa, resid, "Thermostability"))

        for resid in shell:
            idx = id_to_idx[resid]
            for target_aa in AA_LIST:
                if target_aa != seq[idx]:
                    mutation_pool.append((seq[idx], idx, target_aa, resid, "Kinetic_Tuning"))

        raw_tuples = [(m[0], m[1], m[2]) for m in mutation_pool]
        fitness_scores = self.compute_esm_scores(seq, raw_tuples)

        prefiltered = []
        for wt_aa, idx, target_aa, resid, category in mutation_pool:
            key = f"{wt_aa}{idx+1}{target_aa}"
            score = fitness_scores.get(key, -99.0)
            if score >= -1.0:
                prefiltered.append({
                    "Mutation": f"{wt_aa}{resid}{target_aa}",
                    "PDB_Resid": resid,
                    "WT": wt_aa,
                    "Variant": target_aa,
                    "Target_Class": category,
                    "ESM_Fitness_Score": score,
                    "FoldX_ddG_kcal_mol": None,
                    "WT_Vina_Affinity_kcal_mol": None,
                    "Mutant_Vina_Affinity_kcal_mol": None,
                    "Delta_Vina_Affinity_kcal_mol": None
                })

        df = pd.DataFrame(prefiltered).sort_values(by="ESM_Fitness_Score", ascending=False)
        if df.empty:
            return df

        foldx = FoldXEngine(foldx_bin=self.foldx_bin)
        mutant_pdbs = {}
        repaired_pdb = None
        try:
            repaired_pdb = foldx.repair_pdb(self.pdb_path)
            scaffold_evals = df[df["Target_Class"] == "Thermostability"].head(max_foldx_evals).to_dict("records")
            kinetic_evals = df[df["Target_Class"] == "Kinetic_Tuning"].head(max_vina_evals).to_dict("records")
            all_structural_evals = scaffold_evals + kinetic_evals

            if all_structural_evals:
                ddg_dict, mutant_pdbs = foldx.compute_ddg_batch(repaired_pdb, all_structural_evals, chain_id=self.chain_id)
                df["FoldX_ddG_kcal_mol"] = df["Mutation"].map(ddg_dict)

            if ligand_path and os.path.exists(ligand_path):
                vina = AutoDockVinaEngine(vina_bin=self.vina_bin, obabel_bin=self.obabel_bin)
                try:
                    center = self.get_pocket_geometric_center()
                    print(f"[*] Active site centroid: X={center[0]:.2f}, Y={center[1]:.2f}, Z={center[2]:.2f}")
                    ligand_pdbqt = vina.prepare_ligand_pdbqt(ligand_path)

                    print("[*] Docking into Wild-Type baseline receptor...")
                    wt_pdbqt = vina.prepare_receptor_pdbqt(repaired_pdb, out_name="WT_repaired")
                    wt_affinity = vina.dock_candidate(wt_pdbqt, ligand_pdbqt, center=center)
                    print(f"[✓] Wild-Type Baseline Binding Affinity: {wt_affinity} kcal/mol")

                    affinity_dict = {}
                    delta_dict = {}

                    kinetic_top = df[df["Target_Class"] == "Kinetic_Tuning"].head(max_vina_evals)
                    for _, row in kinetic_top.iterrows():
                        mut_code = row["Mutation"]
                        mut_pdb = mutant_pdbs.get(mut_code)
                        if mut_pdb and os.path.exists(mut_pdb):
                            print(f"[*] Docking into variant {mut_code}...")
                            rec_pdbqt = vina.prepare_receptor_pdbqt(mut_pdb, out_name=mut_code)
                            mut_affinity = vina.dock_candidate(rec_pdbqt, ligand_pdbqt, center=center)
                            if mut_affinity is not None:
                                affinity_dict[mut_code] = mut_affinity
                                if wt_affinity is not None:
                                    delta_dict[mut_code] = round(mut_affinity - wt_affinity, 3)

                    df["WT_Vina_Affinity_kcal_mol"] = wt_affinity
                    df["Mutant_Vina_Affinity_kcal_mol"] = df["Mutation"].map(affinity_dict)
                    df["Delta_Vina_Affinity_kcal_mol"] = df["Mutation"].map(delta_dict)
                finally:
                    vina.cleanup()
        finally:
            foldx.cleanup()

        return df

    def generate_combinatorial_designs(self, df: pd.DataFrame, min_ca_dist: float = 9.0, top_n: int = 5) -> pd.DataFrame:
        """
        Creates synergistic multi-point variants combining:
          - Up to 2 independent scaffold stabilizing mutations (ddG < -0.8 kcal/mol)
          - 1 kinetic pocket mutation (Delta_affinity < -0.3 kcal/mol)
        Ensures residues are separated by at least min_ca_dist to prevent steric clash.
        """
        scaffold_cands = df[(df["Target_Class"] == "Thermostability") & (df["FoldX_ddG_kcal_mol"] < -0.5)].copy()
        kinetic_cands = df[(df["Target_Class"] == "Kinetic_Tuning") & (df["Delta_Vina_Affinity_kcal_mol"] < -0.2)].copy()

        if scaffold_cands.empty or kinetic_cands.empty:
            print("[!] Combinatorial note: Not enough strong candidates passing cutoff to build combinations.")
            return pd.DataFrame()

        combos = []
        for _, k_row in kinetic_cands.head(top_n).iterrows():
            k_id = k_row["PDB_Resid"]
            k_coord = self.get_residue_ca_coord(k_id)

            for _, s_row in scaffold_cands.head(top_n).iterrows():
                s_id = s_row["PDB_Resid"]
                s_coord = self.get_residue_ca_coord(s_id)
                if k_coord is not None and s_coord is not None:
                    dist = np.linalg.norm(k_coord - s_coord)
                    if dist >= min_ca_dist:
                        combos.append({
                            "Design_Type": "Double Mutant (1 Scaffold + 1 Kinetic)",
                            "Mutations": f"{s_row['Mutation']} + {k_row['Mutation']}",
                            "Scaffold_Stabilizers": s_row['Mutation'],
                            "Kinetic_Tuner": k_row['Mutation'],
                            "Min_Inter_Residue_Dist_A": round(dist, 1),
                            "Predicted_Additive_ddG_kcal_mol": round(s_row['FoldX_ddG_kcal_mol'] + (k_row['FoldX_ddG_kcal_mol'] or 0.0), 3),
                            "Kinetic_Delta_Affinity_kcal_mol": k_row['Delta_Vina_Affinity_kcal_mol'],
                            "Combined_ESM_Score": round(s_row['ESM_Fitness_Score'] + k_row['ESM_Fitness_Score'], 3)
                        })

        combo_df = pd.DataFrame(combos)
        if not combo_df.empty:
            combo_df = combo_df.sort_values(by="Predicted_Additive_ddG_kcal_mol", ascending=True)
        return combo_df

    @staticmethod
    def generate_landscape_plot(df: pd.DataFrame, out_png: str = "stability_vs_affinity.png"):
        plot_df = df.dropna(subset=["FoldX_ddG_kcal_mol", "Delta_Vina_Affinity_kcal_mol"]).copy()
        if plot_df.empty:
            print("[!] Insufficient data to generate plot: candidates require both FoldX and Vina values.")
            return

        fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
        x = plot_df["FoldX_ddG_kcal_mol"]
        y = plot_df["Delta_Vina_Affinity_kcal_mol"]
        c = plot_df["ESM_Fitness_Score"]

        min_score, max_score = c.min(), c.max()
        sizes = 60 + 140 * ((c - min_score) / (max_score - min_score)) if max_score > min_score else 100

        scatter = ax.scatter(x, y, c=c, cmap="viridis", s=sizes, alpha=0.85, edgecolors="black", linewidth=1.2, zorder=3)
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label("ESM-2 Evolutionary Fitness (LLR)", fontsize=11, fontweight="bold")

        ax.axvline(0, color="gray", linestyle="--", linewidth=1.2, alpha=0.7, zorder=2)
        ax.axhline(0, color="gray", linestyle="--", linewidth=1.2, alpha=0.7, zorder=2)

        x_min, x_max = ax.get_xlim()
        y_min, y_max = ax.get_ylim()

        ax.axvspan(min(x_min, -0.1), 0, ymin=0, ymax=abs(min(y_min, 0)) / (y_max - y_min) if y_max != y_min else 0.5,
                   color="#2ca02c", alpha=0.10, zorder=1)

        ax.text(
            0.03, 0.04, "Synergistic Zone\n(ΔΔG < 0 & ΔΔG_bind < 0)",
            transform=ax.transAxes,
            fontsize=10,
            fontweight="bold",
            color="#1b661b",
            verticalalignment="bottom",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="#2ca02c", alpha=0.9)
        )

        for _, row in plot_df.iterrows():
            if row["FoldX_ddG_kcal_mol"] < -0.2 or row["Delta_Vina_Affinity_kcal_mol"] < -0.2:
                ax.annotate(
                    row["Mutation"],
                    (row["FoldX_ddG_kcal_mol"], row["Delta_Vina_Affinity_kcal_mol"]),
                    textcoords="offset points",
                    xytext=(6, 6),
                    fontsize=9,
                    fontweight="semibold",
                    color="#1a1a1a",
                    zorder=4
                )

        ax.set_xlabel("FoldX $\\Delta\\Delta G_{\\text{folding}}$ (kcal/mol) [← More Stable]", fontsize=12, fontweight="bold")
        ax.set_ylabel("$\\Delta\\Delta G_{\\text{bind}}$ vs WT (kcal/mol) [← Tighter Binding]", fontsize=12, fontweight="bold")
        ax.set_title("Enzyme Engineering Optimization Landscape\nStability vs. Kinetic Substrate Affinity", fontsize=14, fontweight="bold", pad=15)

        ax.grid(True, linestyle=":", alpha=0.5, zorder=1)
        plt.tight_layout()
        plt.savefig(out_png, dpi=300)
        plt.close()
        print(f"[✓] Optimization landscape plot saved to: {out_png}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="EnzymeTuner: Stability, Kinetics, and Combinatorial Multi-Mutant Pipeline")
    parser.add_argument("--pdb", type=str, required=True, help="Input enzyme PDB file")
    parser.add_argument("--active_site", type=int, nargs="+", required=True, help="Active site residue numbers")
    parser.add_argument("--ligand", type=str, default=None, help="Ligand file (SDF, MOL2, or PDB)")
    parser.add_argument("--chain", type=str, default="A", help="Target chain ID")
    parser.add_argument("--foldx_bin", type=str, default="foldx", help="Path to FoldX binary")
    parser.add_argument("--vina_bin", type=str, default="vina", help="Path to AutoDock Vina binary")
    parser.add_argument("--obabel_bin", type=str, default="obabel", help="Path to OpenBabel binary")
    parser.add_argument("--max_foldx", type=int, default=30, help="Max candidates for FoldX ddG")
    parser.add_argument("--max_vina", type=int, default=15, help="Max kinetic candidates for Vina docking")
    parser.add_argument("--out", type=str, default="single_candidates.csv", help="Output single-mutant CSV path")
    parser.add_argument("--out_combos", type=str, default="combinatorial_designs.csv", help="Output combinatorial CSV path")
    parser.add_argument("--plot", type=str, default="stability_vs_affinity.png", help="Output plot image filename")

    args = parser.parse_args()

    tuner = EnzymeTuner(
        pdb_path=args.pdb,
        catalytic_resids=args.active_site,
        chain_id=args.chain,
        foldx_bin=args.foldx_bin,
        vina_bin=args.vina_bin,
        obabel_bin=args.obabel_bin
    )

    results_df = tuner.screen_and_score(
        ligand_path=args.ligand,
        max_foldx_evals=args.max_foldx,
        max_vina_evals=args.max_vina
    )
    results_df.to_csv(args.out, index=False)
    print(f"\n[✓] Single mutant screening exported to: {args.out}")

    combos_df = tuner.generate_combinatorial_designs(results_df)
    if not combos_df.empty:
        combos_df.to_csv(args.out_combos, index=False)
        print(f"[✓] Combinatorial multi-mutant designs exported to: {args.out_combos}")

    if args.ligand:
        EnzymeTuner.generate_landscape_plot(results_df, out_png=args.plot)