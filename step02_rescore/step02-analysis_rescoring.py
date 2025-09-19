#!/usr/bin/env python3
import os
import sys
import numpy as np
from rdkit import Chem
from collections import defaultdict
import pandas as pd
import csv

# CLI args
silcsmcdir = sys.argv[1]      # project/work directory
name = sys.argv[2]            # short name for output grouping
npose = 5                     # number of poses per molecule

class SILCSMCProcessor:
    def __init__(self, minconfdir):
        self.minconfdir = minconfdir

    def merge_sdfs(self, output_sdf, lgfe_csv_path):
        # Read and filter ligands
        lgfes = pd.read_csv(lgfe_csv_path)
        if 'SMILES' in lgfes.columns:
            lgfes.drop(columns='SMILES', inplace=True)
        if 'Ligand' not in lgfes.columns or 'LGFE(kcal/mol)' not in lgfes.columns:
            raise ValueError("Expected columns 'Ligand' and 'LGFE(kcal/mol)' not found in CSV.")

        filtered_df = lgfes[lgfes['LGFE(kcal/mol)'] <= -7]
        allowed_ids = set(filtered_df['Ligand'].astype(str))
        print(f"{len(allowed_ids)} ligands passed LGFE filter")

        count_written = 0
        with open(output_sdf, 'w') as outfile:
            for fn in os.listdir(self.minconfdir):
                if not fn.endswith('.sdf'):
                    continue
                prefix = fn.rsplit('.', 1)[0]
                base_id = prefix.split('_')[0]
                if base_id in allowed_ids:
                    sdf_path = os.path.join(self.minconfdir, fn)
                    with open(sdf_path, 'r') as infile:
                        outfile.write(infile.read())
                    count_written += 1
        print(f"Merged {count_written} molecules into {output_sdf}")

class SDFProcessor:
    def __init__(self, sdf_file):
        # Load molecules (keep Hs)
        self.sdf_supplier = Chem.SDMolSupplier(sdf_file, removeHs=False)
        self.molecules    = [m for m in self.sdf_supplier if m is not None]
        # group into batches of npose
        self.pose_groups  = [ self.molecules[i:i+npose]
                              for i in range(0, len(self.molecules), npose) ]

    def calculate_rmsd(self, ref, target):
        """RMSD between two conformers (same atom ordering)."""
        c1 = [ ref.GetConformer().GetAtomPosition(i) for i in range(ref.GetNumAtoms()) ]
        c2 = [ target.GetConformer().GetAtomPosition(i) for i in range(target.GetNumAtoms()) ]
        return np.sqrt( sum((a.x-b.x)**2+(a.y-b.y)**2+(a.z-b.z)**2
                             for a,b in zip(c1,c2)) / len(c1) )

    def set_properties(self):
        all_pose_data = []
        for group in self.pose_groups:
            ref = group[0]
            rmsds = []
            for mol in group:
                rmsd = self.calculate_rmsd(ref, mol) if mol != ref else 0.0
                rmsds.append(rmsd)
            rmsd_std = np.std(rmsds)

            for mol, rmsd in zip(group, rmsds):
                lgfe = float(mol.GetProp("LGFE")) if mol.HasProp("LGFE") else 0.0
                le = float(mol.GetProp("LE")) if mol.HasProp("LE") else 0.0
                all_pose_data.append({
                    "mol": mol,
                    "Compound_ID": mol.GetProp("_Name"),
                    "LGFE": lgfe,
                    "LE": le,
                    "RMSD": rmsd,
                    "RMSD_STD": rmsd_std,
                    "SMILES": Chem.MolToSmiles(mol)
                })

        # Weighted rescoring (equal weights)
        w_lgfe, w_le, w_rmsd = 1, 1, 1
        all_data = []
        for i, d in enumerate(all_pose_data):
            wtd_LGFE = w_lgfe * d["LGFE"]
            wtd_LE = w_le * d["LE"]
            wtd_RMSD = w_rmsd * d["RMSD"]
            rescore = wtd_LGFE + wtd_LE + wtd_RMSD
            mol = d["mol"]

            mol.SetProp("Weighted_LGFE", f"{wtd_LGFE:.4f}")
            mol.SetProp("Weighted_LE", f"{wtd_LE:.4f}")
            mol.SetProp("Weighted_RMSD", f"{wtd_RMSD:.4f}")
            mol.SetProp("Rescore_Total", f"{rescore:.4f}")
            mol.SetProp("LGFE", f"{d['LGFE']:.4f}")
            mol.SetProp("LE", f"{d['LE']:.4f}")
            mol.SetProp("RMSD", f"{d['RMSD']:.4f}")
            mol.SetProp("RMSD_STD", f"{d['RMSD_STD']:.4f}")

            d.update({
                "index": i,
                "Rescore_Total": f"{rescore:.4f}",
                "Weighted_LGFE": f"{wtd_LGFE:.4f}",
                "Weighted_LE": f"{wtd_LE:.4f}",
                "Weighted_RMSD": f"{wtd_RMSD:.4f}",
                "LGFE": f"{d['LGFE']:.4f}",
                "LE": f"{d['LE']:.4f}",
                "RMSD": f"{d['RMSD']:.4f}"
            })
            all_data.append(d)

        # best pose per compound (by minimal Rescore_Total)
        best_by_compound = {}
        grouped = defaultdict(list)
        for d in all_data:
            cid = d["Compound_ID"].rsplit("_", 1)[0]
            grouped[cid].append(d)

        for cid, poses in grouped.items():
            best_pose = min(poses, key=lambda x: float(x["Rescore_Total"]))
            best_by_compound[cid] = best_pose

        return all_data, list(best_by_compound.values()), None

    def write_csv(self, entries, filename):
        fieldnames = ["Compound_ID", "SMILES", "LGFE", "LE", "RMSD",
                      "Weighted_LGFE", "Weighted_LE", "Weighted_RMSD", "Rescore_Total"]
        with open(filename, 'w', newline='') as outf:
            writer = csv.DictWriter(outf, fieldnames=fieldnames)
            writer.writeheader()
            for e in entries:
                writer.writerow({k: e.get(k, "") for k in fieldnames})

    def write_sdf(self, entries, filename):
        writer = Chem.SDWriter(filename)
        for e in entries:
            mol = self.sdf_supplier[e["index"]]
            for prop in ["LGFE", "LE", "RMSD", "Weighted_LGFE",
                         "Weighted_LE", "Weighted_RMSD", "Rescore_Total"]:
                mol.SetProp(prop, e[prop])
            writer.write(mol)
        writer.close()

    def process_duplicates(self, data):
        # group by base compound ID (strip trailing _pose#)
        groups = {}
        for d in data:
            base = d["Compound_ID"].rsplit('_', 1)[0]
            groups.setdefault(base, []).append(d)
        deduped = []
        for base, lst in groups.items():
            if len(lst) > npose:
                lst = sorted(lst, key=lambda x: float(x["Rescore_Total"]))[:npose]
            deduped.extend(lst)
        return deduped

if __name__ == "__main__":
    # merge SDFs
    minconfdir      = os.path.join(silcsmcdir, "minconfpdb")
    merged_sdf_path = f"silcsmc_{name}_combined.sdf"
    silcs = SILCSMCProcessor(minconfdir)
    lgfe_csv_path = os.path.join(silcsmcdir, "lgfe.csv")
    silcs.merge_sdfs(merged_sdf_path, lgfe_csv_path)

    if not os.path.exists(merged_sdf_path):
        sys.exit(f"Error: '{merged_sdf_path}' not found.")

    # load, rescore, and write outputs
    proc = SDFProcessor(merged_sdf_path)
    all_data, best_data, _ = proc.set_properties()
    all_data = proc.process_duplicates(all_data)

    # write full table and best-poses
    proc.write_sdf(all_data,  f"{name}_all_poses_rescored_scaled.sdf")
    proc.write_csv(all_data,  f"{name}_all_poses_rescored_scaled.csv")
    proc.write_csv(best_data, f"{name}_best_poses_rescored_scaled.csv")
