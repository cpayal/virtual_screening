Pose Rescoring and Best-Pose Selection (Step 02)
===============================================

Usage
-----
python step02-analysis_rescoring.py <project_dir> <output_name>

Inputs in <project_dir>
- minconfpdb/*.sdf
- lgfe.csv with columns: Ligand, LGFE(kcal/mol)

Outputs
-------
- <output_name>_all_poses_rescored_scaled.sdf
- <output_name>_all_poses_rescored_scaled.csv
- <output_name>_best_poses_rescored_scaled.csv

Notes
-----
- LGFE threshold used to pre-filter ligands (≤ -7 kcal/mol).
- Rescore combines LGFE, LE, RMSD (default 1:1:1).
- Best pose per compound chosen by min Rescore_Total.
