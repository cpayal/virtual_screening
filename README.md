Virtual Screening: Pose Rescoring → Filtering/Clustering
=======================================================

This repo contains a 2-step pipeline:

Step 02: Rescoring and Best-Pose Selection
- Input: <project_dir>/minconfpdb/*.sdf, <project_dir>/lgfe.csv
- Output: <name>_all_poses_rescored_scaled.sdf, <name>_best_poses_rescored_scaled.csv

Step 03: Pose Convergence, PhysChem Filters, Clustering, and Substructure 3D Tags
- Input: <name>_all_poses_rescored_scaled.sdf
- Output: <prefix>_filtered_bestposes.sdf/.csv, <prefix>_clusters.csv, <prefix>_cluster_summary.csv

Quickstart
----------
# Step 02
python step02_rescore/step02-analysis_rescoring.py <project_dir> <name>

# Step 03
python step03_filter_cluster/step03-filter_and_cluster.py <name>_all_poses_rescored_scaled.sdf <out_prefix>

Environment
-----------
- Python 3.8+
- RDKit, NumPy, Pandas
(optional) Create with conda:
    conda env create -f env/environment.yml
    conda activate vs-env
