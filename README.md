## Virtual Screening: Pose Rescoring → Filtering/Clustering

### Purpose
This repository is a collection of scripts to triage millions of virtual screening hits down to a handful of best hits with the goal of obtaining a high hit rate.

### Philosophy
Our hits are only as good as how we select them! Ultra-large virtual screening campaigns have become very popular, driven by exceptional algorithms. However, the hit selection criteria still need to be solidified, since these algorithms often cannot differentiate beyond simple scores (e.g., free energy or docking). As a computational chemist, my goal is to prioritize hits that provide the highest confidence of experimental success.

This workflow addresses that by introducing multiple layers of decision-making beyond raw scores:

- #### Rescoring and normalization :
  Docking or LGFE scores are scaled and re-evaluated to reduce bias from different scoring functions.
- #### Best-pose selection :
  Among multiple poses per ligand, the most consistent and representative conformations are chosen.
- #### Pose convergence checks:
  Ensures that clusters of ligands bind in a reproducible way, filtering out unstable or artifactual poses.
- #### Physicochemical filters:
  Removes compounds unlikely to succeed due to poor drug-like properties.
- #### Clustering and substructure tagging:
  Promotes chemical diversity and helps retain structurally meaningful representatives.

By combining these steps, the workflow emphasizes convergence, drug-likeness, and diversity, giving higher confidence that the selected hits will be experimentally validated.

### Why this matters
Traditional large-scale virtual screens often report hit rates below 1%, since most compounds are selected solely on raw docking scores. By layering rescoring, convergence analysis, property filters, and diversity-driven clustering, this workflow consistently enriches for compounds that survive experimental validation. In practice, this approach has achieved hit rates of around 10%, a substantial improvement that saves both synthesis effort and assay resources.

### Pipeline Overview

##### Step 02: Rescoring and Best-Pose Selection
Input: <project_dir>/minconfpdb/*.sdf, <project_dir>/lgfe.csv
Output: <name>_all_poses_rescored_scaled.sdf, <name>_best_poses_rescored_scaled.csv

##### Step 03: Pose Convergence, PhysChem Filters, Clustering, and Substructure 3D Tags
Input: <name>_all_poses_rescored_scaled.sdf
Output: <prefix>_filtered_bestposes.sdf/.csv, <prefix>_clusters.csv, <prefix>_cluster_summary.csv

#### Quickstart
##### Step 02
`python step02_rescore/step02-analysis_rescoring.py <project_dir> <name>`

##### Step 03
`python step03_filter_cluster/step03-filter_and_cluster.py <name>_all_poses_rescored_scaled.sdf <out_prefix>`

##### Environment
`Python 3.8+
RDKit
NumPy
Pandas`

(Optional) create with conda:

`conda env create -f env/environment.yml`
`conda activate vs-env`
