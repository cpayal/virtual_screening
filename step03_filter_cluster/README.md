Filtering, PhysChem Gates, Clustering, and Substructure 3D Tags (Step 03)
========================================================================

Usage
-----
python step03-filter_and_cluster.py <rescored_sdf> <out_prefix>

Input
-----
- <rescored_sdf> = output from Step 02 (e.g., <name>_all_poses_rescored_scaled.sdf)

Pipeline
--------
1) Pose convergence vs best pose:
   - Align to best-scoring pose on heavy atoms; keep poses with RMSD ≤ 3 Å.
   - Drop molecule if it has < 3 converged poses.
2) PhysChem window:
   - MW 300–550, TPSA 40–140, RotB ≤ 10, HBD ≤ 5, (HBD+HBA) ≤ 12
3) Final LE filter:
   - Keep only molecules with LE ≤ -0.3
4) Clustering (chemotype):
   - ECFP (Morgan, r=2, 2048 bits) + Butina, Tanimoto 0.6
5) Substructure 3D checks:
   - Intra-molecule: tag "Substructure Convergence" across poses vs best (≤ 3 Å)
   - Cluster-level: tag "Cluster Substructure Alignment" across best poses and
     store "CSA_MaxRMSD"

Outputs
-------
- <out_prefix>_filtered_bestposes.sdf
- <out_prefix>_filtered_bestposes.csv
- <out_prefix>_clusters.csv
- <out_prefix>_cluster_summary.csv
