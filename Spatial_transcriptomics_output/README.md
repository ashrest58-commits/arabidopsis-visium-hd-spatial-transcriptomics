# Spatial_transcriptomics_output/

This folder is where all RDS checkpoints, CSV tables, and figures from the
`R_code/` scripts are saved (see `results_dir` at the top of each script).

It is intentionally left empty in this repository. Outputs are not version-
controlled here because:

- Seurat `.rds` objects and raw count matrices can be large.
- They are fully reproducible by re-running the numbered scripts in
  `R_code/` against the publicly available dataset (see the main README
  for the download commands).

If you want to share specific output files (e.g., final figures or marker
tables) via GitHub, add them here individually rather than committing the
full output directory, or use Git LFS / Zenodo / institutional storage for
larger files.
