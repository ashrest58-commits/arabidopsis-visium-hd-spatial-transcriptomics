## ============================================================
## 01_setup_and_load.R
## Download the public Visium HD Arabidopsis dataset and load
## the 8 um-binned matrix into a Seurat object.
## ============================================================

## ---- 0. Paths (UPDATE THESE) --------------------------------
data_dir    <- "/path/to/spatial/arabidopsis_visium/data"
results_dir <- "/path/to/spatial/arabidopsis_visium/results"
dir.create(data_dir,    recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

## ---- 1. Download the dataset (run once, in bash/terminal) ---
## wget -O binned_outputs.tar.gz \
##   "https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_3prime_Arabidopsis/Visium_HD_3prime_Arabidopsis_binned_outputs.tar.gz"
##
## wget -O spatial_outputs.tar.gz \
##   "https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_3prime_Arabidopsis/Visium_HD_3prime_Arabidopsis_spatial.tar.gz"
##
## wget -O metrics_summary.csv \
##   "https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_3prime_Arabidopsis/Visium_HD_3prime_Arabidopsis_metrics_summary.csv"
##
## tar -xzf binned_outputs.tar.gz
## tar -xzf spatial_outputs.tar.gz

## ---- 2. Load required packages -------------------------------
library(Seurat)
library(hdf5r)   # required to read the Visium HD .h5 matrix
library(arrow)    # required to read tissue_positions.parquet
library(ggplot2)

## ---- 3. Load the 8 um-binned Visium HD object -----------------
## Three bin sizes are available in the same folder: 2 um, 8 um, 16 um.
## 2 um bins are sparse (often 0 UMIs); 16 um loses spatial resolution;
## 8 um is the commonly used compromise.
obj <- Load10X_Spatial(data.dir = data_dir, bin.size = 8)

obj
## Expected:
## An object of class Seurat
## 31109 features across 7356 samples within 1 assay
## Active assay: Spatial.008um (31109 features, 0 variable features)
## 1 spatial field of view present: slice1.008um

saveRDS(obj, file.path(results_dir, "arabidopsis_visium_raw.rds"))
