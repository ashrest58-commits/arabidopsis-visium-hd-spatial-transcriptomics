## ============================================================
## 02_environment_setup.sh
## Conda environment setup for packages not bundled with a
## default Seurat installation. Run once before R_code scripts.
## ============================================================

conda activate seurat_env

## Check which packages are already present
conda list | grep -iE "seurat|biocmanager|nnsvg|spatialexperiment|ggplot2|bioconductor"

## --- Bioconductor package used for pathway analysis ---
conda install -c bioconda -c conda-forge bioconductor-clusterprofiler
## org.At.tair.db (Arabidopsis TAIR gene ID / GO annotation database)
## install via Bioconductor's R interface if not available through conda:
##   BiocManager::install("org.At.tair.db")

## --- System-level R packages required to read Visium HD output ---
conda install -c conda-forge r-hdf5r     # reads the Visium HD .h5 count matrix
conda install -c conda-forge r-arrow     # reads tissue_positions.parquet
conda install -c conda-forge r-ape       # required by FindSpatiallyVariableFeatures(method = "moransi")

## --- Plotting helpers used in DE volcano plots and GSEA running-score plots ---
## install.packages("ggrepel")          # label repelling on volcano plots
## BiocManager::install("enrichplot")   # gseaplot2() for GSEA running-score plots

## Verify installation
conda list | grep -iE "seurat|biocmanager|nnsvg|spatialexperiment|ggplot2|bioconductor|hdf5r|arrow|ape"

## For long-running steps (clustering, GSEA permutation testing), submit
## an interactive job to a compute node rather than running on a login node:
##   qsub -I -lncpus=16 -lmem=48gb -lwalltime=24:00:00
##   conda activate seurat_env
##   R
