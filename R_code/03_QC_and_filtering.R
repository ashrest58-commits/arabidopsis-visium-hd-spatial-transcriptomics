## ============================================================
## 03_QC_and_filtering.R
## Quality control, on-tissue and low-complexity bin filtering,
## and removal of zero-expression genes.
## ============================================================

library(Seurat)
library(ggplot2)

results_dir <- "/path/to/spatial/arabidopsis_visium/results"
obj <- readRDS(file.path(results_dir, "arabidopsis_visium_raw.rds"))

## ---- 1. Pre-filtering QC ---------------------------------------
summary(obj@meta.data$nCount_Spatial.008um)
summary(obj@meta.data$nFeature_Spatial.008um)

vln.plot   <- VlnPlot(obj, features = "nCount_Spatial.008um", pt.size = 0) + NoLegend()
count.plot <- SpatialFeaturePlot(obj, features = "nCount_Spatial.008um", pt.size.factor = 1.2) +
  theme(legend.position = "right")

png(file.path(results_dir, "QC_UMI_plots.png"), width = 1200, height = 500)
print(vln.plot | count.plot)
dev.off()

## ---- 2. Remove off-tissue bins ----------------------------------
table(obj@meta.data$nCount_Spatial.008um > 0)
obj <- subset(obj, nCount_Spatial.008um > 0)

## ---- 3. Remove low-complexity bins -------------------------------
## NOTE: mitochondrial-content filtering is intentionally NOT applied
## here. Removing high-mitochondrial bins/cells can selectively
## deplete biologically meaningful, metabolically active regions
## rather than removing only technical artifacts (PMC11983838) --
## a risk this case study's metabolic-role question cannot afford.
obj <- subset(obj,
  nCount_Spatial.008um   > 100 &
  nFeature_Spatial.008um > 50
)
dim(obj)  # expect: 31109 features across 1516 bins

## ---- 4. Remove genes with zero expression across remaining bins --
counts_matrix <- GetAssayData(obj, assay = "Spatial.008um", layer = "counts")
genes_keep    <- rowSums(counts_matrix) > 0
obj <- obj[genes_keep, ]
dim(obj)  # expect: 17852 features across 1516 bins

## ---- 5. Post-filtering QC visual check ---------------------------
png(file.path(results_dir, "QC_after_filtering.png"), width = 800, height = 600)
print(
  SpatialFeaturePlot(obj, features = "nCount_Spatial.008um", pt.size.factor = 1.2) +
    theme(legend.position = "right")
)
dev.off()

saveRDS(obj, file.path(results_dir, "arabidopsis_visium_filtered.rds"))
