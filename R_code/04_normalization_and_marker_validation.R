## ============================================================
## 04_normalization_and_marker_validation.R
## Log-normalization, followed by an independent sanity check
## using known TAIR-documented tissue markers, before any
## unsupervised analysis is trusted.
## ============================================================

library(Seurat)
library(ggplot2)

results_dir <- "/path/to/spatial/arabidopsis_visium/results"
obj <- readRDS(file.path(results_dir, "arabidopsis_visium_filtered.rds"))

## ---- 1. Log-normalization ---------------------------------------
DefaultAssay(obj) <- "Spatial.008um"
obj <- NormalizeData(obj)   # LogNormalize, scale.factor = 10,000
saveRDS(obj, file.path(results_dir, "arabidopsis_visium_normalized.rds"))

## ---- 2. Known marker gene validation -----------------------------
## Several candidate marker pairs were tried before settling on
## PHB / ERL2 (shoot meristem markers used in 10x's own
## demonstration of this dataset). Their spatial pattern should
## match the published TAIR pattern BEFORE trusting any clustering.
shoot_meristem_genes <- c("PHB", "ERL2")
present <- shoot_meristem_genes[shoot_meristem_genes %in% rownames(obj)]
print(present)

png(file.path(results_dir, "PHB_ERL2_markers.png"), width = 1200, height = 500)
p1 <- SpatialFeaturePlot(obj, features = "PHB",  pt.size.factor = 1.2) + ggtitle("PHABULOSA")
p2 <- SpatialFeaturePlot(obj, features = "ERL2", pt.size.factor = 1.2) + ggtitle("ERECTA-LIKE 2")
print(p1 | p2)
dev.off()

## Other candidates tried along the way (kept here for reproducibility,
## not used in the final figure): HY5/AP2, STM/SCR
## "HY5" %in% rownames(obj); "AP2" %in% rownames(obj)
## "STM" %in% rownames(obj); "SCR" %in% rownames(obj)
