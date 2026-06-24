## ============================================================
## 06_SVG_detection.R
## Spatially variable gene (SVG) detection via Moran's I,
## independent of cluster identity, and overlap with cluster markers.
## ============================================================

library(Seurat)
library(ape)   # required by FindSpatiallyVariableFeatures(method = "moransi")
library(dplyr)

results_dir <- "/path/to/spatial/arabidopsis_visium/results"
obj <- readRDS(file.path(results_dir, "arabidopsis_visium_annotated.rds"))
DefaultAssay(obj) <- "Spatial.008um"

## ---- 1. Moran's I-based SVG detection -----------------------------
## Moran's I ranges from -1 (perfect dispersion of dissimilar values)
## through 0 (random / no spatial autocorrelation) to +1 (perfect
## clustering of similar values).
obj <- FindSpatiallyVariableFeatures(obj,
  assay = "Spatial.008um",
  features = VariableFeatures(obj),
  selection.method = "moransi")

top_svg <- SpatiallyVariableFeatures(obj, method = "moransi")
svg_info <- SVFInfo(obj, method = "moransi", status = TRUE)
write.csv(svg_info, file.path(results_dir, "spatially_variable_genes.csv"), row.names = TRUE)

head(top_svg, 20)

png(file.path(results_dir, "top6_SVGs.png"), width = 2400, height = 1600, res = 150)
print(SpatialFeaturePlot(obj, features = head(top_svg, 10), pt.size.factor = 8, ncol = 3))
dev.off()

## ---- 2. Overlap between top SVGs and cluster markers ---------------
markers <- read.csv(file.path(results_dir, "cluster_markers.csv"), row.names = 1)
svg_also_markers <- intersect(top_svg, markers$gene)
cat("SVGs that are also cluster markers:", length(svg_also_markers), "\n")
print(svg_also_markers)

svg_marker_details <- markers %>%
  filter(gene %in% svg_also_markers) %>%
  select(gene, cluster, avg_log2FC, pct.1, pct.2, p_val_adj)

write.csv(svg_marker_details, file.path(results_dir, "SVG_and_cluster_markers.csv"), row.names = FALSE)
