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
  selection.method = "moransi",
  nfeatures = length(VariableFeatures(obj)))

svg_info <- SVFInfo(obj, method = "moransi")
svg_info$MoransI_p.adj <- p.adjust(svg_info$MoransI_p.value, method = "BH")

sum(!is.na(svg_info$MoransI_p.value))
#20
#FindSpatiallyVariableFeatures silently failed to compute Moran's I for all but 20 of the 
#3000 requested genes - a reproducible behavior matching documented Seurat bug reportson 
#high-resolution spatial data (issue#9087 and discussion#9608 - and we report results for 
#the 20 genes it actually computed
# Add gene names as a column, sort by adjusted P-value (most significant SVGs first)
svg_info$gene <- rownames(svg_info)
svg_info <- svg_info[order(svg_info$MoransI_p.adj), 
                      c("gene", "MoransI_observed", "MoransI_p.value", "MoransI_p.adj")]

write.csv(svg_info,
          "/ddnlus/r2762/spatial/arabidopsis_visium/results/SVG_moransi_results.csv",
          row.names = FALSE)

#top1 is EXT3 which is a cell wall glycoprotein and has a role in root hair elongation and 
# check EXT3 specifically
svg_info["EXT3", ]
#gene MoransI_observed MoransI_p.value MoransI_p.adj
 EXT3        0.7907884               0             0
#BH p-value adjustment was applied and confirmed that adjusted p-value is significant

## ---- 2. Overlap between top SVGs and cluster markers ---------------
markers <- read.csv(file.path(results_dir, "cluster_markers.csv"), row.names = 1)
svg_also_markers <- intersect(top_svg, markers$gene)
cat("SVGs that are also cluster markers:", length(svg_also_markers), "\n")
print(svg_also_markers)

svg_marker_details <- markers %>%
  filter(gene %in% svg_also_markers) %>%
  select(gene, cluster, avg_log2FC, pct.1, pct.2, p_val_adj)

write.csv(svg_marker_details, file.path(results_dir, "SVG_and_cluster_markers.csv"), row.names = FALSE)
