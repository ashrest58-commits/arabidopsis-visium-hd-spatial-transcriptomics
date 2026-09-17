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
#3000 requested genes - a reproducible behavior consistent with documented Seurat bug reports on 
#high-resolution spatial data (issue#9087 and discussion#9608 - 

### ---- Clean-coordinate diagnostic ------------------------------------
coords <- GetTissueCoordinates(obj)
str(coords)
#'data.frame':   1516 obs. of  3 variables:
 $ x   : num  5179 4918 6898 4918 10491 ...
 $ y   : num  8635 8428 15048 8486 15054 ...
 $ cell: chr  "s_008um_00456_00333-1" "s_008um_00463_00324-1" "s_008um_00237_00394-1" "s_008um_00461_00324-1" ...

spatial_loc <- as.matrix(coords[, c("x", "y")])
rownames(spatial_loc) <- coords$cell

## Sanity check: rownames must exactly match the assay's cell order.
stopifnot(all(rownames(spatial_loc) == Cells(obj)))
cat("Coordinate/cell alignment check passed.\n")
#check passed

assay_obj <- obj[["Spatial.008um"]]
assay_obj <- FindSpatiallyVariableFeatures(
  assay_obj,
  spatial.location = spatial_loc,
  selection.method = "moransi",
  features = VariableFeatures(obj),
  nfeatures = length(VariableFeatures(obj))
)

svg_info_clean <- SVFInfo(assay_obj, method = "moransi")
sum(!is.na(svg_info_clean$MoransI_p.value))
#20

# we report results for 
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

##According to TAIR, EXT3 is expressed in root, but XYP2 is also highly expressed in root nad 
#has positive Moran’s I index and significant adjusted P value

## Do the EXT3/AGP31 positive bins actually sit on your "Root" cluster,
# or somewhere else in the tissue?
root_cells <- WhichCells(obj, idents = "Root")

counts <- GetAssayData(obj, assay = "Spatial.008um", layer = "counts")
for (g in c("EXT3", "AGP31", "XYP2")) {
  pos_cells <- colnames(counts)[counts[g, ] > 0]
  overlap <- length(intersect(pos_cells, root_cells))
  cat(g, ": ", length(pos_cells), " positive bins total, ",
      overlap, " of them (", round(100*overlap/length(pos_cells), 1),
      "%) fall within the Root cluster\n", sep = "")
} 

#EXT3: 255 positive bins total, 119 of them (46.7%) fall within the Root cluster
AGP31: 332 positive bins total, 52 of them (15.7%) fall within the Root cluster
XYP2: 41 positive bins total, 34 of them (82.9%) fall within the Root cluster


#even though EXT3 has highest positive bins, only 46% fall in root cluster whereas 83% of XYP2 
#fall in root cluster so will go with XYP2 as representative


## ---- 2. Overlap between top SVGs and cluster markers ---------------
markers <- read.csv(file.path(results_dir, "cluster_markers.csv"), row.names = 1)
svg_also_markers <- intersect(top_svg, markers$gene)
cat("SVGs that are also cluster markers:", length(svg_also_markers), "\n")
print(svg_also_markers)

svg_marker_details <- markers %>%
  filter(gene %in% svg_also_markers) %>%
  select(gene, cluster, avg_log2FC, pct.1, pct.2, p_val_adj)

write.csv(svg_marker_details, file.path(results_dir, "SVG_and_cluster_markers.csv"), row.names = FALSE)
