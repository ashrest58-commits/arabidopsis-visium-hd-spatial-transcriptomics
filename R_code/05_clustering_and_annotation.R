## ============================================================
## 05_clustering_and_annotation.R
## Unsupervised clustering, marker discovery, and assignment
## of biological tissue-domain names via TAIR cross-reference.
## ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)

results_dir <- "/path/to/spatial/arabidopsis_visium/results"
obj <- readRDS(file.path(results_dir, "arabidopsis_visium_normalized.rds"))

## ---- 1. Variable features, scaling, PCA --------------------------
## vst (not LeverageScore) is used because this object is small
## (1,516 bins) after QC filtering -- LeverageScore subsampling is
## built for much larger bin counts than this dataset has.
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 50)

png(file.path(results_dir, "elbow_plot.png"), width = 800, height = 600)
print(ElbowPlot(obj, ndims = 50))
dev.off()
## Variance flattens around PC 9-10; 15 PCs used to be conservative.

## ---- 2. Neighbors, clustering, UMAP -------------------------------
obj <- FindNeighbors(obj, dims = 1:15)
obj <- FindClusters(obj, resolution = 0.5)
obj <- RunUMAP(obj, dims = 1:15)

table(Idents(obj))

my_colors <- c("0" = "#E64B35", "1" = "#4DBBD5", "2" = "#006400", "3" = "#F39C12",
               "4" = "#800000", "5" = "#3C5488", "6" = "#800080")

png(file.path(results_dir, "UMAP_clusters.png"), width = 1600, height = 1200, res = 150)
print(DimPlot(obj, reduction = "umap", label = TRUE, label.size = 5, pt.size = 2,
               repel = TRUE, cols = my_colors) + ggtitle("Arabidopsis seedling - UMAP"))
dev.off()

png(file.path(results_dir, "spatial_clusters.png"), width = 1200, height = 800, res = 100)
print(SpatialDimPlot(obj, label = FALSE, pt.size.factor = 10, cols = my_colors) +
        ggtitle("Arabidopsis seedling - Spatial clusters"))
dev.off()

saveRDS(obj, file.path(results_dir, "arabidopsis_visium_clustered.rds"))

## ---- 3. Marker discovery -------------------------------------------
markers <- FindAllMarkers(obj, assay = "Spatial.008um", only.pos = TRUE)
write.csv(markers, file.path(results_dir, "cluster_markers.csv"), row.names = TRUE)

top5 <- markers %>% group_by(cluster) %>% filter(avg_log2FC > 1) %>% slice_head(n = 5) %>% ungroup()
top50 <- markers %>% group_by(cluster) %>% filter(avg_log2FC > 1) %>% slice_head(n = 50) %>% ungroup()
write.csv(top50, file.path(results_dir, "top50_markers_per_cluster.csv"), row.names = FALSE)

obj <- ScaleData(obj, assay = "Spatial.008um", features = top5$gene)
png(file.path(results_dir, "cluster_heatmap1.png"), width = 2400, height = 1600, res = 150)
print(
  DoHeatmap(obj, assay = "Spatial.008um", features = top5$gene, size = 4) +
    scale_fill_gradientn(colors = c("blue", "white", "red"), name = "Z-score")
)
dev.off()

## ---- 4. Cluster annotation ------------------------------------------
## Names assigned by cross-referencing top5/top50 marker genes per
## cluster against the TAIR expression atlas, one cluster at a time.
new.cluster.ids <- c(
  "0" = "Cotyledon",
  "1" = "Germinating seed",
  "2" = "Root elongation zone",
  "3" = "Root",
  "4" = "Meristematic tissue",
  "5" = "hypocotyl",
  "6" = "SAM"
)
obj <- RenameIdents(obj, new.cluster.ids)

colors <- c("Cotyledon" = "#E64B35", "Germinating seed" = "#4DBBD5",
            "Root elongation zone" = "#006400", "Root" = "#3C5488",
            "Meristematic tissue" = "#F39C12", "hypocotyl" = "#800000", "SAM" = "#800080")

png(file.path(results_dir, "UMAP_annotated1_nolabel.png"), width = 1200, height = 800)
print(DimPlot(obj, reduction = "umap", label = FALSE, pt.size = 2, repel = TRUE, cols = colors) +
        ggtitle("Arabidopsis seedling - Annotated clusters"))
dev.off()

png(file.path(results_dir, "spatial1_annotated.png"), width = 1200, height = 800, res = 100)
print(SpatialDimPlot(obj, label = FALSE, pt.size.factor = 10, cols = colors) +
        ggtitle("Arabidopsis seedling - Spatial annotation"))
dev.off()

saveRDS(obj, file.path(results_dir, "arabidopsis_visium_annotated.rds"))
