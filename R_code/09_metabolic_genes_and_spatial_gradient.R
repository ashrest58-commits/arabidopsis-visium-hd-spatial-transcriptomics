## ============================================================
## 09_metabolic_genes_and_spatial_gradient.R
## (a) Composite spatial plot of one representative marker per
##     tissue domain.
## (b) Spatial gradient modeling: separate physically distinct
##     tissue pieces, then fit GAMs of gene expression against
##     distance from the SAM centroid.
## ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(mgcv)

results_dir <- "/path/to/spatial/arabidopsis_visium/results"
obj <- readRDS(file.path(results_dir, "arabidopsis_visium_annotated.rds"))

## ============================================================
## PART A: Key metabolic gene spatial mapping
## ============================================================

key_markers <- c(
  "CAB1",     # Cotyledon
  "ICL",      # Germinating seed
  "EXPA9",    # Root elongation zone
  "EXT3",     # Root
  "GRP4",     # Meristematic tissue
  "XYP2",     # Hypocotyl
  "CYP78A5")  # SAM

present <- key_markers[key_markers %in% rownames(obj)]

png(file.path(results_dir, "spatial_key_markers_all_tissues.png"),
    width = 2800, height = 1600, res = 150)
print(SpatialFeaturePlot(obj, features = present, pt.size.factor = 10, ncol = 4))
dev.off()

## ============================================================
## PART B: Spatial gradient modeling
## ============================================================

coords <- GetTissueCoordinates(obj)   # columns: x, y, cell

gradient_genes <- c(
  "CYP78A5", "AFO", "ANT",     # SAM (shoot meristem)
  "GRP4", "NOP10",             # Meristematic (root meristem)
  "EXPA9", "EIR1",             # Root elongation zone
  "EXT3", "AGP14", "XYP2",     # Mature root
  "CAB1", "RBCS1A")            # Cotyledon (shoot maturation)
present_genes <- gradient_genes[gradient_genes %in% rownames(obj)]

expr_data <- FetchData(obj, vars = present_genes)
expr_data$x <- coords[rownames(expr_data), "x"]
expr_data$y <- coords[rownames(expr_data), "y"]
expr_data$cluster <- Idents(obj)[rownames(expr_data)]

## ---- Problem: this Visium HD section contains several physically
## separate tissue pieces. A single whole-slide distance-from-SAM
## calculation would incorrectly mix unrelated tissue fragments.
## Solution: cluster bins by PHYSICAL (x, y) proximity first,
## independent of transcriptional cluster identity.
set.seed(123)
piece_clusters <- kmeans(expr_data[, c("x", "y")], centers = 8)
expr_data$piece <- as.factor(piece_clusters$cluster)

table(expr_data$piece, expr_data$cluster)
## "Piece 1" contains the complete shoot developmental axis:
## Cotyledon (80) + Germinating seed (44) + Hypocotyl (51) + SAM (38)

piece1_data <- expr_data %>% filter(piece == "1")
cat("Piece 1 bins:", nrow(piece1_data), "\n")

png(file.path(results_dir, "piece1_check.png"), width = 800, height = 1000, res = 150)
print(
  ggplot(piece1_data, aes(x = x, y = y, color = cluster)) +
    geom_point(size = 2) + coord_fixed() + scale_y_reverse() +
    theme_minimal() + ggtitle("Piece 1: Cotyledon - Hypocotyl - SAM axis")
)
dev.off()

## Use the SAM centroid (within piece 1) as the developmental
## organizer reference point.
sam_centroid <- piece1_data %>% filter(cluster == "SAM") %>%
  summarise(x_ref = mean(x), y_ref = mean(y))

piece1_data <- piece1_data %>%
  mutate(dist_from_SAM = sqrt((x - sam_centroid$x_ref)^2 + (y - sam_centroid$y_ref)^2))

## Sanity check: distance from SAM should increase along the shoot axis
piece1_data %>%
  group_by(cluster) %>%
  summarise(mean_dist = mean(dist_from_SAM)) %>%
  arrange(mean_dist)
## Expected order: SAM (69.2) < Germinating seed (199) < Cotyledon (413) < Hypocotyl (542)

## ---- Fit GAM models for key genes along the SAM-distance axis ------
gradient_genes_p1 <- c("CYP78A5", "AFO", "ANT", "CAB1", "RBCS1A")
present_p1 <- gradient_genes_p1[gradient_genes_p1 %in% colnames(piece1_data)]

for (gene in present_p1) {
  formula_str <- as.formula(paste0(gene, " ~ s(dist_from_SAM, k = 4)"))
  gam_model <- gam(formula_str, data = piece1_data)
  print(summary(gam_model))

  p <- ggplot(piece1_data, aes(x = dist_from_SAM, y = .data[[gene]])) +
    geom_point(aes(color = cluster), alpha = 0.6) +
    geom_smooth(method = "gam", formula = y ~ s(x, k = 4), color = "black") +
    labs(title = paste(gene, "expression vs distance from SAM"),
         x = "Distance from SAM centroid", y = paste(gene, "expression")) +
    theme_minimal()

  png(file.path(results_dir, paste0("gradient_", gene, ".png")), width = 1000, height = 700, res = 150)
  print(p)
  dev.off()
}

## Findings:
## - AFO: sharply restricted to ~300-500 spatial units of the SAM
##   centroid, consistent with a localized shoot meristem identity gene.
## - ANT: detected in 22% of SAM bins, 0% of Cotyledon bins; significant
##   DE (log2FC = 11.6, padj = 8.1e-26) AND a significant, modest spatial
##   gradient with distance from SAM (GAM p = 0.00026, 9.6% deviance explained).
## - CAB1: bimodal on/off pattern in every cluster, NO significant
##   relationship to distance from SAM (p = 0.135) -- consistent with
##   light exposure / chloroplast content rather than a developmental
##   gradient tied to distance from a shoot organizer.

saveRDS(obj, file.path(results_dir, "arabidopsis_visium_final.rds"))
