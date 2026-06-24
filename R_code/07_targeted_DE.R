## ============================================================
## 07_targeted_DE.R
## Eight pairwise differential expression comparisons between
## specific tissue-domain pairs (Wilcoxon rank-sum test).
## ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)

results_dir <- "/path/to/spatial/arabidopsis_visium/results"
obj <- readRDS(file.path(results_dir, "arabidopsis_visium_annotated.rds"))
table(Idents(obj))

run_de <- function(obj, ident1, ident2, label) {
  de <- FindMarkers(obj,
    assay = "Spatial.008um",
    ident.1 = ident1,
    ident.2 = ident2,
    min.pct = 0.10,
    logfc.threshold = 0.25)
  cat(label, "DEGs:", nrow(de), "\n")
  write.csv(de, file.path(results_dir, paste0("DE_", label, ".csv")), row.names = TRUE)
  de
}

## ---- 8 targeted comparisons -----------------------------------------
de_cot_germ <- run_de(obj, "Cotyledon", "Germinating seed",     "cotyledon_vs_germinating")  # 330 DEGs
de_cot_hyp  <- run_de(obj, "Cotyledon", "hypocotyl",            "cotyledon_vs_hypocotyl")    # 212 DEGs
de_root_hyp <- run_de(obj, "Root",      "hypocotyl",            "Root_vs_hypocotyl")         # 1304 DEGs
de_root_rez <- run_de(obj, "Root",      "Root elongation zone", "Root_vs_REZ")               # 1652 DEGs
de_root_mt  <- run_de(obj, "Root",      "Meristematic tissue",  "Root_vs_MT")                # 5462 DEGs
de_rez_mt   <- run_de(obj, "Root elongation zone", "Meristematic tissue", "REZ_vs_MT")        # 4974 DEGs
de_sam_mt   <- run_de(obj, "SAM",       "Meristematic tissue",  "SAM_vs_MT")                 # 5354 DEGs
de_sam_cot  <- run_de(obj, "SAM",       "Cotyledon",            "SAM_vs_cot")                # 861 DEGs

## ---- Volcano plot helper ---------------------------------------------
make_volcano <- function(de, label_hi, label_lo, color_hi, color_lo, title, out_file) {
  de$significance <- "Not significant"
  de$significance[de$p_val_adj < 0.05 & de$avg_log2FC >  1] <- label_hi
  de$significance[de$p_val_adj < 0.05 & de$avg_log2FC < -1] <- label_lo

  top_hi <- de %>% filter(significance == label_hi) %>% arrange(desc(avg_log2FC)) %>% head(5)
  top_lo <- de %>% filter(significance == label_lo) %>% arrange(avg_log2FC) %>% head(5)
  de$label <- ""
  de$label[rownames(de) %in% c(rownames(top_hi), rownames(top_lo))] <-
    rownames(de)[rownames(de) %in% c(rownames(top_hi), rownames(top_lo))]

  cols <- setNames(c(color_hi, color_lo, "grey70"), c(label_hi, label_lo, "Not significant"))

  p <- ggplot(de, aes(avg_log2FC, -log10(p_val_adj), color = significance, label = label)) +
    geom_point(alpha = 0.6, size = 1.5) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 20) +
    scale_color_manual(values = cols) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    labs(x = "Average log2 Fold Change", y = "-log10 adjusted p-value",
         title = title, color = "Expression") +
    theme_classic()

  png(file.path(results_dir, out_file), width = 1200, height = 800, res = 150)
  print(p)
  dev.off()
}

make_volcano(de_cot_germ, "Higher in Cotyledon", "Higher in Germinating seed",
             "#E64B35", "#4DBBD5", "Cotyledon vs Germinating seed", "volcano_cot_vs_germ.png")
make_volcano(de_sam_cot,  "Higher in SAM", "Higher in Cotyledon",
             "#800080", "#E64B35", "SAM vs Cotyledon", "volcano_SAM_vs_cot.png")
make_volcano(de_cot_hyp,  "Higher in Cotyledon", "Higher in hypocotyl",
             "#800080", "#E64B35", "Cotyledon vs Hypocotyl", "volcano_cot_vs_hyp.png")
make_volcano(de_root_hyp, "Higher in Root", "Higher in Hypocotyl",
             "#3C5488", "#800000", "Root vs Hypocotyl", "volcano_root_vs_hyp.png")

## ---- Spatial feature plots of top genes per comparison ---------------
## Example (Cotyledon vs Germinating seed); repeat pattern for the others.
top_genes_cot_germ <- c("CAB1", "LEJ2", "FER1", "AT5G56670", "AT4G15000", "AT5G56710")
present <- top_genes_cot_germ[top_genes_cot_germ %in% rownames(obj)]

png(file.path(results_dir, "spatial_cot_vs_germ.png"),
    width = 2400, height = 800, res = 150)
print(SpatialFeaturePlot(obj, features = present, pt.size.factor = 10, ncol = length(present)))
dev.off()
