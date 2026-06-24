## ============================================================
## 08_pathway_analysis_ORA_GSEA.R
## GO Biological Process enrichment via ORA (compareCluster/
## enrichGO) across all clusters, and pre-ranked GSEA (gseGO)
## for 3 of the 8 pairwise DE comparisons.
## ============================================================

library(Seurat)
library(clusterProfiler)
library(org.At.tair.db)
library(dplyr)
library(ggplot2)
library(enrichplot)

results_dir <- "/path/to/spatial/arabidopsis_visium/results"
obj <- readRDS(file.path(results_dir, "arabidopsis_visium_annotated.rds"))

## ============================================================
## PART 1: ORA across all 7 clusters
## ============================================================

## Background universe = genes actually detected in this dataset
## (17,852 genes), NOT the whole TAIR genome (~27,000 genes) --
## most undetected genes have no real chance of being a marker
## and would only add noise to the hypergeometric test.
background_universe <- rownames(obj)

markers_annotated <- read.csv(file.path(results_dir, "cluster_markers.csv"), row.names = 1)

## Cluster-specific filtering: a uniform cutoff left Meristematic
## tissue with >1000 candidate genes and Hypocotyl with only 4 --
## too imbalanced for fair cross-cluster comparison.
filtered_standard <- markers_annotated %>%
  filter(p_val_adj < 0.05, avg_log2FC > 1.0, pct.1 > 0.25,
         !cluster %in% c("Meristematic tissue", "hypocotyl"))

filtered_meri <- markers_annotated %>%
  filter(p_val_adj < 0.05, avg_log2FC > 1.5, pct.1 > 0.40,
         cluster == "Meristematic tissue")

filtered_hypo <- markers_annotated %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.5, pct.1 > 0.10,
         cluster == "hypocotyl")

filtered_final <- bind_rows(filtered_standard, filtered_meri, filtered_hypo)
gene_clusters  <- split(filtered_final$gene, filtered_final$cluster)
print(sapply(gene_clusters, length))

cp_ora <- compareCluster(
  geneClusters  = gene_clusters,
  fun           = "enrichGO",
  OrgDb         = org.At.tair.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  universe      = background_universe,
  minGSSize     = 10,
  maxGSSize     = 500)

## Removes semantically redundant GO terms (parent terms inflated
## by overrepresented child terms), keeping the most significant
## representative term per group of similar terms.
cp_ora_simplified <- simplify(cp_ora, cutoff = 0.7, by = "p.adjust", select_fun = min)

write.csv(as.data.frame(cp_ora_simplified),
          file.path(results_dir, "ORA_compareCluster_all.csv"), row.names = FALSE)

png(file.path(results_dir, "ORA_compareCluster_dotplot1.png"), width = 2400, height = 1600, res = 150)
print(
  dotplot(cp_ora_simplified, showCategory = 5, title = "GO Biological Process - All tissue domains") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
)
dev.off()

## ============================================================
## PART 2: Pre-ranked GSEA for 3 pairwise comparisons
## ============================================================
## A pre-ranked approach (genes ranked by avg_log2FC from the
## pairwise DE comparison) is used because biological replicates
## are not available within a single spatial transcriptomics
## section -- this avoids requiring per-group variance estimates.

run_gsea <- function(de_csv, label, top_pathway_id = NULL) {
  de <- read.csv(file.path(results_dir, de_csv), row.names = 1)
  ranked <- sort(setNames(de$avg_log2FC, rownames(de)), decreasing = TRUE)

  set.seed(123)
  gsea_res <- gseGO(
    geneList      = ranked,
    OrgDb         = org.At.tair.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    minGSSize     = 10,
    maxGSSize     = 500,
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    seed          = TRUE,
    verbose       = TRUE)

  gsea_simplified <- simplify(gsea_res, cutoff = 0.7, by = "p.adjust", select_fun = min)

  saveRDS(gsea_res,        file.path(results_dir, paste0("gsea_", label, ".rds")))
  saveRDS(gsea_simplified, file.path(results_dir, paste0("gsea_", label, "_simplified.rds")))
  write.csv(as.data.frame(gsea_simplified),
            file.path(results_dir, paste0("GSEA_", label, ".csv")), row.names = FALSE)

  png(file.path(results_dir, paste0("GSEA_dotplot_", label, ".png")),
      width = 1400, height = 900, res = 150)
  print(
    dotplot(gsea_simplified, showCategory = nrow(as.data.frame(gsea_simplified)), split = ".sign") +
      facet_grid(. ~ .sign) +
      labs(title = paste("GSEA:", label), x = "Gene Ratio", y = NULL) +
      theme_bw(base_size = 12)
  )
  dev.off()

  if (!is.null(top_pathway_id)) {
    png(file.path(results_dir, paste0("GSEA_enrichplot_", label, ".png")),
        width = 1200, height = 800, res = 150)
    print(gseaplot2(gsea_res, geneSetID = top_pathway_id, title = "Top pathway enrichment"))
    dev.off()
  }

  gsea_res
}

## Cotyledon vs Germinating seed -- 22 pathways before simplify, 8 after.
## Top pathway: GO:0015979 (photosynthetic light reaction).
gsea_cot_germ <- run_gsea("DE_cotyledon_vs_germinating.csv", "cot_vs_germ", top_pathway_id = "GO:0015979")

## Cotyledon vs Hypocotyl -- 9 pathways before simplify, 4 after.
## Cotyledon: photosynthesis/metabolic process; Hypocotyl: transport.
gsea_cot_hyp  <- run_gsea("DE_cotyledon_vs_hypocotyl.csv", "cot_vs_hyp")

## SAM vs Cotyledon -- 41 pathways before simplify, 17 after.
## SAM: meristem development, shoot system development, reproductive
## process, flower development, cell cycle.
## Cotyledon: photosynthesis, plastid organization, ROS metabolic
## process, response to light stimulus.
gsea_sam_cot  <- run_gsea("DE_SAM_vs_cot.csv", "SAM_vs_cot")
