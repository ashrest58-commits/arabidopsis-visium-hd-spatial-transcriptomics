# Beginner Spatial Transcriptomics Workflow in R: Arabidopsis Visium HD

This repository contains a beginner-friendly and detailed R workflow for analyzing a 10x Genomics Visium HD 3' spatial transcriptomics dataset using Seurat. The workflow is built around a single *Arabidopsis thaliana* seedling section profiled at 8 µm bin resolution. It covers quality control, filtering, normalization, marker validation, unsupervised clustering, cluster annotation, spatially variable gene (SVG) detection, targeted differential expression, GO pathway analysis (ORA and GSEA), and spatial expression gradient modeling.

This script does not perform cell-type deconvolution or cell–cell communication analysis, because the dataset lacks cell segmentation (see *Important Notes* below).

---

## Project Goal

The goal of this project is to provide a clear, modular, reproducible spatial transcriptomics workflow for beginners who want to analyze Visium HD data without cell segmentation, using only a standard Seurat-based pipeline.

The workflow helps answer questions such as:

- Is the Visium HD section usable based on bin-level QC metrics?
- Can unsupervised clustering recover known seedling tissue domains without cell segmentation?
- Which marker genes support the biological identity of each cluster?
- Which genes show spatially structured expressions independent of cluster identity?
- How do tissue domains differ transcriptionally, and what pathways distinguish them?
- Does gene expression change continuously or sharply with distance from a developmental organizer (the shoot apical meristem)?

---

## Original Study and Data Source

This workflow uses the publicly available Visium HD 3' Arabidopsis seedling dataset from:

**10x Genomics. (2025). Visium HD 3' Spatial Gene Expression — *Arabidopsis thaliana* seedling [dataset].**
<https://www.10xgenomics.com/datasets/visium-hd-three-prime-arabidopsis>

Raw and binned outputs were generated with **Space Ranger 4.0.1** against the **TAIR10** reference genome and downloaded directly from 10x Genomics' public CDN (see `R_code/01_setup_and_load.R` for the exact download commands).

This repository is a beginner-oriented re-analysis and teaching workflow based on that publicly available dataset. The original data belong to 10x Genomics; please cite the dataset above if you use this workflow.

---

## Dataset Design

Unlike multi-sample case-control workflows, this repository analyzes a **single Visium HD section**, with no integration step and no condition comparison.

| Metric | Value |
|---|---|
| Bin size used | 8 µm |
| Total sequencing reads | 141,748,868 |
| Bins under tissue (8 µm) | 7,356 |
| Bins retained after QC filtering | 1,516 |
| Genes retained after filtering | 17,852 |
| Tissue domains recovered | 7 (Cotyledon, Germinating seed, Root elongation zone, Root, Meristematic tissue, Hypocotyl, SAM) |

---

## Main R Packages

```
Seurat
ggplot2
dplyr
hdf5r
arrow
ape
clusterProfiler
org.At.tair.db
mgcv
enrichplot
ggrepel
```

Most of these are not bundled with a default Seurat installation — see `R_code/02_environment_setup.sh` for the conda install commands used to add the missing Bioconductor package (`clusterProfiler`) and system dependencies (`hdf5r`, `arrow`, `ape`).

---

## Expected Input Folder Structure

Before running the scripts, update the `data_dir` path in `R_code/01_setup_and_load.R` to point to your downloaded and extracted Visium HD output.

```
spatial/arabidopsis_visium/
├── data/
│   ├── binned_outputs/
│   │   └── square_008um/
│   │       ├── filtered_feature_bc_matrix.h5
│   │       └── spatial/
│   └── metrics_summary.csv
└── results/
```

Space Ranger's `binned_outputs` folder contains pre-binned matrices at 2 µm, 8 µm, and 16 µm resolution in the same directory; `Load10X_Spatial(bin.size = 8)` selects the 8 µm resolution used throughout this workflow.

---

## How to Run the Workflow

1. Download the public dataset using the `wget` commands in `R_code/01_setup_and_load.R`.
2. Set up the `seurat_env` conda environment and install the additional packages listed in `R_code/02_environment_setup.sh`.
3. Update `data_dir` and `results_dir` paths at the top of each script.
4. Run the scripts in numeric order:

```
01_setup_and_load.R
02_QC_and_filtering.R
03_normalization_and_marker_validation.R
04_clustering_and_annotation.R
05_SVG_detection.R
06_targeted_DE.R
07_pathway_analysis_ORA_GSEA.R
08_spatial_gradient_analysis.R
```

Each script saves an `.rds` checkpoint that the next script reads back in, so the workflow can be stopped and resumed between steps — useful on a shared HPC cluster with walltime limits.

All output files are saved to:

```
Spatial_transcriptomics_output/
```

---

## Workflow Modules

### 1. Setup and Data Loading
Downloads the public dataset, sets up the working directory, and loads the 8 µm-binned matrix into a Seurat object with `Load10X_Spatial()`.

### 2. Quality Control and Filtering
Examines per-bin UMI and feature counts, removes off-tissue and low-complexity bins (nCount > 100, nFeature > 50), and removes genes with zero expression across all retained bins. Mitochondrial content is deliberately not used as a filtering criterion (see script comments for rationale).

### 3. Normalization and Marker Validation
Log-normalizes the filtered object, then checks known TAIR-documented tissue markers (e.g., `PHB`, `ERL2`) as an independent sanity check before any unsupervised analysis.

### 4. Clustering and Annotation
Performs variable feature selection, PCA, graph-based clustering, and UMAP, then assigns biological tissue-domain names to each numeric cluster using `FindAllMarkers()` output cross-referenced against the TAIR expression atlas.

### 5. Spatially Variable Gene (SVG) Detection
Uses Moran's I (`FindSpatiallyVariableFeatures(method = "moransi")`) to identify genes with non-random spatial expression independent of cluster identity, and checks their overlap with cluster marker genes.

### 6. Targeted Differential Expression
Runs 8 pairwise Wilcoxon rank-sum comparisons between specific tissue-domain pairs, with volcano plots and spatial feature plots of the leading genes per comparison.

### 7. Pathway Analysis (ORA and GSEA)
Runs Over-Representation Analysis (`compareCluster()` + `enrichGO()`) across all 7 clusters and pre-ranked GSEA (`gseGO()`) for 3 of the 8 pairwise DE comparisons, both restricted to a detected-gene background universe rather than the whole genome.

### 8. Spatial Gradient Modeling
Separates the slide's multiple physically distinct tissue pieces using k-means on spatial coordinates, isolates the piece spanning the full shoot developmental axis, and fits GAM models of gene expression against distance from the SAM centroid to distinguish sharply bounded vs. continuously graded expression.

---

## Main Output Files

### Seurat Objects (RDS checkpoints)

| File | Description |
|---|---|
| `arabidopsis_visium_filtered.rds` | Object after on-tissue and QC-threshold filtering |
| `arabidopsis_visium_normalized.rds` | Log-normalized object |
| `arabidopsis_visium_clustered.rds` | Object after PCA, UMAP, clustering (numeric cluster IDs) |
| `arabidopsis_visium_annotated.rds` | Final object with tissue-domain names assigned |

### Tables

| File | Description |
|---|---|
| `cluster_markers.csv` | Full positive marker table from `FindAllMarkers` |
| `top50_markers_per_cluster.csv` | Top 50 markers per cluster, used for annotation |
| `spatially_variable_genes.csv` | Moran's I statistics, all tested genes |
| `SVG_and_cluster_markers.csv` | Overlap between top SVGs and cluster markers |
| `DE_<clusterA>_vs_<clusterB>.csv` (×8) | Pairwise DE results |
| `ORA_compareCluster_all.csv` | Simplified GO BP ORA results, all clusters |
| `GSEA_<comparison>.csv` (×3) | Simplified GSEA results per comparison |

### Figures

| File | Description |
|---|---|
| `QC_UMI_plots.png` / `QC_after_filtering.png` | Pre- and post-filtering QC plots |
| `UMAP_clusters.png` / `spatial_clusters.png` | Numeric cluster UMAP and spatial overlay |
| `UMAP_annotated1_nolabel.png` / `spatial1_annotated.png` | Final annotated UMAP and spatial overlay |
| `cluster_heatmap1.png` | Heatmap of top 5 markers per cluster |
| `top6_SVGs.png`, `top11_16_SVGs.png`, `top17_20_SVGs.png` | Spatial plots of top 20 SVGs |
| `volcano_<comparison>.png` (×4) | Volcano plots |
| `spatial_<comparison>.png` | Spatial plots of top DE genes |
| `ORA_compareCluster_dotplot1.png` / `GSEA_dotplot_<comparison>.png` | Pathway dotplots |
| `spatial_key_markers_all_tissues.png` | One marker per tissue domain, plotted together |
| `gradient_<gene>.png` (×5) | GAM fit of expression vs. distance from SAM |

---

## Important Notes

- Update `data_dir` and `results_dir` in each script before running.
- This workflow assumes Space Ranger 4.0.1 binned outputs are already generated and extracted.
- Nuclei-based cell segmentation from the H&E image performed poorly on this plant section, due to plant cell morphology differing from the human/mouse tissue the segmentation model was trained on, and background staining from plant cell walls. As a result, this dataset is listed as a challenging tissue type on 10x Genomics' Visium HD 3' Tested Tissue List, and this workflow proceeds entirely at the bin level rather than the cell level.
- Cluster numbers are not stable identifiers; they may change if resolution, PCA dimensions, or neighbor parameters are modified. Biological identity should always be re-derived from marker genes and spatial position, not assumed from a cluster number alone.
- This Visium HD section contains multiple physically separate tissue pieces on one capture area. Any spatial-distance-based analysis (see Module 8) must first separate pieces by physical coordinates before computing distances, or risks mixing unrelated tissue fragments into one model.
- Mitochondrial-content filtering was intentionally skipped (see Module 2) because it can deplete biologically meaningful, metabolically active regions rather than removing only technical artifacts.
- Large raw data files and `.rds` Seurat objects are not included in this repository. Re-download the public dataset using the commands in `R_code/01_setup_and_load.R`, or use Git LFS / institutional storage for large derived files.

---

## GitHub Description

Beginner-friendly Seurat workflow for 10x Visium HD spatial transcriptomics analysis of an Arabidopsis seedling, including QC, filtering, marker validation, clustering, SVG detection, targeted DE, GO pathway analysis (ORA/GSEA), and spatial gradient modeling: without cell segmentation.

---

## Citation

10x Genomics. (2025). Visium HD 3' Spatial Gene Expression — *Arabidopsis thaliana* seedling [dataset]. https://www.10xgenomics.com/datasets/visium-hd-three-prime-arabidopsis

10x Genomics. Visium HD 3' Tested Tissue List.
Satija Lab. (2025). Visium HD analysis vignette. https://satijalab.org/seurat/articles/visiumhd_commands_intro
Satija Lab. Spatial transcriptomics vignette (archive v3.2). https://satijalab.org/seurat/archive/v3.2/spatial_vignette.html
The Arabidopsis Information Resource (TAIR) — Klepikova expression atlas.
Kim, H.-Y. (2014). Statistical notes for clinical researchers: Nonparametric statistical methods 1 — Nonparametric methods for comparing two groups.
Kim, H.-Y. (2014). Statistical notes for clinical researchers: Nonparametric statistical methods 2 — Nonparametric methods for the Wilcoxon rank-sum test.
Yates, J., et al. (2025). Filtering cells with high mitochondrial content depletes viable metabolically altered malignant cell populations in cancer single-cell studies. PMC11983838.
Topping, J. F., et al. (1997). Promoter trap markers differentiate structural and positional components of polar development in Arabidopsis. PMC157016.
Subramanian, A., et al. (2005). Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles.
Reimand, J., et al. (2019). Pathway enrichment analysis and visualization of omics data using g:Profiler, GSEA, Cytoscape and EnrichmentMap.
Wijesooriya, K., et al. (2022). Urgent need for consistent standards in functional enrichment analysis.
Wu, T., et al. (2021). clusterProfiler 4.0: A universal enrichment tool for interpreting omics data.
Jafari, M., & Ansari-Pour, N. (2019). Why, when and how to adjust your p values?
Chicco, D., & Agapito, G. (2022). Nine quick tips for pathway enrichment analysis.
Yu, G. (2024). clusterProfiler documentation / methods notes.
Zhou, et al. (2024). Spatial transcriptomics reveals unique metabolic profile and key oncogenic regulators of cervical squamous cell carcinoma.
de Oliveira, et al. (2026). Benchmarking multiple gene ontology enrichment tools reveals high biological significance, ranking, and stringency heterogeneity among datasets.
Xu, et al. (2024). Using clusterProfiler to characterize multiomics data.
Sagendorf, T. (2022). Over-Representation Analysis.
Matsuoka, et al. (2025). Single-cell and spatial transcriptomic characterization of pulmonary pleomorphic carcinoma.
Tokura, et al. (2026). Spatial Transcriptomics and Bulk RNA-Seq Analysis Revealed Molecular Classification of Invasive Lobular Carcinoma.
Zuo, et al. (2025). Spatial transcriptomic analysis of tumor microenvironment in esophageal squamous cell carcinoma with HIV infection.
Chen, et al. (2023). [clusterProfiler/ORA/GSEA application and methods reference — complete citation pending.]
ESRI. How Spatial Autocorrelation (Global Moran's I) works. https://doc.esri.com/en/arcgis-pro/latest/tool-reference/spatial-statistics/h-how-spatial-autocorrelation-moran-s-i-spatial-st.html
StatisticsHowTo. Moran's I: Definition, Examples. https://www.statisticshowto.com/morans-i/
