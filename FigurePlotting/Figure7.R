# This script will create figures 7A-G from the paper and save them to files.
# Author: Jaclyn Beck
# Final script used for paper as of Sep 09, 2022

library(Seurat)
library(monocle3)
library(SeuratWrappers) #remotes::install_github('satijalab/seurat-wrappers')
source(file.path("functions", "Figure_PlottingFunctions.R"))
source("Filenames.R")

##### Colors #####

# WT (green), 5X (red), PS-19 (blue), PS-5X (purple)
geno.colors <- c("#51B400", "#DE3923", "#4B5BDD", "#6B3D92")

# Should be colorbind/black and white printing safe:
# Light red, mustard, light blue, green, purple, magenta, cyan, orange, light purple
clust.colors <- c("#EF766D", "#E5B920", "#009BED", "#099600", "#6560BF", 
                  "#DF588F", "#00CBFB", "#F7962B", "#A5A0FF")


##### Setup #####

scRNA <- readRDS(file_seurat_analyzed_cd8)

DefaultAssay(scRNA) <- "RNA"

scRNA$genotype <- factor(scRNA$genotype, 
                         levels = c("WT", "5XFAD", "PS19", "PS-5X"))

excluded_genes <- readRDS(file_excluded_genes)

# We don't need the integrated assay
scRNA <- DietSeurat(scRNA, assays = c("RNA"), dimreducs = c("pca", "umap"))
gc()

Idents(scRNA) <- scRNA$clusters


##### Figure 7A #####

# UMAP with labels in the legend only
plt <- DimPlot(scRNA, reduction = "umap", shuffle = TRUE, label = FALSE) +
  scale_color_discrete(type = clust.colors)

plt

ggsave(file.path(dir_figures_paper, "fig7a.png"),
       plot = plt, width = 5.9, height = 3.0, units = "in", dpi = "print")


##### Figure 7B #####

# UMAP of clusters split by genotype, downsampled to have even cell numbers
plt <- umapSplitDownsampled(scRNA, ident.name = "clusters", clust.colors)
plt

ggsave(file.path(dir_figures_paper, "fig7b.png"),
       plot = plt, width = 4.9, height = 4.0, units = "in", dpi = "print")


##### Figure 7C #####

# Dotplot of top 5 upregulated genes by cluster
plt <- dotplotTop5Markers(scRNA, file_markers_cd8_all, excluded_genes, FDR = 0.01)
plt

ggsave(file.path(dir_figures_paper, "fig7c.png"),
       plot = plt, width = 5, height = 6, units = "in", dpi = "print")


##### Figure 7D #####

# Population bar graph. Significance values are from DPA analysis.
# A lot of values in this section are hard-coded and would need to be changed
# if the data changes. 

# centers of each bar relative to the center of the group (-0.33 to +0.33)
bars <- c("WT" = -0.23, "5XFAD" = -0.08, "PS19" = 0.08, "PS-5X" = 0.23)
sig.df <- 
  data.frame(cluster = rep("Naive", 2), 
             group1 = c("5XFAD", "PS19"), 
             group2 = c("PS19",  "PS-5X"),
             p.adj = c(0.01918,  0.02992),
             xmin = 2 + c(bars["5XFAD"], bars["PS19"]),
             xmax = 2 + c(bars["PS19"],  bars["PS-5X"]), 
             y.position = c(31.5, 31.5)) %>% 
  rstatix::add_significance(p.col = "p.adj")

plt <- populationBarGraph(scRNA, geno.colors, sig.df)

plt 

ggsave(file.path(dir_figures_paper, "fig7d.png"),
       plot = plt, width = 6, height = 3, units = "in", dpi = "print")


##### Figure 7E #####

# UMAPs with the top 10 expanded clonotypes, split by genotype
scRNA <- assignTop10Clonos( scRNA )
plt <- umapTop10ClonotypesSplit(scRNA, ncol = 2)

plt

ggsave(file.path(dir_figures_paper, "fig7e.png"),
       plot = plt, width = 7.8, height = 4.0, units = "in", dpi = "print")


##### Figure 7F #####

# Pseudotime UMAP
cds <- as.cell_data_set(scRNA)

cds <- cluster_cells(cds = cds, reduction_method = "UMAP", random_seed = 1234,
                     resolution = 1e-3)
cds <- learn_graph(cds)
cds <- order_cells(cds, reduction_method = "UMAP", root_pr_nodes = c("Y_199"))

plt <- plot_cells(cds,
                  color_cells_by = "pseudotime",
                  label_cell_groups=FALSE,
                  label_leaves=FALSE,
                  label_branch_points=FALSE,
                  label_principal_points = FALSE,
                  graph_label_size=1.5)

plt

ggsave(file.path(dir_figures_paper, "fig7f.png"),
       plot = plt, width = 5, height = 2.8, units = "in", dpi = "print")


##### Figure 7G #####

# Exhaustion genes in pseudotime

rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name

cds <- estimate_size_factors(cds)
cds <- order_cells(cds, reduction_method = "UMAP", root_pr_nodes = c("Y_192"))

# Get only the clusters located where cytotoxic 1&2 and exhausted cells are
ok <- clusters(cds)
ok <- ok[ok %in% c(1, 3, 7, 11, 4)]
ok <- ok[names(ok) %in% colnames(scRNA)[scRNA$clusters %in% 
                                          c("Cytotoxic 1", "Cytotoxic 2", "Exhausted")]]
cds_subset <- cds[, names(ok)]

exh <- c("Pdcd1", "Ctla4", "Tigit", "Ifng", "Irf4", "Lag3", "Ccl3", "Ccl4")

cds_subset <- cds_subset[exh,]
colData(cds_subset)$clusters <- factor(colData(cds_subset)$clusters, 
                                       levels = c("Cytotoxic 1", "Cytotoxic 2", "Exhausted"))

plt <- plot_genes_in_pseudotime(cds_subset,
                                color_cells_by = "clusters",
                                min_expr = 0.5,
                                ncol = 4,
                                cell_size = 0.1,
                                horizontal_jitter = 0.1,
                                vertical_jitter = 0.1) + 
  scale_color_discrete(type = clust.colors[c(1, 5, 4)])

plt

ggsave(file.path(dir_figures_paper, "fig7g.png"),
       plot = plt, width = 8, height = 4, units = "in", dpi = "print")


rm(list = ls())
gc()
