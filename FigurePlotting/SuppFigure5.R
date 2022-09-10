# This script will create supplemental figures 5A-F from the paper and save them
# to files.
# Author: Jaclyn Beck
# Final script used for paper as of Sep 09, 2022

library(Seurat)
source(file.path("functions", "Figure_PlottingFunctions.R"))
source("Filenames.R")


##### Colors #####

# WT (green), 5X (red), PS-19 (blue), PS-5X (purple)
geno.colors <- c("#51B400", "#DE3923", "#4B5BDD", "#6B3D92") 

# Should be colorbind/black and white printing safe:
# Light red, mustard, light blue, green, purple, cyan
clust.colors <- c("#FF867D", "#E5B920", "#009BED", "#099600", "#6560BF", "#00CBFB")


##### Setup #####

scRNA <- readRDS(file_seurat_analyzed_cd4)

DefaultAssay(scRNA) <- "RNA"

scRNA$genotype <- factor(scRNA$genotype, 
                         levels = c("WT", "5XFAD", "PS19", "PS-5X"))

excluded_genes <- readRDS(file_excluded_genes)

# We don't need the integrated assay
scRNA <- DietSeurat(scRNA, assays = c("RNA"), dimreducs = c("pca", "umap"))
gc()

Idents(scRNA) <- scRNA$clusters


##### Supplemental Figure 5A #####

# UMAP with labels in the legend only
plt <- DimPlot(scRNA, reduction = "umap", shuffle = TRUE, label = FALSE) +
  scale_color_discrete(type = clust.colors)

plt

ggsave(file.path(dir_figures_paper, "suppfig5a.png"),
       plot = plt, width = 5.25, height = 3.8, units = "in", dpi = "print")


##### Supplemental Figure 5B #####

# UMAP of clusters split by genotype, downsampled to have even cell numbers
plt <- umapSplitDownsampled(scRNA, ident.name = "clusters", clust.colors)
plt

ggsave(file.path(dir_figures_paper, "suppfig5b.png"),
       plot = plt, width = 4.0, height = 4.3, units = "in", dpi = "print")


##### Supplemental Figure 5C #####

# Dotplot of top 5 upregulated genes by cluster
plt <- dotplotTop5Markers(scRNA, file_markers_cd4_all, excluded_genes, FDR = 0.01)
plt

ggsave(file.path(dir_figures_paper, "suppfig5c.png"),
       plot = plt, width = 5, height = 4.5, units = "in", dpi = "print")


##### Supplemental Figure 5D #####

# Population bar graph. Not enough cells to do DPA, so no significance markers.
plt <- populationBarGraph(scRNA, geno.colors, sig.df = NULL)

plt 

ggsave(file.path(dir_figures_paper, "suppfig5d.png"),
       plot = plt, width = 5, height = 3, units = "in", dpi = "print")


##### Supplemental Figure 5E #####

# UMAPs of cytotoxic genes
plt <- FeaturePlot(scRNA, features = c("Prf1", "Gzma", "Gzmb", "Gzmc", 
                                       "Gzmk", "Gzmm"),
                   ncol = 3, order = TRUE, cols = c("lightgray", "red"))
plt

ggsave(file.path(dir_figures_paper, "suppfig5e.png"),
       plot = plt, width = 9, height = 5, units = "in", dpi = "print")


##### Supplemental Figure 5F #####

# UMAPs with the top 10 expanded clonotypes, split by genotype, wide format
scRNA <- assignTop10Clonos( scRNA )
plt <- umapTop10ClonotypesSplit(scRNA, ncol = 4)

plt

ggsave(file.path(dir_figures_paper, "suppfig5f.png"),
       plot = plt, width = 10.1, height = 2.5, units = "in", dpi = "print")


rm(list = ls())
gc()
