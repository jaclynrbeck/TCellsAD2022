# This script will create supplemental figures 6A-E from the paper and save them
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
# Magenta, mustard, purple
clust.colors <- c("#DF588F", "#E5B920", "#6560BF")


##### Setup #####

scRNA <- readRDS(file_seurat_analyzed_gd)

DefaultAssay(scRNA) <- "RNA"

scRNA$genotype <- factor(scRNA$genotype, 
                         levels = c("WT", "5XFAD", "PS19", "PS-5X"))

excluded_genes <- readRDS(file_excluded_genes)

# We don't need the integrated assay
scRNA <- DietSeurat(scRNA, assays = c("RNA"), dimreducs = c("pca", "umap"))
gc()

Idents(scRNA) <- scRNA$clusters_fine


##### Supplemental Figure 6A #####

# UMAP with labels in the legend only. Labels are manually added to the figure
plt <- DimPlot(scRNA, reduction = "umap", shuffle = TRUE, label = FALSE,
               label.size = 3, pt.size = 0.1) +
  scale_color_discrete(type = clust.colors)

plt

ggsave(file.path(dir_figures_paper, "suppfig6a.png"),
       plot = plt, width = 5.1, height = 3.0, units = "in", dpi = "print")


##### Supplemental Figure 6B #####

# UMAP colored by genotype
Idents(scRNA) <- scRNA$genotype
plt <- DimPlot(scRNA, reduction = "umap", shuffle = TRUE, pt.size = 0.1) +
  scale_color_discrete(type = geno.colors)

plt

ggsave(file.path(dir_figures_paper, "suppfig6b.png"),
       plot = plt, width = 4.35, height = 3.0, units = "in", dpi = "print")


##### Supplemental Figure 6C #####

# Dotplot of top 5 upregulated genes by cluster
Idents(scRNA) <- scRNA$clusters_fine
plt <- dotplotTop5Markers(scRNA, file_markers_gd_all, excluded_genes, FDR = 0.01)
plt

ggsave(file.path(dir_figures_paper, "suppfig6c.png"),
       plot = plt, width = 4, height = 3.1, units = "in", dpi = "print")


##### Supplemental Figure 6D #####

# Population bar graph. Not enough cells to do DPA, so no significance markers.
plt <- populationBarGraph(scRNA, geno.colors, sig.df = NULL)

plt 

ggsave(file.path(dir_figures_paper, "suppfig6d.png"),
       plot = plt, width = 4, height = 3, units = "in", dpi = "print")


##### Supplemental Figure 6E #####

# Gene expression UMAPS + projection of clusters onto all cells UMAP. The
# two sub-figures are manually combined with nice labels. 

genes <- c("Sell", # Naive/CM
           "Klrb1c", "Tyrobp", "Nkg7", "Ccl5", # Gamma delta-1 / NK-like
           "Rorc", "Blk", "Il17a") # Gamma delta-17

# 8 genes but ncol=5 leaves 2 spaces at the end for cluster UMAPs
plt <- umapGeneExprPlots(scRNA, genes, coord.fixed = FALSE)

plt

ggsave(file.path(dir_figures_paper, "suppfig6e_1.png"),
       plot = plt, width = 7.5, height = 3, units = "in", dpi = "print")

# Project clusters onto all cells UMAP

Idents(scRNA) <- scRNA$clusters_fine
plt <- umapProjectClustsToAllCells(scRNA, clust.colors)
plt

ggsave(file.path(dir_figures_paper, "suppfig6e_2.png"),
       plot = plt, width = 5, height = 2.7, units = "in", dpi = "print")


rm(list = ls())
gc()

