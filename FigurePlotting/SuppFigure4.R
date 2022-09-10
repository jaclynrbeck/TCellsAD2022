# This script will create supplemental figures 4A-B from the paper and save them
# to files.
# Author: Jaclyn Beck
# Final script used for paper as of Sep 09, 2022

library(Seurat)
source(file.path("functions", "Figure_PlottingFunctions.R"))
source("Filenames.R")

##### Colors #####

# Should be colorbind/black and white printing safe:
# Light red, mustard, light blue, green, purple, magenta, cyan, orange, light purple
cd8.clust.colors <- c("#EF766D", "#E5B920", "#009BED", "#099600", "#6560BF", 
                      "#DF588F", "#00CBFB", "#F7962B", "#A5A0FF")

# Should be colorbind/black and white printing safe:
# Light red, mustard, light blue, green, purple, cyan
cd4.clust.colors <- c("#FF867D", "#E5B920", "#009BED", "#099600", "#6560BF", 
                      "#00CBFB")


##### Supplemental Figure 4A #####

scRNA.cd8 <- readRDS(file_seurat_analyzed_cd8)

DefaultAssay(scRNA.cd8) <- "RNA"

# We don't need the integrated assay
scRNA.cd8 <- DietSeurat(scRNA.cd8, assays = c("RNA"), dimreducs = c("pca", "umap"))
gc()

genes <- c("Sell", "Ccr7", "Cd44", # Naive/CM vs effector
           "Pdcd1", "Tigit", # Exhaustion
           "Ifng", "Gzmk", "Ccl5", # Cytotoxic
           "Cx3cr1", "Klrg1", # Activation
           "Stat1", "Ifit3", # IFN response
           "Mki67" # Proliferating
) # 13 genes

# 13 genes but ncol=5 leaves 2 spaces at the end for cluster UMAPs
plt <- umapGeneExprPlots(scRNA.cd8, genes, coord.fixed = TRUE)
plt

ggsave(file.path(dir_figures_paper, "suppfig4a_1.png"),
       plot = plt, width = 7.5, height = 10, units = "in", dpi = "print")

# Project clusters onto all cells UMAP
Idents(scRNA.cd8) <- scRNA.cd8$clusters

plt <- umapProjectClustsToAllCells(scRNA.cd8, cd8.clust.colors)
plt

ggsave(file.path(dir_figures_paper, "suppfig4a_2.png"),
       plot = plt, width = 4.9, height = 2.7, units = "in", dpi = "print")

rm(scRNA.cd8, plt)
gc()

##### Supplemental Figure 4B #####

scRNA.cd4 <- readRDS(file_seurat_analyzed_cd4)

DefaultAssay(scRNA.cd4) <- "RNA"

# We don't need the integrated assay
scRNA.cd4 <- DietSeurat(scRNA.cd4, assays = c("RNA"), dimreducs = c("pca", "umap"))
gc()

genes <- c("Sell", "Cd44", # Naive/CM vs effector
           "Mki67", # Proliferating
           "Ifng", "Cxcr3", # Th1
           "Il4", "Il13", # Th2
           "Rorc", "Il17a", # Th17
           "Foxp3", # T-reg
           "Gzmk", # Cytotoxic
           "Cx3cr1", # Activation
           "Ifit3" # IFN response
) # 13 genes

# 13 genes but ncol=5 leaves 2 spaces at the end for cluster UMAPs
plt <- umapGeneExprPlots(scRNA.cd4, genes, coord.fixed = TRUE)
plt

ggsave(file.path(dir_figures_paper, "suppfig4b_1.png"),
       plot = plt, width = 7.5, height = 10, units = "in", dpi = "print")

# Project clusters onto all cells UMAP
Idents(scRNA.cd4) <- scRNA.cd4$clusters

plt <- umapProjectClustsToAllCells(scRNA.cd4, cd4.clust.colors)
plt

ggsave(file.path(dir_figures_paper, "suppfig4b_2.png"),
       plot = plt, width = 4.6, height = 2.7, units = "in", dpi = "print")


rm(list = ls())
gc()
