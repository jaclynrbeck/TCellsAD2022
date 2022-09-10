# This script will create supplemental figure 3 from the paper and save it
# to a file.
# Author: Jaclyn Beck
# Final script used for paper as of Sep 09, 2022

library(Seurat)
source(file.path("functions", "Figure_PlottingFunctions.R"))
source("Filenames.R")

##### Setup #####

scRNA <- readRDS(file_seurat_analyzed_allcells)

DefaultAssay(scRNA) <- "RNA"

# We don't need the integrated assay
scRNA <- DietSeurat(scRNA, assays = c("RNA"), dimreducs = c("pca", "umap"))
gc()


##### Supplemental Figure 3 #####

genes <- c("Cd8b1", "Cd4", "Trdc",  # Cell type
           "Ifng", # Type 1
           "Il17a", # Type 17
           "Sell", "Ccr7", "Cd44", # Naive/CM vs effector
           "Il7r", # Memory 
           "Mki67", # Proliferating
           "Prf1", "Gzmk", "Ccl5", # Cytotoxic
           "Tyrobp", "Klrb1c", # Gamma delta-1
           "Gzma", "Cx3cr1", "Klrg1", # Activation
           "Cd69", "Itga1", # Residence
           "Stat1", "Bst2", "Ifit3", # Interferon response
           "Pdcd1", "Tigit" # Exhaustion
) # 25 genes

plt <- umapGeneExprPlots(scRNA, genes, coord.fixed = TRUE)
plt

ggsave(file.path(dir_figures_paper, "suppfig3.png"),
       plot = plt, width = 7.5, height = 10, units = "in", dpi = "print")


rm(list = ls())
gc()
