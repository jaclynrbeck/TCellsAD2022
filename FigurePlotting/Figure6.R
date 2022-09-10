# This script will create figures 6A-B from the paper and save them to files.
# Author: Jaclyn Beck
# Final script used for paper as of Sep 09, 2022

library(Seurat)
source(file.path("functions", "Figure_PlottingFunctions.R"))
source("Filenames.R")


##### Setup #####

scRNA <- readRDS(file_seurat_analyzed_allcells)

DefaultAssay(scRNA) <- "RNA"

scRNA$genotype <- factor(scRNA$genotype, 
                         levels = c("WT", "5XFAD", "PS19", "PS-5X"))

# We don't need the integrated assay
scRNA <- DietSeurat(scRNA, assays = c("RNA"), dimreducs = c("pca", "umap"))
gc()

scRNA <- assignTop10Clonos( scRNA )


##### Figure 6A #####

# UMAP with the top 10 expanded clonotypes, by percentage of sample
plt <- umapTop10Clonotypes(scRNA)
plt

ggsave(file.path(dir_figures_paper, "fig6a.png"),
       plot = plt, width = 7.3, height = 3.8, units = "in", dpi = "print")


##### Figure 6B #####

# UMAPs with the top 10 expanded clonotypes, split by genotype

plt <- umapTop10ClonotypesSplit(scRNA, ncol = 2)

plt

ggsave(file.path(dir_figures_paper, "fig6b.png"),
       plot = plt, width = 7.2, height = 4.25, units = "in", dpi = "print")


rm(list = ls())
gc()

