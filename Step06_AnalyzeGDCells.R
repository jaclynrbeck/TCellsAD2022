# Scales Gamma Delta cell data using SCTransform, then clusters cells and 
# generates the following data:
#   Differential genes between clusters
#   Differential genes between genotypes
#   List of clonotypes per cluster
#   List of clonotypes per genotype
# This script was written to be run line-by-line because of exploratory plots
# to determine parameters. The values as entered in the script are those used
# for the paper.
# This script is very similar to Step03 but has enough modifications to
# require a new script.

# Author: Jaclyn Beck
# Final script used for paper as of Sep 08, 2022

library(Seurat)
library(dplyr)
library(stringr)
library(writexl)
library(ggplot2)
library(gprofiler2)
source(file.path("functions", "Analysis_HelperFunctions.R"))
source(file.path("functions", "General_HelperFunctions.R"))
source("Filenames.R")


##### Get Gamma-Delta cells only #####

scRNA <- readRDS(file_seurat_unnorm)

single.pos <- getSinglePositiveCells(scRNA)
gd.pos <- single.pos[["GD"]]
scRNA <- scRNA[,gd.pos]


##### Scale w/ SCT and regression #####

# For gamma delta cells we need to include orig.ident in the regression
# to avoid cells clustering by sample
DefaultAssay(scRNA) <- "RNA"
#scRNA <- DietSeurat(scRNA, assays = c("RNA"))
scRNA <- removeLowExpressedGenes(scRNA, threshold = 1)
scRNA <- SCTransform(scRNA, variable.features.n = 4000, 
                     vars.to.regress = c("nCount_RNA", "StressScoreVariable",
                                         "orig.ident"))

# Save intermediate step
saveRDS(scRNA, file = file_seurat_norm_gd)


##### Run PCA and cluster #####

scRNA <- RunPCA(scRNA, npcs = 100)

# Examining how many PCs we really need
ElbowPlot(scRNA, ndims = 100)

# For the paper: need 30 dimensions for this set
scRNA <- RunUMAP(scRNA, dims = 1:30)
scRNA <- FindNeighbors(scRNA, reduction="pca", dims = 1:30)
scRNA <- FindClusters(scRNA, resolution = 0.2)

# Separate Naive/CM-like cells from gamma delta-1 cells
scRNA <- FindSubCluster(scRNA, 1, graph.name = "SCT_snn", resolution = 0.2)

DimPlot(scRNA, reduction = "umap", label=TRUE) + NoLegend()

Idents(scRNA) <- scRNA$sub.cluster
DimPlot(scRNA, reduction = "umap", label=TRUE) + NoLegend()


##### Annotate clusters manually #####

clusters <- list("0" = "\u03b3\u03b4-17",
                 "1" = "\u03b3\u03b4-1")

clusters.sub <- list("0" = "\u03b3\u03b4-17",
                     "1_0" = "Naive/CM-like \u03b3\u03b4",
                     "1_1" = "\u03b3\u03b4-1")

Idents(scRNA) <- scRNA$seurat_clusters
scRNA <- RenameIdents(scRNA, clusters)
scRNA$clusters <- Idents(scRNA)
DimPlot(scRNA, reduction = "umap", label=TRUE) + NoLegend()

Idents(scRNA) <- scRNA$sub.cluster
scRNA <- RenameIdents(scRNA, clusters.sub)
scRNA$clusters_fine <- Idents(scRNA)
DimPlot(scRNA, reduction = "umap", label=TRUE) + NoLegend()

# Save integrated/annotated object
saveRDS(scRNA, file_seurat_analyzed_gd)


##### Find cluster markers #####

# Get all diff genes in each cluster (including subclusters)
Idents(scRNA) <- scRNA$clusters_fine
DefaultAssay(scRNA) <- "RNA"
scRNA <- NormalizeData(scRNA)

all.markers <- FindAllMarkers(scRNA, test.use = "MAST", 
                              logfc.threshold = 0.25, min.pct = 0.1,
                              latent.vars = c("nCount_RNA", "StressScoreVariable"))

saveRDS(all.markers, file = file_markers_gd_all)

# Read in clonotype data
tcr.anno <- read.csv(file_clonotypes_processed)
table(tcr.anno$Genotype)

tcr.match <- subset(tcr.anno, Seurat.Barcode %in% colnames(scRNA))
table(tcr.match$Genotype)

# Save significant diff genes in each cluster
FDR = 0.01
sig.markers = filter(all.markers, p_val_adj <= FDR & pct.1 >= 0.1)
table(sig.markers$cluster)

writeDifferentialGenes(sig.markers, file_markers_gd_clusters)

top10 <- printTop10Markers(sig.markers, pos.only = TRUE)

# Which clonotypes are in each cluster
writeClusterClonotypes(scRNA, tcr.anno, file_clonotypes_summary,
                       file_clonotypes_gd_clusters)

# Diff genes and clonotypes by genotype
genotypes <- unique(scRNA$genotype)
writeGenotypeDifferentialGenes(scRNA, genotypes, FDR, file_markers_gd_genotypes)
writeGenotypeClonotypes(scRNA, tcr.anno, genotypes, 
                        file_clonotypes_summary,
                        file_clonotypes_gd_genotypes)

# Fixes invalid character in sheet names for this method
Idents(scRNA) <- str_replace(Idents(scRNA), "/", "+")
writeClusterVGenotypeDiffGenes(scRNA, genotypes, FDR, file_markers_gd_clusters_vs_genotype)

# Done
