# Scales CD4+ cell data using SCTransform, then clusters cells and generates the
# following data:
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


##### Get CD4+ cells only #####

# For some reason I used the analyzed Seurat object instead of seurat_unnorm
# when originally analyzing this data, and the cells are in a different order 
# between the two. Changing the order alters UMAP/clustering, so I'm keeping
# this as using the analyzed object to be consistent with the paper.
scRNA <- readRDS(file_seurat_analyzed_allcells)

single.pos <- getSinglePositiveCells(scRNA)
cd4.pos <- single.pos[["CD4"]]
scRNA <- scRNA[,cd4.pos]


##### Scale w/ SCT and regression #####

DefaultAssay(scRNA) <- "RNA"
scRNA <- DietSeurat(scRNA, assays = c("RNA"))
scRNA <- removeLowExpressedGenes(scRNA)
scRNA <- SCTransform(scRNA, variable.features.n = 4000, 
                     vars.to.regress = c("nCount_RNA", "StressScoreVariable"))

# Save intermediate step
saveRDS(scRNA, file = file_seurat_norm_cd4)


##### Run PCA and cluster #####

scRNA <- RunPCA(scRNA, npcs = 100)

# Examining how many PCs we really need
ElbowPlot(scRNA, ndims = 100)

# For the paper: need 25 dimensions for this set
scRNA <- RunUMAP(scRNA, dims = 1:25)
scRNA <- FindNeighbors(scRNA, reduction="pca", dims = 1:25)
scRNA <- FindClusters(scRNA, resolution = 0.3)

DimPlot(scRNA, reduction = "umap", shuffle = TRUE, label = TRUE) + NoLegend()

# Cluster 5 is almost entirely WT-2, re-assign it to cluster 0 (Thelper)
scRNA$seurat_clusters[scRNA$seurat_clusters == 5] <- 0
scRNA$seurat_clusters <- factor(scRNA$seurat_clusters)
Idents(scRNA) <- scRNA$seurat_clusters


##### Annotate clusters manually #####

clusters <- list("0" = "T Helper",
                 "1" = "Naive",
                 "2" = "Cytotoxic",
                 "3" = "Activated", 
                 "4" = "T-regs",
                 "6" = "Proliferating")

scRNA <- RenameIdents(scRNA, clusters)
scRNA$clusters <- Idents(scRNA)

DimPlot(scRNA, reduction = "umap", shuffle = TRUE, label = TRUE) + NoLegend()

# Save integrated/annotated object
saveRDS(scRNA, file_seurat_analyzed_cd4)


##### Find cluster markers #####

# Get all diff genes in each cluster
DefaultAssay(scRNA) <- "RNA"
scRNA <- NormalizeData(scRNA)

all.markers <- FindAllMarkers(scRNA, test.use = "MAST", 
                              logfc.threshold = 0.25, min.pct = 0.1,
                              latent.vars = c("nCount_RNA", "StressScoreVariable"))

saveRDS(all.markers, file = file_markers_cd4_all)

# Read in clonotype data
tcr.anno <- read.csv(file_clonotypes_processed)
table(tcr.anno$Genotype)

tcr.match <- subset(tcr.anno, Seurat.Barcode %in% colnames(scRNA))
table(tcr.match$Genotype)

# Save significant diff genes in each cluster
FDR = 0.01
sig.markers = filter(all.markers, p_val_adj <= FDR & pct.1 >= 0.1)
table(sig.markers$cluster)

writeDifferentialGenes(sig.markers, file_markers_cd4_clusters)

top10 <- printTop10Markers(sig.markers, pos.only = TRUE)

# Which clonotypes are in each cluster
writeClusterClonotypes(scRNA, tcr.anno, file_clonotypes_summary,
                       file_clonotypes_cd4_clusters)

# Diff genes and clonotypes by genotype
genotypes <- unique(scRNA$genotype)
writeGenotypeDifferentialGenes(scRNA, genotypes, FDR, file_markers_cd4_genotypes)
writeGenotypeClonotypes(scRNA, tcr.anno, genotypes, 
                        file_clonotypes_summary,
                        file_clonotypes_cd4_genotypes)

writeClusterVGenotypeDiffGenes(scRNA, genotypes, FDR, file_markers_cd4_clusters_vs_genotype)

# Done
