# Integrates CD8+ cells using the Seurat Integration pipeline, then 
# clusters cells and generates the following data:
#   Differential genes between clusters
#   Differential genes between genotypes
#   List of clonotypes per cluster
#   List of clonotypes per genotype
#   Significant GO terms based on cluster DGE (gprofiler, not used in paper)
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


##### Get CD8+ cells only #####

scRNA <- readRDS(file_seurat_unnorm)

single.pos <- getSinglePositiveCells(scRNA)
cd8.pos <- single.pos[["CD8"]]
scRNA <- scRNA[,cd8.pos]


##### Scale and integrate w/ regression #####

DefaultAssay(scRNA) <- "RNA"
scRNA <- generateIntegratedData(scRNA, method  = "SCT",
                                regress = c("nCount_RNA", "StressScoreVariable"))

# Save intermediate step
saveRDS(scRNA, file = file_seurat_norm_cd8)

# Get rid of SCT assay to save memory / space
scRNA <- DietSeurat(scRNA, assays = c("RNA", "integrated"), scale.data = TRUE)
gc()


##### Run PCA and cluster #####

scRNA <- RunPCA(scRNA, npcs = 100)

# Examining how many PCs we really need
ElbowPlot(scRNA, ndims = 100)

# For the paper: 50 dimensions needed to separate out some clusters
scRNA <- RunUMAP(scRNA, dims = 1:50)
scRNA <- FindNeighbors(scRNA, reduction="pca", dims = 1:50)
scRNA <- FindClusters(scRNA, resolution = 0.4)

DimPlot(scRNA, reduction = "umap", shuffle = TRUE, label = TRUE) + NoLegend()


##### Annotate clusters manually #####

clusters <- list("0" = "Cytotoxic 2",
                 "1" = "Naive",
                 "2" = "CM",
                 "3" = "Exhausted",
                 "4" = "Cytotoxic 1",
                 "5" = "Activated",
                 "6" = "Early Activation",
                 "7" = "IFN Responding",
                 "8" = "Proliferating")

scRNA <- RenameIdents(scRNA, clusters)
scRNA$clusters <- Idents(scRNA)

DimPlot(scRNA, reduction = "umap", shuffle = TRUE, label = TRUE) + NoLegend()

# Save integrated/annotated object
saveRDS(scRNA, file_seurat_analyzed_cd8)


##### Find cluster markers #####

# Print most expressed genes in each cluster
highest.expressed <- getHighestExpressedGenes(scRNA)
highest.expressed

# Get all diff genes in each cluster
DefaultAssay(scRNA) <- "RNA"
scRNA <- NormalizeData(scRNA)

all.markers <- FindAllMarkers(scRNA, test.use = "MAST", 
                              logfc.threshold = 0.25, min.pct = 0.1,
                              latent.vars = c("nCount_RNA", "StressScoreVariable"))

saveRDS(all.markers, file = file_markers_cd8_all)

# Read in clonotype data
tcr.anno <- read.csv(file_clonotypes_processed)
table(tcr.anno$Genotype)

tcr.match <- subset(tcr.anno, Seurat.Barcode %in% colnames(scRNA))
table(tcr.match$Genotype)

# Save significant diff genes in each cluster
FDR = 0.01
sig.markers = filter(all.markers, p_val_adj <= FDR & pct.1 >= 0.1)
table(sig.markers$cluster)

writeDifferentialGenes(sig.markers, file_markers_cd8_clusters)

top10 <- printTop10Markers(sig.markers, pos.only = TRUE)

# Which clonotypes are in each cluster
writeClusterClonotypes(scRNA, tcr.anno, file_clonotypes_summary,
                       file_clonotypes_cd8_clusters)

# Diff genes and clonotypes by genotype
genotypes <- unique(scRNA$genotype)
writeGenotypeDifferentialGenes(scRNA, genotypes, FDR, file_markers_cd8_genotypes)
writeGenotypeClonotypes(scRNA, tcr.anno, genotypes, 
                        file_clonotypes_summary,
                        file_clonotypes_cd8_genotypes)

writeClusterVGenotypeDiffGenes(scRNA, genotypes, FDR, file_markers_cd8_clusters_vs_genotype)


##### Look at diff between Cytotoxic 1, Cytotoxic 2, and Exhausted Clusters #####

cyto.markers <- FindMarkers(scRNA, ident.1 = "Cytotoxic 1", 
                            ident.2 = "Cytotoxic 2", test.use = "MAST",
                            logfc.threshold = 0.25, min.pct = 0.1,
                            latent.vars = c("nCount_RNA", "StressScoreVariable"))

cyto.markers <- filter(cyto.markers, p_val_adj <= FDR & pct.1 >= 0.1)
cyto.markers <- cyto.markers[order(-cyto.markers$avg_log2FC),]
cyto.markers$Ensembl.Id <- geneNameToEnsembl(rownames(cyto.markers))
write.csv(cyto.markers, file.path(dir_cd8, "DiffGenes_Cytotoxic1v2.csv"))

cyto.markers <- FindMarkers(scRNA, ident.1 = "Cytotoxic 2", 
                            ident.2 = "Exhausted", test.use = "MAST",
                            logfc.threshold = 0.25, min.pct = 0.1,
                            latent.vars = c("nCount_RNA", "StressScoreVariable"))

cyto.markers <- filter(cyto.markers, p_val_adj <= FDR & pct.1 >= 0.1)
cyto.markers <- cyto.markers[order(-cyto.markers$avg_log2FC),]
cyto.markers$Ensembl.Id <- geneNameToEnsembl(rownames(cyto.markers))
write.csv(cyto.markers, file.path(dir_cd8, "DiffGenes_Cytotoxic2vExhausted.csv"))


##### Everything beyond this point was not used in the paper #####

##### GO Analysis - By cluster to help annotate. Not used for paper. #####

all.markers <- readRDS(file_markers_cd8_all)
sig.genes.df <- filter(all.markers, p_val_adj <= 0.01)
sig.genes.df$Ensembl.ID <- geneNameToEnsembl(sig.genes.df$gene)

for (C in unique(sig.genes.df$cluster)) {
  markers <- subset(sig.genes.df, cluster == C)
  
  go.res <- runGOAnalysis(markers, sources = c("GO:BP"))
  saveRDS(go.res, 
          file = file.path(dir_cd8_go,
                           paste0("gprofiler_cd8_", C, ".rds")))
  
  goResultsToGem(go.res, file.path(dir_cd8_go, 
                                   paste0("gProfiler_gem_", C, ".txt")))
}

# GEM file is used with the EnrichmentMap plugin for Cytoscape
