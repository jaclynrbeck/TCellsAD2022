# This script was written to be run line-by-line because of exploratory plots
# to determine parameters. The values as entered in the script are those used
# for the paper.

# Author: Jaclyn Beck
# Final script used for paper as of Feb 26, 2022

library(Seurat)
library(dplyr)
library(stringr)
library(writexl)
library(ggplot2)
source(file.path("functions", "Analysis_HelperFunctions.R"))
source(file.path("functions", "Analysis_PlottingFunctions.R"))
source("Filenames.R")
source("GeneMarkers.R")

scRNA <- readRDS(file_seurat_norm)
DefaultAssay(scRNA)

tcr.anno <- read.csv(file_clonotypes_processed)
table(tcr.anno$Genotype)

tcr.match <- matchTcrs(scRNA, tcr.anno)
table(tcr.match$Genotype)

scRNA <- RunPCA(scRNA, npcs = 100)

# Examining how many PCs we really need

plotPCACurves(scRNA)
ElbowPlot(scRNA, ndims = 100)

# For paper: 25 dimensions is sufficient
scRNA <- RunUMAP(scRNA, dims = 1:25)
scRNA <- FindNeighbors(scRNA, reduction="pca", dims = 1:25)
scRNA <- FindClusters(scRNA, resolution = 0.1)

VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 
                            "percent.rb"), ncol=2, pt.size = 0)

DimPlot(scRNA, reduction = "umap", shuffle = TRUE, label = TRUE) + NoLegend()

# Look at the 20 most highly variable genes
top20 <- head(VariableFeatures(scRNA), 20)
top20

# Print cluster distribution information
printClusterDistributions(scRNA)

clusters <- list("0" = "CD8 Cytotoxic",
                 "1" = "Naive",
                 "2" = "CD4 Helper",
                 "3" = "CD8 CM + \u03b3\u03b4",
                 "4" = "\u03b3\u03b4", 
                 "5" = "CD8 Activated",
                 "6" = "CD8 IFN Response",
                 "7" = "Proliferating")
scRNA <- RenameIdents(scRNA, clusters)
scRNA$clusters <- Idents(scRNA)

saveRDS(scRNA, file_seurat_analyzed_allcells)

##### Find cluster markers #####

DefaultAssay(scRNA) <- "RNA"
scRNA <- NormalizeData(scRNA)

all.markers <- FindAllMarkers(scRNA)
saveRDS(all.markers, file = file_markers_combined_all)

#all.markers <- readRDS(file_markers_combined_all)

# Most expressed genes in each cluster
highest.expressed <- getHighestExpressedGenes(scRNA)
highest.expressed

# Get all diff genes in each cluster
FDR = 0.05
sig.markers = filter(all.markers, p_val_adj <= FDR & pct.1 >= 0.1)
table(sig.markers$cluster)

writeDifferentialGenes(sig.markers, file_markers_combined_clusters)

genotypes <- unique(scRNA$genotype)
writeGenotypeDifferentialGenes(scRNA, genotypes, FDR, file_markers_combined_genotypes)
writeGenotypeClonotypes(scRNA, tcr.anno, genotypes, file_clonotypes_raw,
                        file_clonotypes_combined_genotypes)

writeClusterVGenotypeDiffGenes(scRNA, genotypes, FDR, file_markers_combined_clusters_vs_genotype)

# Which clonotypes are in each cluster
writeClusterClonotypes(scRNA, tcr.anno, file_clonotypes_raw,
                       file_clonotypes_combined_clusters)

##### Some visualization #####

top10 <- printTop10Markers(sig.markers, pos.only = TRUE)

all.genes <- rownames(scRNA)

FeaturePlot(scRNA, features = c("Cd8a", "Cd8b1", "Cd4", "Trdv4"), min.cutoff = 0)
FeaturePlot(scRNA, features = c("Cd8a", "Cd4"), blend = TRUE, 
            cols = c("lightgrey", "red", "blue"), blend.threshold = 0)
FeaturePlot(scRNA, features = c("nFeature_RNA"))

VlnPlot(scRNA, features = c("Cd8a", "Cd8b1", "Cd4", "Cd3e", "Cd3d", "Cd3g"), 
        ncol=2, pt.size = 0)

# Cluster identification -- This code is old an may not work anymore
markers.pos <- getCellTypeMarkers(rownames(scRNA))
markers.pos <- convertMarkerListToDf(markers.pos)

sig.filtered.pos = filter(sig.markers, (gene %in% unique(markers.pos$Gene) & avg_log2FC > 0))
sig.filtered.neg = filter(sig.markers, (gene %in% unique(markers.pos$Gene) & avg_log2FC < 0))

printCellTypeAssignments(sig.markers.pos, sig.markers.neg)

DimPlot(scRNA, reduction = "umap", shuffle = TRUE, label = TRUE, repel = TRUE) + NoLegend()
DimPlot(scRNA, reduction = "umap", shuffle = TRUE, label = FALSE, repel = TRUE) 
DimPlot(scRNA, reduction = "umap", group.by="orig.ident", shuffle = TRUE)
DimPlot(scRNA, reduction = "umap", group.by="orig.ident", split.by = "orig.ident") + NoLegend()

single.pos <- getSinglePositiveCells(scRNA)

types <- data.frame(Cell = single.pos[["CD8"]], Type = "CD8")
types <- rbind(types, data.frame(Cell = single.pos[["CD4"]], Type = "CD4"))
types <- rbind(types, data.frame(Cell = single.pos[["GD"]], Type = "GD"))
types <- rbind(types, data.frame(Cell = single.pos[["Double"]], Type = "Double Pos"))
rownames(types) <- types$Cell
id <- types[colnames(scRNA), 2]
id[is.na(id)] <- "Unknown"

Idents(scRNA) <- id


##### GO Analysis - Genotypes #####

sig.genes.df <- readSigGenesGenotype(file_markers_combined_genotypes, "vs WT")

for (C in unique(sig.genes.df$cluster)) {
  markers <- subset(sig.genes.df, cluster == C)
  
  go.res <- runGOAnalysis(markers, sources = c("GO:BP"))
  saveRDS(go.res, 
          file = file.path(dir_allcells_go, 
                           paste0("gprofiler_allcells_", C, ".rds")))

  goResultsToGem(go.res, file.path(dir_allcells_go, 
                                   paste0("gProfiler_gem_", C, ".txt")))
}

# GEM file is used with the EnrichmentMap plugin for Cytoscape


##### TODO: DPA #####

