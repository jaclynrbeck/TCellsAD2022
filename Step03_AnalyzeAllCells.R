# Integrates all samples using the Seurat Integration pipeline, then 
# clusters cells and generates the following data:
#   Differential genes between clusters
#   Differential genes between genotypes
#   List of clonotypes per cluster
#   List of clonotypes per genotype
#   List of AD risk genes significantly changed between clusters
#   DPA to look for sig population differences (takes a long time)
#   Significant GO terms based on genotype DGE (gprofiler, not used in paper)
# This script was written to be run line-by-line because of exploratory plots
# to determine parameters. The values as entered in the script are those used
# for the paper.

# Author: Jaclyn Beck
# Final script used for paper as of Sep 09, 2022

library(Seurat)
library(dplyr)
library(stringr)
library(writexl)
library(ggplot2)
library(gprofiler2)
source(file.path("functions", "Analysis_HelperFunctions.R"))
source(file.path("functions", "General_HelperFunctions.R"))
source(file.path("functions", "DPA_HelperFunctions.R"))
source("Filenames.R")


##### Scale and integrate w/ regression #####

scRNA <- readRDS(file_seurat_unnorm)
scRNA <- generateIntegratedData(scRNA, method = "SCT",
                                regress = c("nCount_RNA", "StressScoreVariable"))

saveRDS(scRNA, file = file_seurat_norm) # Save intermediate step

DefaultAssay(scRNA)

# Get rid of SCT assay to save memory / space
scRNA <- DietSeurat(scRNA, assays = c("RNA", "integrated"), scale.data = TRUE)
gc()


##### Run PCA and cluster #####

scRNA <- RunPCA(scRNA, npcs = 100)

# Examining how many PCs we really need
ElbowPlot(scRNA, ndims = 100)

# For paper: 40 dimensions is sufficient
scRNA <- RunUMAP(scRNA, dims = 1:40)
scRNA <- FindNeighbors(scRNA, reduction="pca", dims = 1:40)
scRNA <- FindClusters(scRNA, resolution = 0.1)

VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 
                            "percent.rb"), ncol=2, pt.size = 0)

DimPlot(scRNA, reduction = "umap", shuffle = TRUE, label = TRUE) + NoLegend()

# Look at the 20 most highly variable genes
top20 <- head(VariableFeatures(scRNA), 20)
top20

# Print cluster distribution information
printClusterDistributions(scRNA)


##### Annotate clusters manually #####

clusters <- list("0" = "CD8 Cytotoxic",
                 "1" = "Naive",
                 "2" = "CD4 Helper",
                 "3" = "CD8 CM + \u03b3\u03b4-1", # unicode for gamma-delta symbols
                 "4" = "\u03b3\u03b4-17", 
                 "5" = "CD8 Activated",
                 "6" = "Proliferating")
scRNA <- RenameIdents(scRNA, clusters)
scRNA$clusters <- Idents(scRNA)

DimPlot(scRNA, reduction = "umap", shuffle = TRUE, label = TRUE) + NoLegend()

# Save integrated/annotated object
saveRDS(scRNA, file_seurat_analyzed_allcells)


##### Find cluster markers #####

DefaultAssay(scRNA) <- "RNA"
scRNA <- NormalizeData(scRNA)

all.markers <- FindAllMarkers(scRNA, test.use = "MAST", 
                              logfc.threshold = 0.25, min.pct = 0.1,
                              latent.vars = c("nCount_RNA", "StressScoreVariable"))

saveRDS(all.markers, file = file_markers_combined_all)


# Print most expressed genes in each cluster
highest.expressed <- getHighestExpressedGenes(scRNA)
highest.expressed

# Read in clonotype data
tcr.anno <- read.csv(file_clonotypes_processed)
table(tcr.anno$Genotype)

tcr.match <- subset(tcr.anno, Seurat.Barcode %in% colnames(scRNA))
table(tcr.match$Genotype)

# Get all diff genes in each cluster
FDR = 0.01
sig.markers = filter(all.markers, p_val_adj <= FDR & pct.1 >= 0.1)
table(sig.markers$cluster)

writeDifferentialGenes(sig.markers, file_markers_combined_clusters)

top10 <- printTop10Markers(sig.markers, pos.only = TRUE)

# Which clonotypes are in each cluster
writeClusterClonotypes(scRNA, tcr.anno, file_clonotypes_summary,
                       file_clonotypes_combined_clusters)

# Diff genes and clonotypes by genotype
genotypes <- unique(scRNA$genotype)
writeGenotypeDifferentialGenes(scRNA, genotypes, FDR, file_markers_combined_genotypes)
writeGenotypeClonotypes(scRNA, tcr.anno, genotypes, file_clonotypes_summary,
                        file_clonotypes_combined_genotypes)

writeClusterVGenotypeDiffGenes(scRNA, genotypes, FDR, file_markers_combined_clusters_vs_genotype)


##### AD risk genes - by cluster #####

ad.genes <- read.csv(file.path(dir_data_external, "20220214 Tanzi AD Gene List.csv"))

homologs <- getHumanMouseHomologs(unique(ad.genes$GENE))
ms.genes <- unique(homologs$MGI.symbol)
ms.genes <- intersect(ms.genes, rownames(scRNA))
ms.genes <- ms.genes[order(ms.genes)]

write_xlsx(homologs, path = file.path(dir_data, "AD_Gene_Homologs.xlsx"))

Idents(scRNA) <- scRNA$clusters

diff.genes <- lapply(excel_sheets(file_markers_combined_clusters), read_excel, 
                     path = file_markers_combined_clusters)

names(diff.genes) <- excel_sheets(file_markers_combined_clusters)

sig.genes <- lapply(diff.genes, function(D) {
  match <- intersect(D$Gene, ms.genes)
  filter(D, Gene %in% match)
})

write_xlsx(sig.genes, path = file_ad_risk_combined)


##### DPA #####

# This takes a long time with 10,000 iterations. 
Idents(scRNA) <- scRNA$clusters
runDPA(scRNA, file_dpa_allcells_glm_summary, file_dpa_allcells_pairwise,
       file.path(dir_allcells_dpa, "res_intermediate.rds"),
       file.path(dir_allcells_dpa, "fit_intermediate.rds"),
       file.path(dir_allcells_dpa, "pairwise_intermediate.rds"),
       nIter = 10000)


##### Everything beyond this point was not used in the paper #####


##### Identify naive/CM/EM/Eff cells like in flow cytometry. Not used for paper #####

DefaultAssay(scRNA) <- "RNA"
scRNA <- NormalizeData(scRNA)

scRNA <- UCell::AddModuleScore_UCell(scRNA, features = list("Naive" = c("Sell", "Ccr7", "Il7r", "Cd44-", "Klrg1-", "Itgal-", "Cd69-", "Cx3cr1-", "Il2rb-"),
                                                            "CM" = c("Sell", "Ccr7", "Il7r", "Cd44", "Klrg1-", "Itgal", "Cd69-", "Cx3cr1-", "Il2rb"),
                                                            "EM" = c("Sell-", "Ccr7-", "Il7r", "Cd44", "Klrg1", "Itgal", "Cd69-", "Cx3cr1", "Il2ra-"),
                                                            "RM" = c("Sell-", "Ccr7-", "Il7r", "Cd44", "Klrg1-", "Itgal", "Cd69", "Cx3cr1-"),
                                                            "Effector" = c("Sell-", "Ccr7-", "Il7r-", "Cd44-", "Klrg1", "Itgal-", "Cd69", "Cx3cr1", "Il2ra")),
                                     maxRank = 3000)

FeaturePlot(scRNA, features = c("Naive_UCell", "CM_UCell", "EM_UCell", "RM_UCell", "Effector_UCell"),
            cols = c("purple", "yellow"), ncol = 3)

FeaturePlot(scRNA, features = c("Naive_UCell", "CM_UCell", "EM_UCell", "RM_UCell", "Effector_UCell"),
            cols = c("purple", "yellow"), order = TRUE, ncol = 3)

types <- paste0(c("Naive", "CM", "EM", "RM", "Effector"), "_UCell")
scores <- scRNA@meta.data[,types]
colnames(scores) <- str_replace(colnames(scores), "_UCell", "")

tmp <- apply(scores, 1, which.max)
scRNA <- AddMetaData(scRNA, colnames(scores)[tmp], col.name = "Assignment")

Idents(scRNA) <- scRNA$Assignment
DimPlot(scRNA, shuffle = TRUE)
DimPlot(scRNA, split.by = "Assignment")

table(scRNA$Assignment)
table(scRNA$Assignment) / ncol(scRNA) * 100


##### GO Analysis - By genotype. Not used for paper. #####

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

# Done
