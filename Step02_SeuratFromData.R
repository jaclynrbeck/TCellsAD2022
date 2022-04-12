# Reads 10X data from the output of cellranger aggr and turns it into a Seurat
# object that has been normalized and scaled. 

# This script contains multiple methods for creating the Seurat object from the
# matrix, but ultimately only the integration/SCT method is used for downstream 
# analysis.

# Author: Jaclyn Beck
# Final script used for paper as of April 11, 2022

library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
source(file.path("functions", "SeuratFromData_HelperFunctions.R"))
source(file.path("functions", "Analysis_HelperFunctions.R"))
source("Filenames.R")

# Look at UMIs vs cell counts
data <- Read10X(dir_filtered_counts, gene.column=2)

plotUMICurve(data, combine_types=FALSE)
plotUMICurve(data, combine_types=TRUE)

rm(data)

# Create the Seurat object
scRNA.all <- makeSeurat(dir_filtered_counts)
scRNA <- scRNA.all

head(scRNA@meta.data, 5)
VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 
                            "percent.rb", "complexity"), ncol=3, pt.size = 0)

# Identify potential doublets from TCR data -- looking for multiple beta-chains
# only. Some T cells can express two alpha chains so we ignore that. 
tcr.anno <- read.csv(file_clonotypes_processed)
table(tcr.anno$Genotype)

genes <- sapply(c("TRB.V", "TRB.D", "TRB.J", "TRB.C"),
                function(X) { str_split(tcr.anno[,X], pattern = ", ") })

tcr.anno$Singlet <- sapply(1:nrow(genes), function(X) {
  lengths <- sapply(genes[X,], length)
  all(lengths <= 1)
})
(table(tcr.anno$Singlet) / nrow(tcr.anno)) * 100 # Doublet rate is ~4%

scRNA$Singlet <- "Unknown"
scRNA$Singlet[tcr.anno$Sample] <- tcr.anno$Singlet
(table(scRNA$Singlet) / ncol(scRNA)) * 100 

# Only take Cd3+ cells or cells with a recognized TCR, remove doublets
pos.cells <- filterOnGenePositive(scRNA, "Cd3[edg]$", 1)
scRNA <- subset(scRNA, cells = unique(c(pos.cells, tcr.anno$Sample)))
scRNA <- subset(scRNA, Singlet != FALSE)

scRNA <- removeLowExpressedGenes(scRNA, 10)

# Plots to examine the quality of the remaining data

VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 
                            "percent.rb", "complexity"), ncol=3, pt.size = 0)

FeatureScatter(scRNA, feature1 = "nFeature_RNA", feature2 = "nCount_RNA", shuffle=TRUE, pt.size = 0.1) +
  scale_y_log10(n.breaks=16) + scale_x_log10(n.breaks=16) + theme_bw()
FeatureScatter(scRNA, feature1 = "complexity", feature2 = "nCount_RNA", shuffle=TRUE, pt.size = 0.1) +
  scale_y_log10(n.breaks=16) + scale_x_continuous(n.breaks=16) + theme_bw()
FeatureScatter(scRNA, feature1 = "nFeature_RNA", feature2 = "percent.mt", shuffle=TRUE, pt.size = 0.1) +
  scale_y_log10(n.breaks=16) + scale_x_log10(n.breaks=16) + theme_bw()
FeatureScatter(scRNA, feature1 = "nFeature_RNA", feature2 = "percent.rb", shuffle=TRUE, pt.size = 0.1) +
  scale_y_continuous(n.breaks=16) + scale_x_log10(n.breaks=16) + theme_bw()
FeatureScatter(scRNA, feature1 = "percent.rb", feature2 = "percent.mt", shuffle=TRUE, pt.size = 0.1) +
  scale_y_log10(n.breaks=16) + scale_x_continuous(n.breaks=16) + theme_bw()


ggplot(data.frame(y=sort(scRNA[["percent.mt"]]$percent.mt), x=1:length(scRNA[["percent.mt"]]$percent.mt)), aes(x,y)) + 
  geom_line() + scale_y_log10(n.breaks=16)
ggplot(data.frame(y=sort(scRNA[["percent.rb"]]$percent.rb), x=1:length(scRNA[["percent.rb"]]$percent.rb)), aes(x,y)) + 
  geom_line() + scale_y_log10(n.breaks=16)
ggplot(data.frame(y=sort(scRNA[["nFeature_RNA"]]$nFeature_RNA), x=1:length(scRNA[["nFeature_RNA"]]$nFeature_RNA)), aes(x,y)) + 
  geom_line() + scale_y_log10(n.breaks=16) + geom_abline(slope = 0, intercept = log10(2*median(scRNA$nFeature_RNA)), color = "red")
ggplot(data.frame(y=sort(scRNA[["nCount_RNA"]]$nCount_RNA), x=1:length(scRNA[["nCount_RNA"]]$nCount_RNA)), aes(x,y)) + 
  geom_line() + scale_y_log10(n.breaks=16) + geom_abline(slope = 0, intercept = log10(2*median(scRNA$nCount_RNA)), color = "red")

median(scRNA$nFeature_RNA)
mean(scRNA$nFeature_RNA)
sd(scRNA$nFeature_RNA)
median(scRNA$nCount_RNA)
mean(scRNA$nCount_RNA)
sd(scRNA$nCount_RNA)

# Get rid of the set of low-quality cells
scRNA <- subset(scRNA, subset = nFeature_RNA >= 800 & 
                  percent.mt <= 5 &   # Clear cutoff at 4 or 5
                  percent.rb >= 8 &   # Clear cutoff at 8-10
                  percent.rb <= 55)   # Removes cells with low RNA + high rb

VlnPlot(scRNA, features = "nFeature_RNA", pt.size = 0) + 
  geom_hline(yintercept = 2*median(scRNA$nFeature_RNA))

VlnPlot(scRNA, features = "nCount_RNA", pt.size = 0.1) + 
  geom_hline(yintercept = 1000) +
  geom_hline(yintercept = 2*median(scRNA$nCount_RNA)) +
  geom_hline(yintercept = 3*median(scRNA$nCount_RNA))

# Remove potential doublets by thresholding RNA count
scRNA <- subset(scRNA, nCount_RNA <= 3*median(scRNA$nCount_RNA))

VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 
                            "percent.rb", "complexity"), ncol=3, pt.size = 0)

table(scRNA$orig.ident)
scRNA <- removeLowExpressedGenes(scRNA, 10)
scRNA <- DietSeurat(scRNA, scale.data = FALSE)
scRNA$genotype <- str_replace(scRNA$orig.ident, "-[12]", "")
scRNA$batch <- "Batch1"
scRNA$batch[scRNA$orig.ident == "WT-2" | scRNA$orig.ident == "5XFAD-2"] <- "Batch2"

# At this point scRNA contains all cells that passed QC
saveRDS(scRNA, file = file_seurat_unnorm)

# Looking for batch effect of cellular stress
scRNA <- readRDS(file_seurat_unnorm)
scRNA <- NormalizeData(scRNA)
scRNA <- FindVariableFeatures(scRNA, nfeatures = 3000)

cells = colnames(scRNA)[scRNA$genotype == "WT" | scRNA$genotype == "5XFAD"]

scRNA.sub <- subset(scRNA, cells = cells)
markers <- FindMarkers(scRNA.sub, ident.1 = "Batch2", ident.2 = "Batch1",
                       group.by = "batch", min.pct = 0.05, test.use = "MAST",
                       latent.vars = c("nCount_RNA", "genotype"))

markers <- markers[order(markers$avg_log2FC, decreasing = TRUE),]
markers <- subset(markers, p_val_adj < 0.01)
markers$Ensembl.Id <- geneNameToEnsembl(rownames(markers))

write.csv(markers, file = "DiffGenes_Batch2v1.csv")

rm(scRNA.sub)

# After input of diff genes into PANTHER, scoring on Reactome pathways, the top 
# results all relate to cellular stress. 9 genes overlap with the pathway
# "Cellular response to heat stress", so we use genes in that pathway to
# create a "stress score"
reac <- read.table("~/Downloads/reactome_ResponseToHeatStress.tsv", sep = "\t",
                   skip = 4, header = TRUE)

reac_stress <- reac$Gene.Name[reac$Gene.Name %in% rownames(scRNA)]

assay <- GetAssayData(scRNA, slot = "data")
assay <- assay[reac_stress,]
assay <- rowSums(assay > 0)

reac_stress <- reac_stress[assay > ncol(scRNA) * 0.01]

# Only score on genes that are variable, to help discriminate between clusters
reac_v <- reac_stress[reac_stress %in% VariableFeatures(scRNA)]

scRNA <- AddModuleScore(scRNA, features = list(reac_stress, reac_v), name = "Stress")

RidgePlot(scRNA, features = c("Stress1", "Stress2"), 
          group.by = "orig.ident", same.y.lims = TRUE)

# This data is used for the paper 
#scRNA <- readRDS(file_seurat_unnorm)
scRNA <- generateIntegratedData(scRNA, method = "LogNormalize",
                                regress = c("nCount_RNA", "Stress2", "percent.mt", "percent.rb"))

saveRDS(scRNA, file = file_seurat_norm)

# At this point we are done. The rest of this script is just more data 
# exploration and visualization for debugging


# Looking at correlations between metadata variables, including the first 10 PC's

scRNA <- readRDS(file_seurat_norm)
scRNA <- RunPCA(scRNA, npcs = 75)
ElbowPlot(scRNA, ndims=75)
DimPlot(scRNA, reduction = "pca", group.by = "orig.ident", shuffle = TRUE)
DimPlot(scRNA, reduction = "pca", group.by = "orig.ident", 
        split.by = "orig.ident", shuffle = TRUE)

scRNA <- FindNeighbors(scRNA, reduction="pca", dims = 1:25)
scRNA <- FindClusters(scRNA, resolution = 0.2)

scRNA <- RunUMAP(scRNA, dims = 1:25)
DimPlot(scRNA, reduction = "umap", shuffle=TRUE, label=TRUE)
DimPlot(scRNA, group.by="orig.ident", reduction = "umap", shuffle=TRUE)
DimPlot(scRNA, group.by="orig.ident", split.by="orig.ident", reduction = "umap", shuffle=TRUE)

# Cell Cycle
cc_genes <- getCellCycleGenes()
scRNA <- CellCycleScoring(scRNA,
                          g2m.features = cc_genes$g2m_genes,
                          s.features = cc_genes$s_genes)

DimPlot(scRNA, group.by="Phase", reduction = "umap", shuffle=TRUE)
DimPlot(scRNA, group.by="Phase", split.by="Phase", reduction = "umap", shuffle=TRUE)

FeatureScatter(scRNA, feature1 = "G2M.Score", feature2 = "nCount_RNA", shuffle = TRUE, pt.size = 0.1) + 
  scale_y_log10(n.breaks=16) + scale_x_continuous(n.breaks=16) + theme_bw()

good.cells <- scRNA$nCount_RNA <= 2*median(scRNA$nCount_RNA)
prolif.cells <- scRNA$G2M.Score > 0 & scRNA$nCount_RNA <= 3*median(scRNA$nCount_RNA)

scRNA <- subset(scRNA, cells = colnames(scRNA)[good.cells | prolif.cells])


library(pheatmap)

head(Embeddings(scRNA, reduction="pca")[,1:10])

merged <- cbind(scRNA@meta.data, 
                Embeddings(scRNA, reduction="pca")[rownames(scRNA@meta.data),1:10])
merged$orig.ident <- as.numeric(factor(merged$orig.ident))
merged$integrated_snn_res.0.2 <- as.numeric(merged$integrated_snn_res.0.2)
merged$RNA_snn_res.1 <- as.numeric(merged$RNA_snn_res.1)
merged$seurat_clusters = as.numeric(merged$seurat_clusters)
merged$Phase <- as.numeric(factor(merged$Phase))

correlation = cor(merged)
map = pheatmap::pheatmap(correlation, display_numbers = TRUE)

# Visualize the PCA, grouping by cell cycle phase
DimPlot(scRNA, reduction = "pca", group.by = "Phase", shuffle = TRUE)

DimPlot(scRNA, reduction = "pca", group.by = "Phase", split.by = "Phase")

DimPlot(scRNA, reduction = "pca", group.by = "orig.ident", 
        split.by = "orig.ident", shuffle = TRUE)

FeatureScatter(scRNA, feature1 = "PC_1", feature2 = "PC_2", 
               shuffle=TRUE, raster=TRUE) 
FeatureScatter(scRNA, feature1 = "PC_1", feature2 = "S.Score", 
               shuffle=TRUE, raster=TRUE)
FeatureScatter(scRNA, feature1 = "PC_1", feature2 = "G2M.Score", 
               shuffle=TRUE, raster=TRUE)
FeatureScatter(scRNA, feature1 = "PC_2", feature2 = "nCount_RNA", 
               shuffle=TRUE, raster=TRUE) + scale_y_log10()
FeatureScatter(scRNA, feature1 = "PC_2", feature2 = "nFeature_RNA", 
               shuffle=TRUE, raster=TRUE) + scale_y_log10()
FeatureScatter(scRNA, feature1 = "PC_1", feature2 = "complexity", 
               shuffle=TRUE, raster=TRUE)

DimPlot(scRNA, dims = c(1,2), reduction = "pca", shuffle = TRUE)
DimHeatmap(scRNA, dims = c(1,2), cells = 500, balanced = TRUE)
ElbowPlot(scRNA, ndims=75)

rm(list = ls())

# Test code: DoubletFinder
library(DoubletFinder)
splitS <- SplitObject(scRNA, split.by = "orig.ident")
sweep.res <- lapply(splitS, FUN = paramSweep_v3, PCs = 1:30, sct = FALSE)
sweep.stats <- lapply(sweep.res, FUN = summarizeSweep, GT = FALSE)
bcmvn <- lapply(sweep.stats, FUN = find.pK)
est_pk <- lapply(bcmvn, function(X) { X$pK[which.max(X$BCmetric)] })

nExp_poi <- lapply(splitS, function(X) { round(0.05*ncol(X)) })

splitS <- lapply(1:length(splitS), function(X) {
  doubletFinder_v3(splitS[[X]], PCs = 1:30, pN = 0.25, pK = as.numeric(as.character(est_pk[[X]])), nExp = nExp_poi[[X]], reuse.pANN = FALSE, sct = FALSE)
})

splitS <- lapply(splitS, function(X) {
  colnames(X@meta.data)[grep("DF.class", colnames(X@meta.data))] = "DF.Class"
  colnames(X@meta.data)[grep("pANN", colnames(X@meta.data))] = "DF.pANN"
  X
})

tmp <- merge(splitS[[1]], unlist(splitS[2:6]))
