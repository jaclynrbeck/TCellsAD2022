# Reads 10X data from the output of cellranger aggr and turns it into a Seurat
# object that has been normalized and scaled. 

# This script contains multiple methods for creating the Seurat object from the
# matrix, but ultimately only the integration/SCT method is used for downstream 
# analysis.

# Author: Jaclyn Beck
# Final script used for paper as of Feb 03, 2022

library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(stringr)
source(file.path("functions", "SeuratFromData_HelperFunctions.R"))
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

# Find and remove bad clusters of cells
scRNA <- NormalizeData(scRNA)
scRNA <- FindVariableFeatures(scRNA)
scRNA <- ScaleData(scRNA, features = rownames(scRNA), model.use = "negbinom")

scRNA <- RunPCA(scRNA, npcs = 50, features = rownames(scRNA))
ElbowPlot(scRNA, ndims=50)
DimPlot(scRNA, reduction = "pca", group.by = "orig.ident", shuffle = TRUE)

scRNA <- FindNeighbors(scRNA, reduction="pca", dims = 1:10)
scRNA <- FindClusters(scRNA, resolution = 1.0)

scRNA <- RunUMAP(scRNA, dims = 1:10)
DimPlot(scRNA, reduction = "umap", shuffle=TRUE, label=TRUE)
DimPlot(scRNA, group.by="orig.ident", reduction = "umap", shuffle=TRUE)

VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 
                            "percent.rb", "complexity"), ncol=3, pt.size = 0)

VlnPlot(scRNA, features = "nFeature_RNA", pt.size = 0) + 
  geom_abline(slope = 0, intercept = 800) +
  geom_abline(slope = 0, intercept = 1000)

VlnPlot(scRNA, features = "nCount_RNA", pt.size = 0.1) + 
  geom_abline(slope = 0, intercept = 2*median(scRNA$nCount_RNA)) +
  geom_abline(slope = 0, intercept = 3*median(scRNA$nCount_RNA))

# Clusters 12, 16, and 19 have low complexity/RNA. 15 is borderline but we'll leave it in
clust.remove <- c(12, 16, 19)
clust.ok <- setdiff(unique(scRNA$seurat_clusters), clust.remove)
scRNA <- subset(scRNA, seurat_clusters %in% clust.ok)

Idents(scRNA) <- scRNA$orig.ident

# Plots to examine the quality of the remaining data

VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 
                            "percent.rb", "complexity"), ncol=3, pt.size = 0.1)

FeatureScatter(scRNA, feature1 = "nFeature_RNA", feature2 = "nCount_RNA", shuffle=TRUE, pt.size = 0.1) +
  scale_y_log10(n.breaks=16) + scale_x_log10(n.breaks=16) + theme_bw()
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

# Final QC filtering
scRNA <- subset(scRNA, subset = nCount_RNA <= 3*median(scRNA$nCount_RNA) &
                  nFeature_RNA >= 800 &   # Possible cutoff between ~800-1000
                  percent.mt <= 5 &   # Clear cutoff at 4 or 5
                  percent.rb >= 10 &  # Clear cutoff at 10
                  percent.rb <= 55)   # Removes cells with low RNA + high rb

VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 
                            "percent.rb", "complexity"), ncol=3, pt.size = 0)

table(scRNA[["orig.ident"]])
scRNA <- DietSeurat(scRNA, scale.data = FALSE)
scRNA$genotype <- str_replace(scRNA$orig.ident, "-[12]", "")

# At this point scRNA contains only cells with UMI >= 800
saveRDS(scRNA, file = file_seurat_unnorm)

# Experimenting with different methods of combining/normalizing data. 
# Commented out because they're unused. 
#scRNA <- readRDS(file_seurat_unnorm)
#scRNA <- simpleNormalizeAndScale(scRNA)
#saveRDS(scRNA, file = file.path(data_dir, "seurat_umi800_normSimple_2022-01-04.rds"))

#scRNA <- readRDS(file_seurat_unnorm)
#scRNA <- regressionNormalizeAndScale(scRNA, 
#                                     regress = c("percent.mt", "percent.rb"))
#saveRDS(scRNA, file = file.path(data_dir, "seurat_umi800_norm_regMtRb_2022-01-04.rds"))

#scRNA <- readRDS(file_seurat_unnorm)
#scRNA <- SCTransformNormalize(scRNA, regress = c())
#saveRDS(scRNA, file = file.path(data_dir, "seurat_umi800_sct_2022-01-04.rds"))

#scRNA <- readRDS(file_seurat_unnorm)
#scRNA <- SCTransformNormalize(scRNA, regress = c("percent.mt", "percent.rb"))
#saveRDS(scRNA, file = file.path(data_dir, "scRNA_Seurat_umi800_sct_regMtRb_2022-01-04.rds"))

#scRNA <- readRDS(file_seurat_unnorm)
#scRNA <- generateIntegratedData(scRNA)
#saveRDS(scRNA, file = file.path(data_dir, "seurat_umi800_integrated_2022-01-04.rds"))

# This data is used for the paper -- only includes cells with > 1 count of a 
# CD3 gene (Cd3d, Cd3e, or Cd3g)
scRNA <- readRDS(file_seurat_unnorm)
scRNA <- filterOnGenePositive(scRNA, "Cd3[edg]$", 1)
scRNA <- generateIntegratedData(scRNA, 
                                regress = c("percent.rb", "percent.mt", "nCount_RNA"))
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


# Perform cell cycle scoring
cc_genes <- getCellCycleGenes()
scRNA <- CellCycleScoring(scRNA,
                          g2m.features = cc_genes$g2m_genes,
                          s.features = cc_genes$s_genes)

DimPlot(scRNA, group.by="Phase", reduction = "umap", shuffle=TRUE)
DimPlot(scRNA, group.by="Phase", split.by="Phase", reduction = "umap", shuffle=TRUE)

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
map = pheatmap(correlation, display_numbers = TRUE)

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




