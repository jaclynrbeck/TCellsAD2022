# This script will create figures 5B-5H from the paper and save them to files.
# Author: Jaclyn Beck
# Final script used for paper as of Sep 09, 2022

library(Seurat)
library(readxl)
library(pheatmap)
library(Matrix)
source(file.path("functions", "General_HelperFunctions.R"))
source(file.path("functions", "Figure_PlottingFunctions.R"))
source("Filenames.R")

##### Colors #####

# WT (green), 5X (red), PS-19 (blue), PS-5X (purple)
geno.colors <- c("WT" = "#51B400", "5XFAD" = "#DE3923", 
                 "PS19" = "#4B5BDD", "PS-5X" = "#6B3D92") 

# Should be colorbind/black and white printing safe:
# Light red, mustard, light blue, green, magenta, purple, cyan
clust.colors <- c("#FF867D", "#E5B920", "#009BED", "#099600", "#DF588F", 
                  "#6560BF", "#00CBFB")


##### Setup #####

scRNA <- readRDS(file_seurat_analyzed_allcells)

DefaultAssay(scRNA) <- "RNA"

scRNA$genotype <- factor(scRNA$genotype, 
                         levels = c("WT", "5XFAD", "PS19", "PS-5X"))

excluded_genes <- readRDS(file_excluded_genes)

# We don't need the integrated assay
scRNA <- DietSeurat(scRNA, assays = c("RNA"), dimreducs = c("pca", "umap"))
gc()

Idents(scRNA) <- scRNA$clusters


##### Figure 5B #####

# UMAP with labels in the legend -- labels manually placed afterward in figure
plt <- DimPlot(scRNA, reduction = "umap", shuffle = TRUE, label = FALSE) +
  scale_color_discrete(type = clust.colors)

plt

ggsave(file.path(dir_figures_paper, "fig5b.png"),
       plot = plt, width = 6.1, height = 3.8, units = "in", dpi = "print")


##### Figure 5C #####

# Dotplot of top 5 upregulated genes by cluster
plt <- dotplotTop5Markers(scRNA, file_markers_combined_all, excluded_genes, 
                          FDR = 0.01)
plt

ggsave(file.path(dir_figures_paper, "fig5c.png"),
       plot = plt, width = 5.25, height = 5.25, units = "in", dpi = "print")


##### Figure 5D #####

# UMAP of clusters split by genotype, downsampled to have even cell numbers
plt <- umapSplitDownsampled(scRNA, ident.name = "clusters", clust.colors)
plt

ggsave(file.path(dir_figures_paper, "fig5d.png"),
       plot = plt, width = 4.5, height = 4.25, units = "in", dpi = "print")


##### Figure 5E #####

# Population bar graph. Significance values are from DPA analysis.
# A lot of values in this section are hard-coded and would need to be changed
# if the data changes. 

# centers of each bar relative to the center of the group (-0.33 to +0.33)
bars <- c("WT" = -0.23, "5XFAD" = -0.08, "PS19" = 0.08, "PS-5X" = 0.23)
sig.df <- rbind(
  data.frame(cluster = rep("CD8 Cytotoxic", 4), 
             group1 = c("WT",    "WT",    "5XFAD", "PS19"), 
             group2 = c("5XFAD", "PS-5X", "PS19",  "PS-5X"),
             p.adj = c(0.01062,  0.00834, 0.00542, 0.00542),
             xmin = 1 + c(bars["WT"], bars["WT"], bars["5XFAD"], bars["PS19"]),
             xmax = 1 + c(bars["5XFAD"], bars["PS-5X"], bars["PS19"], bars["PS-5X"]), 
             y.position = c(59.5, 62.5, 59.5, 59.5)),
  data.frame(cluster = rep("Naive", 3), 
             group1 = c("WT",    "5XFAD", "PS19"), 
             group2 = c("5XFAD", "PS19",  "PS-5X"),
             p.adj = c(0.01021,  0.00027, 0.00309),
             xmin = 2 + c(bars["WT"], bars["5XFAD"], bars["PS19"]),
             xmax = 2 + c(bars["5XFAD"], bars["PS19"], bars["PS-5X"]), 
             y.position = c(27.3, 33.3, 33.3)),
  data.frame(cluster = rep("CD8 Activated", 2), 
             group1 = c("WT",    "WT"), 
             group2 = c("5XFAD", "PS-5X"),
             p.adj = c(0.03486,  0.03486),
             xmin = 6 + c(bars["WT"], bars["WT"]),
             xmax = 6 + c(bars["5XFAD"],  bars["PS-5X"]), 
             y.position = c(11.75, 14.75))
) %>%
  rstatix::add_significance(p.col = "p.adj") 

plt <- populationBarGraph(scRNA, geno.colors, sig.df)
plt

ggsave(file.path(dir_figures_paper, "fig5e.png"),
       plot = plt, width = 5, height = 3, units = "in", dpi = "print")


##### Figure 5F #####

# Top 100 differential genes for each genotype vs. WT

# Get significant genes vs WT for all genotypes
sig.genes.df <- readSigGenesGenotype(file_markers_combined_genotypes, pattern = "vs WT")

# Remove mitochondrial and ribosomal genes plus batch effect genes for display.
sig.genes.df <- subset(sig.genes.df, !(Gene %in% c(unlist(excluded_genes), "Xist")))

sig.genes.df <- arrange(sig.genes.df, desc(abs(avg_log2FC)))
features <- unique(sig.genes.df$Gene)[1:100]

assay <- GetAssayData(scRNA, slot = "data", assay = "RNA")
assay <- assay[features,]

df <- melt(cbind(rownames(assay), as.data.frame(assay)))
colnames(df) <- c("Gene", "Cell", "Expression")
df$Genotype <- str_replace(df$Cell, ".*_", "")
df$Genotype <- str_replace(df$Genotype, "-[1|2]", "")

# Average expression per genotype
df.means <- aggregate(df$Expression, by = list(Gene = df$Gene, Genotype = df$Genotype), 
                      FUN = mean)
colnames(df.means) <- c("Gene", "Genotype", "Expression")

df.means <- dcast(df.means, Genotype ~ Gene)
rownames(df.means) <- df.means$Genotype
df.means <- t(df.means[,-1])

pheatmap(df.means, fontsize_row = 5, 
         treeheight_col = 0, treeheight_row = 0, 
         color = viridis(50, option = "D", begin = 0), border_color = NA,
         scale = "row", 
         filename = file.path(dir_figures_allcells, "allcells_heatmap_siggenes_genotype_average.png"),
         width = 2.2, height = 7.5)

if (dev.cur() != 1) {
  dev.off()
}


##### Figure 5G #####

# Heatmap of key cell type-specific markers
markers <- list(CellType = c("Cd8b1", "Cd4", "Trdc"),
                Type1 = c("Ifng", "Tnf"),
                Type2 = c("Il4", "Il13"),
                Type17 = c("Rorc", "Il17a"),
                Treg = c("Il2ra", "Foxp3"),
                NaiveCM = c("Tcf7", "Ccr7", "Sell"),
                Memory = c("Il7r", "Itgal", "Cd44"),
                CytotoxicEffector = c("Gzma", "Gzmb", "Gzmk", "Nkg7", "Prf1",
                                      "Ccl5", "Cx3cr1"),
                Residence = c("Itga1", "Itgae", "Cd69"),
                Exhaustion = c("Ctla4", "Pdcd1", "Lag3", "Tigit"),
                IFN = c("Ifit3", "Irf7", "Stat1"),
                Proliferation = c("Mki67", "Top2a")
)

assay <- GetAssayData(scRNA, slot = "data", assay = "RNA")
assay <- assay[unlist(markers),]

assay <- t(scale(t(assay))) # Center / scale
assay <- pmin(assay, 2.5)
assay <- pmax(assay, -2.5)

genes <- c()

# Cluster each category of genes by expression so similar genes are near 
# each other.
for (N in names(markers)) {
  d <- dist(assay[markers[[N]],])
  ord <- hclust(d)
  
  genes <- c(genes, markers[[N]][ord$order])
}

# Cluster each cluster of cells so similar cells are near each other, but
# remain separated by Seurat cluster
cells <- c()

for (C in levels(scRNA$clusters)) {
  meta <- subset(scRNA@meta.data, clusters == C)
  dist_cells <- dist(t(assay[,rownames(meta)]))
  order_cells <- hclust(dist_cells)
  
  cells <- c(cells, rownames(meta)[order_cells$order])
}

rm(dist_cells, order_cells)

assay <- assay[genes, cells]

anno <- data.frame(row.names = cells, Cluster = scRNA@meta.data[cells, "clusters"])
anno.cols <- clust.colors
names(anno.cols) <- levels(scRNA$clusters)
anno.cols <- list(Cluster = anno.cols)

anno.row <- data.frame(row.names = rownames(assay), Marker = stack(markers)$ind)

gaps_row <- cumsum(sapply(markers, length))

pheatmap(assay, fontsize_row = 8, cluster_cols = FALSE, cluster_rows = FALSE,
         treeheight_col = 0, treeheight_row = 0, 
         color = viridis(50, option = "D"),
         border_color = NA,
         annotation_col = anno, 
         #annotation_row = anno.row, 
         annotation_colors = anno.cols,
         annotation_names_col = FALSE,
         show_colnames = FALSE, gaps_row = gaps_row,
         filename = file.path(dir_figures_paper, "fig5g.png"),
         width = 8, height = 4.7)
if (dev.cur() != 1) {
  dev.off()
}


##### Figure 5H #####

# AD risk genes by cluster

Idents(scRNA) <- scRNA$clusters

# Read significant genes by cluster
sig.genes <- lapply(excel_sheets(file_ad_risk_combined), read_excel, 
                     path = file_ad_risk_combined)

names(sig.genes) <- excel_sheets(file_ad_risk_combined)

# Do some re-shaping
sig.genes <- Map(function(i, x) {
  x$cluster <- i
  x
}, names(sig.genes), sig.genes)

sig.genes.df <- do.call(rbind, sig.genes)

sig.genes.pos <- subset(sig.genes.df, avg_log2FC > 0)

plt <- doViridisDotPlot(scRNA, unique(sig.genes.pos$Gene))

plt

ggsave(file.path(dir_figures_paper, "fig5h.png"),
       plot = plt, width = 5.25, height = 5.5, units = "in", dpi = "print")


rm(list = ls())
gc()
