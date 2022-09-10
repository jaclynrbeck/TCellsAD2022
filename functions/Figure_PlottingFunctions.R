# Plotting functions commonly used for multiple figures
# Author: Jaclyn Beck
# Final script used for paper as of Sep 09, 2022

library(ggplot2)
library(dplyr)
library(viridis)
library(reshape2)
library(stringr)
library(ggprism)

# Creates a dot plot of the top 5 most upregulated markers for each cluster
# scRNA: Seurat object
# markers_file: full file path to an RDS file that contains the output of FindAllMarkers
# excluded genes: a list of sets of genes to exclude from display
# FDR: p value cutoff for significance
dotplotTop5Markers <- function( scRNA, markers_file, excluded_genes, FDR = 0.01 ) {
  all.markers <- readRDS(markers_file)
  sig.markers <- subset(all.markers, p_val_adj <= FDR & pct.1 >= 0.1)
  sig.markers <- subset(sig.markers, !(gene %in% unlist(excluded_genes)))
  
  markers.to.plot <- sig.markers %>% group_by(cluster) %>%
    top_n(5, wt = avg_log2FC)
  
  plt <- doViridisDotPlot(scRNA, unique(markers.to.plot$gene))
}


# Creates a dotplot with specific viridis colors and a few other tweaks
# scRNA: Seurat object
# feats: vector of gene names to include in the plot
doViridisDotPlot <- function( scRNA, feats ) {
  plt <- DotPlot(scRNA, features = feats, dot.scale = 3) + 
    RotatedAxis() + coord_flip() + 
    scale_color_viridis(option = "A", direction=-1, begin=0.6) +
    theme(axis.text=element_text(size=8), axis.title = element_blank(),
          plot.background = element_rect(fill = "white"))
}


# Creates a bar graph of population distributions, with the 4 genotypes
# side by side for each cluster.
# scRNA: Seurat object
# geno.colors: vector of colors to use for genotype bars. Should be in the same
#              order as levels(scRNA$genotype)
# sig.df: a data frame containing details about which comparisons are 
#         significant and where to draw significance bars and stars. Can be
#         null if no significance. See "Figure5.R" for an example format, but
#         briefly it should have the following fields: 
#           cluster: cluster name, repeated to fill out the column
#           group1: the first genotype in each comparison
#           group2: the second genotype in each comparison, paired with group1
#           p.adj: the p value for each comparison
#           xmin: x coordinate of the left side of each significance bar
#           xmax: x coordinate of the right side of each significance bar
#           y.position: y coordinate of each significance bar
populationBarGraph <- function( scRNA, geno.colors, sig.df = NULL) {
  pop <- table(Idents(scRNA), scRNA$orig.ident)
  cell.counts <- colSums(pop)
  cell.counts <- do.call(rbind, lapply(1:nrow(pop), function(x) cell.counts))
  pct <- melt(pop / cell.counts)
  
  pop <- melt(pop)
  pop <- cbind(pop, pct$value*100)
  colnames(pop) <- c("Cluster", "Sample", "Count", "Percent")
  pop$Genotype <- str_replace(pop$Sample, "-[12]", "")
  
  # Just for ordering of clusters in the graph
  means = aggregate(pop, list(Clust = pop$Cluster), mean)
  pop$Cluster <- factor(pop$Cluster, 
                        levels = means$Clust[order(means$Percent, decreasing = TRUE)])
  pop$Genotype <- factor(pop$Genotype, levels = c("WT", "5XFAD", "PS19", "PS-5X"))
  
  # Used to make error bars for the genotypes w/ 2 samples
  data <- pop %>% group_by(Cluster, Genotype) %>% 
    summarize(Pct = mean(Percent), SD = sd(Percent), Max = max(Percent), Min = min(Percent))
  data$Max[is.na(data$SD)] <- NA
  data$Min[is.na(data$SD)] <- NA
  
  # Part 1 of plot: bar graph with separate bars for each genotype, grouped
  # by cluster
  plt <- ggplot(data, aes(x = Cluster, y = Pct)) +
    geom_bar(aes(fill = Genotype), stat = "identity", color = "black", 
             position = "dodge", width = 0.6,
             size = 0.2) + 
    theme_classic() +
    xlab("Cluster") + ylab("Percent of Cell Population") + 
    scale_fill_discrete(type = geno.colors) +
    scale_x_discrete(guide = guide_axis(angle = 45)) + 
    theme(axis.text.x = element_text(size=8), axis.title.x = element_blank()) + 
    geom_errorbar(aes(group = Genotype, ymin = Min, ymax = Max), 
                  position = position_dodge(width = 0.6), width = 0.25, size = 0.2)
  
  # Part 2: Add significance bars
  if (!is.null(sig.df)) {
    plt <- plt + add_pvalue(sig.df, label = "p.adj.signif", label.size = 2.5,
                            xmin = "xmin", xmax = "xmax", 
                            tip.length = 0.01, bracket.size = 0.3, 
                            bracket.shorten = 0.04)
  }
  
  plt 
}


# Downsamples each genotype to have the same number of cells. Returns a new
# Seurat object with the downsampled data.
downsample <- function( scRNA ) {
  counts <- table(scRNA$genotype)
  downsamp <- min(counts)
  cells <- c()
  
  set.seed(10)
  for (N in names(counts)) {
    subs <- filter(scRNA@meta.data, genotype == N)
    cells <- c(cells, sample(rownames(subs), downsamp))
  }
  
  scRNA.down <- subset(scRNA, cells = cells)
  scRNA.down
}


# Draws UMAPs split by genotype, colored by cluster, where the cells have been
# downsampled to be equal across genotypes
umapSplitDownsampled <- function( scRNA, ident.name = "clusters", clust.colors ) {
  scRNA.down <- downsample(scRNA)
  Idents(scRNA.down) <- scRNA.down@meta.data[,ident.name]
  
  plt <- DimPlot(scRNA.down, reduction = "umap", split.by = "genotype", 
                 ncol = 2, pt.size = 0.01) + 
    NoLegend() + scale_color_discrete(type = clust.colors)
  
  rm(scRNA.down)
  
  plt
}


# Finds the top 10 most expanded clonotypes, by percentage of sample, then
# gives scRNA a new metadata field where cells with those clonotypes have
# the clonotype name assigned, and cells without those clonotypes have "N/A"
assignTop10Clonos <- function ( scRNA ) {
  # Get the top 10 most expanded clonotypes, by percentage of sample
  tcr.anno <- read.csv(file_clonotypes_processed)
  tcr.match <- subset(tcr.anno, Seurat.Barcode %in% colnames(scRNA))
  
  samp_counts <- table(tcr.match$ClonotypeId, tcr.match$Sample)
  samp_counts <- melt(samp_counts)
  colnames(samp_counts) <- c("Clonotype", "Sample", "Count")
  
  samp_counts <- subset(samp_counts, Count > 0)
  samp_counts$Genotype <- str_replace(samp_counts$Sample, "-[12]", "")
  
  # Calculates clonotype percent for each individual sample
  samp_counts <- samp_counts %>% group_by(Sample) %>% 
    mutate(Percent = 100 * Count / sum(Count)) %>% arrange(desc(Percent))
  
  top10 <- top_n(ungroup(samp_counts), n = 10, wt = Percent)
  
  scRNA[["Clonotype"]] <- "N/A"
  
  # Adds padding so single digit vs double digit clonotype IDs line up
  top10$ClonoName <- as.character(top10$Clonotype)
  short <- sapply(top10$ClonoName, nchar)
  top10$ClonoName[short < 11] <- str_pad(top10$ClonoName[short < 11], 
                                         width = 12, side = "right")
  top10$ClonoName <- paste0(top10$ClonoName, " (", 
                            format(round(top10$Percent, digits = 1), nsmall = 1), 
                            "% of ", top10$Sample,  ")")
  
  # Assigns the top 10 clonotype IDs to cells that have those clonotypes
  for (R in 1:nrow(top10)) {
    cells <- filter(tcr.match, ClonotypeId == top10$Clonotype[R])
    scRNA@meta.data[cells$Seurat.Barcode, "Clonotype"] <- top10$ClonoName[R]
  }
  
  scRNA$Clonotype <- factor(scRNA$Clonotype, 
                            levels = c(as.character(top10$ClonoName), "N/A"))
  
  scRNA
}


# Color-blind safe palette: light gray, orange, light blue, green, 
# light purple, dark blue, red, pink, light green, teal, dark purple
clono.colors <- c("#DDDDDD", "#E69F00", "#56B4E9", "#009E73", "#A461B4FF", 
                  "#0072B2", "#D55E00", "#CC79A7", 
                  viridis(3, end = 0.8, direction = -1))

# Draws a UMAP with the locations of the cells with the top 10 clonotypes.
umapTop10Clonotypes <- function( scRNA ) {
  Idents(scRNA) <- scRNA$Clonotype
  plt <- DimPlot(scRNA, reduction = "umap", pt.size = 0.1, shuffle = TRUE,
                 order = c(rev(levels(scRNA$Clonotype)[1:10]), "N/A")) +  
    scale_color_discrete(type = clono.colors) 
}


# Draws UMAPs with the locations of the cells with the top 10 clonotypes, 
# split by genotype.
# ncol: how many columns per row. For 4 genotypes, ncol = 2 gives a 2x2 grid,
#       ncol = 4 gives a 1x4 grid.
umapTop10ClonotypesSplit <- function( scRNA, ncol = 2 ) {
  Idents(scRNA) <- scRNA$Clonotype
  plt <- DimPlot(scRNA, reduction = "umap", pt.size = 0.1,
                 split.by = "genotype", ncol=ncol,
                 order = c(rev(levels(scRNA$Clonotype)[1:10]), "N/A")) +  
    scale_color_discrete(type = clono.colors)
}


# Draws UMAPs of gene expression with 5 columns per row, using the viridis
# color scheme and a few visual tweaks
# scRNA: Seurat object
# genes: vector of genes to make UMAPs for
# coord.fixed: whether to plot so the x and y axes of the UMAP have the same
#              scale, or not.
umapGeneExprPlots <- function( scRNA, genes, coord.fixed = TRUE ) {
  plt <- FeaturePlot(scRNA, features = genes,
                     order = TRUE, ncol = 5, pt.size = 0.01, 
                     coord.fixed = coord.fixed) &
    scale_color_viridis(begin = 0.1, option = "D") & 
    theme(axis.title = element_blank(), axis.text = element_blank(),
          title = element_text(size = 10), legend.position = "none",
          axis.line = element_blank(), axis.ticks = element_blank())
}


# Projects clusters from subsets (i.e. CD8 only) onto the UMAP of all cells
# together. Cells in the subset will be colored their assigned cluster color,
# and cells not in the subset will be gray.
# scRNA: Seurat object
# clust.colors: vector of colors to use for clusters. Should be in the same 
#               order as levels(Idents(scRNA))
umapProjectClustsToAllCells <- function( scRNA, clust.colors ) {
  scRNA.all <- readRDS(file_seurat_analyzed_allcells)
  
  DefaultAssay(scRNA.all) <- "RNA"
  scRNA.all <- DietSeurat(scRNA.all, assays = c("RNA"), dimreducs = c("pca", "umap"))
  gc()
  
  scRNA.all$NewClusts <- "N/A"
  scRNA.all$NewClusts[colnames(scRNA)] <- as.character(Idents(scRNA)[colnames(scRNA)])
  scRNA.all$NewClusts <- factor(scRNA.all$NewClusts, levels = c("N/A", levels(Idents(scRNA))))
  
  Idents(scRNA.all) <- scRNA.all$NewClusts
  
  new.colors <- c("#DDDDDDFF", clust.colors)
  plt <- DimPlot(scRNA.all, reduction = "umap", shuffle = FALSE, label = FALSE,
                 pt.size = 0.01, order = rev(levels(scRNA.all$NewClusts))) + 
    scale_color_discrete(type = new.colors) + 
    theme(axis.line = element_blank(), axis.text = element_blank(), 
          axis.title = element_blank(), axis.ticks = element_blank())
  
  plt
}


