# Helper functions used by "Step02_SeuratFromData".

library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)

# Plots a cumulative line graph of UMI vs number of cells
# combine_types = TRUE: all genotypes are merged into one line
# combine_types = FALSE: each genotype gets its own line
plotUMICurve <- function(data, combine_types = FALSE) {
  barcodes <- data.frame(Barcode = colnames(data),
                         Sample = str_replace(colnames(data), ".*_", ""))
  
  umis <- data.frame(Cells=c(), UMIs=c(), Sample=c())
  
  if (combine_types == TRUE) {
    cs <- colSums(data > 0)
    umis <- data.frame(Cells=1:length(cs),
                       UMIs=sort(cs, decreasing=TRUE), 
                       Sample="Combined")
  }
  else {
    for (s in levels(factor(barcodes$Sample))) {
      samp_bcs <- subset(barcodes, Sample == s)
      cs <- colSums(data[, samp_bcs$Barcode] > 0)
      scs <- data.frame(Cells=1:length(cs), 
                        UMIs=sort(cs, decreasing=TRUE), 
                        Sample=s)
      
      umis <- rbind(umis, scs)
    }
  }
  
  plt = ggplot(umis, aes(x=Cells, y=UMIs, color=Sample)) + 
    geom_line() +
    scale_x_log10(n.breaks=12) + scale_y_log10(n.breaks=12) +
    ylab("UMIs") + theme_bw() +
    scale_color_brewer(palette="Set1")
  plt
}


# Create a Seurat object from 10X data. Adds percent mitochondrial genes and
# percent ribosomal genes as metadata. We immediately filter out any cells
# that have < 200 distinct genes counted, and each gene must have at least 10 
# cells that express it
makeSeurat <- function(file_dir, gene_column = 2, min_cells = 10, 
                       min_features = 200) {
  data <- Read10X(file_dir, gene.column=gene_column)
  scRNA <- CreateSeuratObject(data, 
                              project="T_Cells", 
                              min.cells = min_cells, 
                              min.features = min_features, 
                              names.field = 2, 
                              names.delim = "_")
  
  rm(data)
  
  # We are counting mitochondrial ribosomal proteins here. These patterns will 
  # recognize genes starting with "mt-" or "Mrps" or "Mrpl"
  scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^mt-")
  scRNA[["percent.mt"]] <- scRNA[["percent.mt"]] + PercentageFeatureSet(scRNA, pattern = "^Mrp[s|l]")
  
  # This pattern matches genes starting with "Rpl" or "Rps"
  scRNA[["percent.rb"]] <- PercentageFeatureSet(scRNA, pattern = "^Rp[s|l]")
  
  scRNA[["complexity"]] = scRNA[["nCount_RNA"]] / scRNA[["nFeature_RNA"]]

  scRNA
}


# Finds the confidence interval for nCount_RNA using its log-normal distribution.
# 2 sds gives a ~95% confidence interval, 3 gives ~99.7
getConfidenceInterval <- function(scRNA, sds = 2) {
  log.rna <- log10(scRNA$nCount_RNA)
  ci <- sd(log.rna)*sds
  high.cutoff <- 10^(mean(log.rna) + ci)
  low.cutoff <- 10^(mean(log.rna) - ci)
  
  list("high" = high.cutoff, "low" = low.cutoff)
}


# Returns a Seurat object that only contains cells that are "positive" for a 
# specific gene (i.e. RNA count > threshold). 
filterOnGenePositive <- function(scRNA, pattern, threshold = 1) {
  all.genes <- rownames(scRNA)
  
  assay <- GetAssayData(scRNA, "counts")
  assay.pos <- assay[grep(pattern, all.genes, value=TRUE),]
  
  if (class(assay.pos) != "numeric") {
    assay.pos <- colSums(assay.pos)
  }
  
  # Require >= threshold count
  pos.cells <- subset(names(assay.pos), assay.pos >= threshold) 
}

# Removes genes that are expressed in < (threshold) cells
removeLowExpressedGenes <- function(scRNA, threshold = 10) {
  assay <- GetAssayData(scRNA, slot = "counts")
  assay <- rowSums(assay > 0)
  genes <- assay[assay >= threshold]
  scRNA <- subset(scRNA, features = names(genes))
  scRNA
}


# Gets the set of mouse genes associated with cell cycle, used for Seurat's
# CellCycleScoring function
getCellCycleGenes <- function() {
  library(RCurl)
  
  cc_file <- getURL("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Mus_musculus.csv") 
  cell_cycle_genes <- read.csv(text = cc_file)
  
  cell_cycle_genes$geneSymbol <- ensemblToGeneName(cell_cycle_genes$geneID)
  
  # Acquire the S phase genes
  s_genes <- cell_cycle_genes$geneSymbol[cell_cycle_genes$phase == "S"]
  g2m_genes <- cell_cycle_genes$geneSymbol[cell_cycle_genes$phase == "G2/M"]
  
  cc_genes <- list("s_genes" = s_genes, "g2m_genes" = g2m_genes)
  cc_genes
}


# Plots an NxN grid of PC scatterplots. Quick way to plot multiple PC dimensions
# against each other.
plotPcaGrid <- function(scRNA, dims = 1:4, group = "orig.ident") {
  plots <- list()
  
  for (d1 in dims) {
    for (d2 in dims) {
      plt <- DimPlot(scRNA, dims=c(d1,d2), reduction = "pca", group.by = group, 
                     shuffle = TRUE, raster = TRUE) + 
        theme(legend.position = "none") + 
        labs(title = NULL)
      plots[[(d1-dims[1])*length(dims) + (d2-dims[1]+1)]] = plt
    }
  }
  
  plot_grid(plotlist = plots, nrow = length(dims), labels = c("A", "B", "C"))
}


# Uses Seurat's standard NormalizeData and ScaleData pipeline, except that we
# use 3000 features instead of 2000 to be consistent with the integration 
# pipeline.
simpleNormalizeAndScale <- function(scRNA, model = "negbinom") {
  scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize")
  scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 3000)
  scRNA <- ScaleData(scRNA, features = VariableFeatures(scRNA), 
                     model.use = model)
  scRNA
}


# Uses Seurat's standard NormalizeData and ScaleData pipeline but regresses out
# selected features from the metadata. Also uses 3000 features instead of 2000.
regressionNormalizeAndScale <- function(scRNA, model = "negbinom", 
                                        regress = c("percent.rb", "percent.mt", "nCount_RNA")) {
  scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize")
  scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 3000)
  scRNA <- ScaleData(scRNA, features = VariableFeatures(scRNA), 
                     model.use = model, 
                     vars.to.regress = regress)
  scRNA
}


# Uses Seurat's SCTransform pipeline to normalize and scale counts
SCTransformNormalize <- function (scRNA, regress = c("percent.mt", "percent.rb")) {
  scRNA <- SCTransform(scRNA, assay = "RNA", new.assay.name = "SCT",
                       vars.to.regress = regress)
  scRNA
}

# Uses Seurat's integration pipeline to integrate individual samples into one
# Seurat object. This is the method used for the paper, as we have two samples 
# that were run on a different day than the rest.
# Method can be "SCT" or "LogNormalize". Note: The LogNormalize code may not
# work, it has been awhile since I tried that method. SCT works.
generateIntegratedData <- function(scRNA, method = "SCT", 
                                   regress = c("percent.rb", "percent.mt", "nCount_RNA")) {
  split_seurat <- SplitObject(scRNA, split.by = "orig.ident")
  
  for (i in 1:length(split_seurat)) {
    if (method == "LogNormalize") {
      split_seurat[[i]] <- NormalizeData(split_seurat[[i]], 
                                         normalization.method = "LogNormalize")
      split_seurat[[i]] <- FindVariableFeatures(split_seurat[[i]], 
                                                selection.method = "vst", 
                                                nfeatures = 4000)
    }
    else { # SCT
      split_seurat[[i]] <- SCTransform(split_seurat[[i]], 
                                       vars.to.regress = regress)
    }
  }
  
  # Select the most variable features to use for integration
  integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                              nfeatures = 4000) 
  
  if (method == "SCT") {
    split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                       anchor.features = integ_features)
  }
  
  # Find anchors
  integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                          normalization.method = method,
                                          anchor.features = integ_features)
  
  # Integrate across conditions
  seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                     normalization.method = method)
  
  # Final scaling if using LogNormalize -- has to be linear, not negbinom
  if (method == "LogNormalize") {
    seurat_integrated <- ScaleData(seurat_integrated, 
                                   features = VariableFeatures(seurat_integrated),
                                   vars.to.regress = regress)
  }
  
  seurat_integrated
}
