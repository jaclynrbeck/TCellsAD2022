# Helper functions that are used in the analysis steps only (Steps 3-6)

# Author: Jaclyn Beck
# Final script used for paper as of Sep 02, 2022

library(readxl)
library(writexl)
library(biomaRt)

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


# Determines which cells are CD8+, CD4+, or gamma-delta single-positive and
# returns a list.
# Criteria:
# CD8+: Must express at least 1 count of either Cd8a or Cd8b1 AND doesn't 
#       express Cd4 or any gamma-delta receptor gene
# CD4+: Must express at least 1 count of Cd4 AND doesn't express Cd8a, Cd8b1, or
#       any gamma-delta receptor gene
# Gamma-delta: More stringent criteria are used to get "true" gamma-deltas. Must
#       express at least one count of a Trd gene or Sox13 or Blk, AND doesn't
#       express Cd4, Cd8a, or Cd8b1. 
# Returns a list containing fields for "CD8", "CD4", "GD" (gamma-delta), 
#       "Double" (double-positive for one of the above), "Negative" (cells that
#       didn't meet any of the above criteria)
getSinglePositiveCells <- function( scRNA ) {
  assay <- GetAssayData(scRNA, slot = "counts", assay = "RNA")
  all.genes <- rownames(assay)
  
  assay.cd4 <- assay["Cd4",]
  cd4.pos <- subset(names(assay.cd4), assay.cd4 >= 1)
  
  assay.cd8 <- assay[c("Cd8a", "Cd8b1"),]
  cd8.pos <- subset(colnames(assay.cd8), assay.cd8["Cd8a",] >= 1 | assay.cd8["Cd8b1",] >= 1)
  
  # Get all cells that express any gamma-delta marker except Trgv2 and Tcrg-C2,
  # which seem to be anomalously expressed in many cells
  assay.gd <- assay[grep("Trdc|Trdv|Tcrg-V|Tcrg-C1|Tcrg-C3|Tcrg-C4|Sox13|Blk", 
                         all.genes, value=TRUE),]
  assay.gd <- colSums(assay.gd)
  gd.pos <- subset(names(assay.gd), assay.gd >= 1)
  
  # Get cells most likely to be true gamma-deltas: cells must express a Trd 
  # gene OR Sox13 or Blk, which are specific to gamma delta 17's
  assay.gd <- assay[grep("Trdc|Trdv|Sox13|Blk", all.genes, value=TRUE),]
  assay.gd <- colSums(assay.gd)
  gd.delta <- subset(names(assay.gd), assay.gd >= 1)
  
  cd4cd8 <- intersect(cd4.pos, cd8.pos)
  cd4gd <- intersect(cd4.pos, gd.pos) # use gd.pos and not gd.delta here
  cd8gd <- intersect(cd8.pos, gd.pos) # use gd.pos and not gd.delta here
  
  double.pos <- unique(c(cd4cd8, cd4gd, cd8gd))
  
  # Remove double positives
  cd8.pos <- cd8.pos[!(cd8.pos %in% double.pos)]
  cd4.pos <- cd4.pos[!(cd4.pos %in% double.pos)]
  gd.delta <- gd.delta[!(gd.delta %in% double.pos)]
  
  neg <- setdiff(colnames(scRNA), c(cd8.pos, cd4.pos, gd.delta, double.pos))
  
  list("CD8" = cd8.pos, "CD4" = cd4.pos, "GD" = gd.delta, "Double" = double.pos, 
       "Negative" = neg)
}


matchTcrs <- function( scRNA, tcr.anno, tcr.epi = NULL ) {
  tcr.match <- subset(tcr.anno, Seurat.Barcode %in% colnames(scRNA))

  if (is.null(tcr.epi)) {
    return(tcr.match)
  }
  
  merged <- merge(tcr.match, tcr.epi, by = c("ClonotypeId"), 
                  all.x = TRUE, all.y = FALSE)

  rownames(merged) <- merged$Sample
  return(merged)
}


##### File I/O #####

# Writes output of FindAllMarkers to an Excel file. Automatically replaces
# gene names that have been altered with "--1" at the end, with the actual
# gene name. Adds the Ensembl ID of the gene as a column.
# sig.markers should be a data frame as output by FindAllMarkers.
# out.file should be a string with the full file path of the output file.
writeDifferentialGenes <- function( sig.markers, out.file ) {
  sheets <- c()
  
  sig.markers$Ensembl.ID <- geneNameToEnsembl(sig.markers$gene)
  sig.markers$Gene <- str_replace(sig.markers$gene, "--.*", "")
  
  for (i in levels(sig.markers$cluster)) {
    sub <- filter(sig.markers, cluster==i)
    sub <- sub[order(sub$avg_log2FC, decreasing=TRUE),]
    
    # cut out "p_val", "cluster", and old "gene" columns
    sub <- sub[, c("Ensembl.ID", "Gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
    
    # Some clusters have a "/" in their names (i.e. Naive/CM), which is an
    # illegal character for sheet names. Replace it with "+".
    sheets[[str_replace(i, "/", "+")]] <- sub
  }
  
  write_xlsx(sheets, path=out.file)
}


# Writes information about which clonotypes are in each cluster, and adds
# what percentage of cells in the cluster have each clonotype vs percent of
# cells outside the cluster. 
# scRNA: Seurat object
# tcr.anno: data frame read from file_clonotypes_processed
# clono.file: string with full file path to the CellRanger "clonotypes.csv" 
#             file. Use file_clonotypes_summary. 
# out.file: string with full file path to output file.
writeClusterClonotypes <- function( scRNA, tcr.anno, clono.file, out.file ) {
  sheets <- c()
  
  clono <- read.csv(clono.file)
  rownames(clono) <- clono$clonotype_id
  
  all_counts <- table(tcr.anno$ClonotypeId)
  
  for (i in levels(scRNA)) {
    sub <- subset(scRNA, idents=i)
    cells <- colnames(sub)
    
    clust_tcr <- filter(tcr.anno, Seurat.Barcode %in% cells)
    cl_counts <- as.data.frame(table(clust_tcr$ClonotypeId, clust_tcr$Sample), 
                               stringsAsFactors = FALSE)
    cl_counts <- subset(cl_counts, Freq > 0)
    cl_df <- data.frame(ClonotypeId = cl_counts$Var1, 
                        Clonotype.Sequence = clono[cl_counts$Var1, "cdr3s_aa"],
                        Count = cl_counts$Freq,
                        Sample = cl_counts$Var2,
                        iNKT = (clono[cl_counts$Var1, "inkt_evidence"] != ""),
                        MAIT = (clono[cl_counts$Var1, "mait_evidence"] != ""))
    cl_df$Clonotype.Sequence <- str_replace_all(cl_df$Clonotype.Sequence, ";", "; ")
    
    cl_df$Percent.Cells.In.Clust <- cl_df$Count / sum(cl_df$Count) * 100
    
    total_counts <- all_counts[cl_df$ClonotypeId]
    cl_df$Percent.Cells.Out.Of.Clust <- (total_counts - cl_df$Count) / sum(total_counts) * 100
    
    cl_df <- cl_df[order(cl_df$Percent.Cells.In.Clust, decreasing=TRUE),]
    sheets[[str_replace(i, "/", "+")]] <- cl_df
  }
  
  write_xlsx(sheets, path=out.file)
}


# Calculates differential genes between genotypes, running comparisons on 
# <genotype> vs all other genotypes, and <genotype> vs each individual genotype.
# Automatically replaces gene names that have been altered with "--1" at the 
# end, with the actual gene name. Adds the Ensembl ID of the gene as a column.
# scRNA: Seurat object
# genotypes: list or vector of genotype names, corresponding to what is in
#            scRNA$genotype. i.e. c("WT", "5XFAD", "PS19", "PS-5X")
# FDR: threshold for significance, i.e. 0.01
# out.file: string with full file path to output file
writeGenotypeDifferentialGenes <- function( scRNA, genotypes, FDR, out.file ) {
  sheets <- c()
  
  for (G1 in genotypes) {
    # G1 vs All
    markers <- FindMarkers(scRNA, ident.1 = G1, group.by = "genotype",
                           test.use = "MAST", min.pct = 0.05,
                           latent.vars = c("nCount_RNA", "StressScoreVariable"))

    sig.markers <- filter(markers, p_val_adj <= FDR & pct.1 >= 0.05)
    sig.markers$Ensembl.ID <- geneNameToEnsembl(rownames(sig.markers))
    sig.markers$Gene <- str_replace(rownames(sig.markers), "--.*", "")
    
    sub <- sig.markers[order(sig.markers$avg_log2FC, decreasing=TRUE),]
    
    # cut out "p_val"
    sub <- sub[, c("Ensembl.ID", "Gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]

    sheets[[paste(G1, "vs All")]] <- sub
    
    # G1 vs G2 for each G2 in <genotypes>
    for (G2 in genotypes[genotypes != G1]) {
      markers <- FindMarkers(scRNA, ident.1 = G2, ident.2 = G1, min.pct = 0.05,
                             group.by = "genotype", test.use = "MAST",
                             latent.vars = c("nCount_RNA", "StressScoreVariable"))
      
      sig.markers <- filter(markers, p_val_adj <= FDR & pct.1 >= 0.05)
      sig.markers$Ensembl.ID <- geneNameToEnsembl(rownames(sig.markers))
      sig.markers$Gene <- str_replace(rownames(sig.markers), "--.*", "")
      
      sub <- sig.markers[order(sig.markers$avg_log2FC, decreasing=TRUE),]
      
      # cut out "p_val"
      sub <- sub[, c("Ensembl.ID", "Gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]

      sheets[[paste(G2, "vs", G1)]] <- sub
    }
    
    write_xlsx(sheets, path=out.file)
  }
}


# Calculates differential genes between genotypes, on a per-cluster basis. 
# Only runs <genotype> vs All comparison. Automatically accounts for gene names
# with "--1" at the end. 
# scRNA: Seurat object
# genotypes: list or vector of genotype names, corresponding to what is in
#            scRNA$genotype. i.e. c("WT", "5XFAD", "PS19", "PS-5X")
# FDR: threshold for significance, i.e. 0.01
# out.file: string with full file path to output file
writeClusterVGenotypeDiffGenes <- function( scRNA, genotypes, FDR, out.file ) {
  sheets <- c()
  clust <- unique(Idents(scRNA))
  
  for (C in clust) {
    for (G in genotypes) {
      markers <- FindMarkers(scRNA, ident.1 = G, subset.ident = C, 
                             group.by = "genotype", 
                             test.use = "MAST",
                             latent.vars = c("nCount_RNA", "StressScoreVariable"))
      
      sig.markers <- filter(markers, p_val_adj <= FDR & pct.1 >= 0.1)
      sig.markers$Ensembl.ID <- geneNameToEnsembl(rownames(sig.markers))
      sig.markers$Gene <- str_replace(rownames(sig.markers), "--.*", "")
      
      sub <- sig.markers[order(sig.markers$avg_log2FC, decreasing=TRUE),]
      
      # cut out "p_val"
      sub <- sub[, c("Ensembl.ID", "Gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
      
      sheets[[paste(G, "vs All -", C)]] <- sub
    }
  }
  
  write_xlsx(sheets, path=out.file)
}


# Writes information about which clonotypes are in each genotype, and adds
# what percentage of cells in the genotype have each clonotype.
# scRNA: Seurat object
# tcr.anno: data frame read from file_clonotypes_processed
# genotypes: list or vector of genotype names, corresponding to what is in
#            scRNA$genotype. i.e. c("WT", "5XFAD", "PS19", "PS-5X")
# clono.file: string with full file path to the CellRanger "clonotypes.csv" 
#             file. Use file_clonotypes_summary. 
# out.file: string with full file path to output file.
writeGenotypeClonotypes <- function( scRNA, tcr.anno, genotypes, clono.file, out.file ) {
  sheets <- c()
  
  clono <- read.csv(clono.file)
  rownames(clono) <- clono$clonotype_id
  
  all_counts <- table(tcr.anno$ClonotypeId)
  
  for (i in genotypes) {
    sub <- scRNA$genotype[scRNA$genotype == i]
    cells <- names(sub)
    
    clust_tcr <- filter(tcr.anno, Seurat.Barcode %in% cells)
    cl_counts <- as.data.frame(table(clust_tcr$ClonotypeId, clust_tcr$Sample), 
                               stringsAsFactors = FALSE)
    cl_counts <- subset(cl_counts, Freq > 0)
    cl_df <- data.frame(ClonotypeId = cl_counts$Var1, 
                        Clonotype.Sequence = clono[cl_counts$Var1, "cdr3s_aa"],
                        Count = cl_counts$Freq,
                        Sample = cl_counts$Var2,
                        iNKT = (clono[cl_counts$Var1, "inkt_evidence"] != ""),
                        MAIT = (clono[cl_counts$Var1, "mait_evidence"] != ""))
    cl_df$Clonotype.Sequence <- str_replace_all(cl_df$Clonotype.Sequence, ";", "; ")
    
    cl_df$Percent.Cells.In.Geno <- cl_df$Count / sum(cl_df$Count) * 100
    cl_df <- cl_df[order(cl_df$Percent.Cells.In.Geno, decreasing=TRUE),]
    
    sheets[[str_replace(i, "/", "+")]] <- cl_df
  }
  
  write_xlsx(sheets, path=out.file)
}


# Reads significant genes from specific tabs of the genotype diff genes file.
# filename: full file path to the genotype diff genes file, generated from
#           writeGenotypeDifferentialGenes().
# pattern: regular expression for grep search on tab names. Example: "All"
#          will find all sheets with "<genotype> vs All". "vs WT" would find
#          all sheets with "<genotype> vs WT". 
readSigGenesGenotype <- function ( filename, pattern = "All" ) {
  diff.genes <- lapply(excel_sheets(filename), read_excel, 
                       path = filename)
  names(diff.genes) <- excel_sheets(filename)
  
  sig.genes <- Map(function(i, x) {
    x$cluster <- i
    x
  }, names(diff.genes), diff.genes)
  
  sig.genes.df <- do.call(rbind, sig.genes)
  
  comparisons <- unique(grep(pattern, sig.genes.df$cluster, value = TRUE))
  sig.genes.df <- subset(sig.genes.df, cluster %in% comparisons)
  
  sig.genes.df
}


##### Things that print to the console #####

# Prints the top 10 markers for each cluster
# sig.markers: data frame output by FindAllMarkers
# pos.only: TRUE = display only positively-changed genes. 
#           FALSE = display both positively and negatively changed genes
printTop10Markers <- function ( sig.markers, pos.only = TRUE ) {
  if (pos.only) {
    sorted <- sig.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  }
  else {
    sorted <- sig.markers %>% group_by(cluster) %>% top_n(n = 10, wt = abs(avg_log2FC))
  }
  visual_sorted <- data.frame(Cluster=unique(sorted$cluster),
                              Genes=" ")
  genes <- sapply(visual_sorted$Cluster, 
                  function(x,y) 
                    list(subset(y, cluster==x)$gene),
                  sorted)
  visual_sorted$Genes <- genes
  
  tmp <- data.frame(visual_sorted$Genes[1])
  for (r in 2:nrow(visual_sorted)) {
    tmp <- cbind(tmp, visual_sorted$Genes[r])
  }
  colnames(tmp) <- paste("Cluster", 1:nrow(visual_sorted)-1)
  print(tmp)
  tmp
}


# Prints the percentage of each genotype in each cluster, and the percentage
# of each cluster in each genotype. 
printClusterDistributions <- function(scRNA) {
  summ_ident_table = table(scRNA$orig.ident, Idents(scRNA))
  
  # Distribution across clusters of each group separately (sum of columns = 100)
  clust_per_group = round(t(summ_ident_table / rowSums(summ_ident_table))*100, 2)
  print("Dist. of each cluster within each genotype (columns sum to 100)")
  print(clust_per_group)
  
  # Distribution of groups in each cluster (sum of rows = 100)
  print("")
  print("Dist. of each genotype within each cluster (rows sum to 100)")
  print(round((clust_per_group / rowSums(clust_per_group))*100, 2))
}


# Gets the 10 highest expressed genes (by scaled count) in each cluster.
# Uses the "integrated" assay. 
getHighestExpressedGenes <- function( scRNA, ngenes = 10 ) {
  DefaultAssay(scRNA) <- "integrated"
  highest.expressed = data.frame()
  for (i in levels(scRNA)) {
    sub <- subset(scRNA, idents=i)
    data <- GetAssayData(object = sub, slot = "scale.data")
    top.genes <- sort(rowMeans(data), decreasing=TRUE)[1:ngenes]
    if (length(highest.expressed) == 0) {
      highest.expressed <- data.frame("0"=names(top.genes))
    }
    else {
      highest.expressed <- cbind(highest.expressed, names(top.genes))
    }
  }
  
  colnames(highest.expressed) <- paste("Cluster", levels(scRNA), sep=" ")
  DefaultAssay(scRNA) <- "RNA"
  highest.expressed
}





##### Not refactored yet #####


# Queries BioMart to get mouse homologs for human genes
# human.genes: vector of gene names
getHumanMouseHomologs <- function(human.genes) {
  # Temporarily using archive because new ensembl update crashes biomart
  ens.human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", 
                       host = "https://dec2021.archive.ensembl.org/")
  ens.mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", 
                       host = "https://dec2021.archive.ensembl.org/")
  
  homologs <- getLDS(attributes = c('hgnc_symbol'),
                     filters = 'hgnc_symbol', 
                     values = unique(human.genes), 
                     mart = ens.human,
                     attributesL = c('mgi_symbol'), 
                     martL = ens.mouse)
}


plotGeneExpressionUMAP <- function(scRNA, gene.list, zero.cutoff = TRUE) {
  for (M in seq(1,length(gene.list), 8)) {
    end = min(7, length(gene.list) - M)
    
    if (zero.cutoff == TRUE) {
      plt = FeaturePlot(scRNA, features = gene.list[M:(M+end)], pt.size = 0.5,
                        ncol = 4, label = FALSE, cols = c("#eeeeee", "red"), 
                        order = TRUE, min.cutoff = 0)
    }
    else {
      plt = FeaturePlot(scRNA, features = gene.list[M:(M+end)], pt.size = 0.5,
                        ncol = 4, label = FALSE, #cols = c("blue", "red"), 
                        order = TRUE) & 
            scale_color_viridis(option = "plasma")
    }
    
    print(plt)
  }
}


plotGeneExpressionViolin <- function(scRNA, gene.list, zero.cutoff = TRUE, 
                                     group.by = "orig.ident", colors = NULL) {
  for (M in seq(1,length(gene.list), 4)) {
    end = min(3, length(gene.list) - M)
    
    if (zero.cutoff == TRUE) {
      plt1 <- FeaturePlot(scRNA, features = gene.list[M:(M+end)], pt.size = 0.5,
                         ncol = 4, label = FALSE, cols = c("#eeeeee", "red"), 
                         order = TRUE, min.cutoff = 0, combine = FALSE)
    }
    else {
      plt1 <- FeaturePlot(scRNA, features = gene.list[M:(M+end)], pt.size = 0.5,
                         ncol = 4, label = FALSE, order = TRUE, combine = FALSE)
      plt1 <- lapply(plt1, function(P) P + 
                       scale_color_viridis(option = "plasma"))
    }
    
    plt2 <- VlnPlot(scRNA, features = gene.list[M:(M+end)], pt.size = 0, 
                    ncol = 4, cols = colors, group.by = group.by, 
                    combine = FALSE)
    plt2 <- lapply(plt2, function(P) P + NoLegend())
    
    
    print(Reduce('|', plt1) / Reduce('|', plt2))
  }
}


##### GO analysis functions #####

# Runs GProfiler to get statisticaly significant GO terms.
# markers: data frame of significant markers as output by FindMarkers
# sources: list of sources compatible with the "sources" argument of gost.
#          See "?gost" for possible values. 
runGOAnalysis <- function( markers, sources = c("GO", "REAC") ) {
  m.up <- subset(markers, avg_log2FC > 0)
  m.down <- subset(markers, avg_log2FC < 0)
  
  # Order by -log10(p value) : most significant first
  m.up <- m.up[order(-log10(m.up$p_val_adj), decreasing = TRUE),]
  m.down <- m.down[order(-log10(m.down$p_val_adj), decreasing = TRUE),]
  
  go.res.up <- list()
  go.res.down <- list()
  
  if (nrow(m.up) > 0) {
    go.res.up <- gost(query = m.up$Ensembl.ID, 
                      organism = "mmusculus", ordered_query = TRUE, 
                      multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                      measure_underrepresentation = FALSE, evcodes = TRUE, 
                      user_threshold = 0.05, correction_method = "g_SCS", 
                      domain_scope = "annotated", custom_bg = NULL, 
                      numeric_ns = "", sources = sources, 
                      as_short_link = FALSE)
  }
  
  if (nrow(m.down) > 0) {
    go.res.down <- gost(query = m.down$Ensembl.ID, 
                        organism = "mmusculus", ordered_query = TRUE, 
                        multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                        measure_underrepresentation = FALSE, evcodes = TRUE, 
                        user_threshold = 0.05, correction_method = "g_SCS", 
                        domain_scope = "annotated", custom_bg = NULL, 
                        numeric_ns = "", sources = sources, 
                        as_short_link = FALSE)
  }
  
  if (!is.null(go.res.up)) {
    go.res.up$result$Phenotype <- "+1"
  }
  if (!is.null(go.res.down)) {
    go.res.down$result$Phenotype <- "-1"
  }
  
  return(list("Up" = go.res.up, "Down" = go.res.down))
}


# Converts gost output to a Gem file for use with Cytoscape
# go.res.list: a list with data frames output by gost. Items in list should
#              be named "Up" and "Down"
# outfile: string with the full file path to the output gem file
# use.ensembl: TRUE = use Ensembl IDs in the gem file. (recommended)
#              FALSE = use gene names in the gem file.
goResultsToGem <- function( go.res.list, outfile, use.ensembl = TRUE ) {
  go.res.up <- go.res.list[["Up"]]
  go.res.down <- go.res.list[["Down"]]
  
  go.res <- NULL
  
  if (!is.null(go.res.up) & !is.null(go.res.down)) {
    go.res <- rbind(go.res.up$result, go.res.down$result)
  }
  else if (!is.null(go.res.up)) {
    go.res <- go.res.up
  }
  else if (!is.null(go.res.down)) {
    go.res <- go.res.down
  }
  
  if (!is.null(go.res)) {
    go.res <- subset(go.res, term_size <= 1000)
    gem <- go.res[,c("term_id", "term_name", "p_value", "intersection", "Phenotype")]  
    colnames(gem) <- c("GO.ID", "Description", "p.Val", "Genes", "Phenotype")
    gem$FDR <- gem$p.Val
    gem <- gem[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Genes")]
    
    if (!use.ensembl) {
      tmp <- sapply(str_split(gem$Genes, pattern = ","), ensemblToGeneName)
      gem$Genes <- sapply(tmp, str_c, collapse = ",")
    }
    
    write.table(gem, file = outfile, sep = "\t", quote = F, row.names = F)
  }
}


