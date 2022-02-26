# Helper functions that are used more than once in the code

matchTcrs <- function( scRNA, tcr.anno, tcr.epi = NULL ) {
  matches <- intersect(colnames(scRNA), tcr.anno$Sample)
  tcr.match <- subset(tcr.anno, Sample %in% matches)

  if (is.null(tcr.epi)) {
    return(tcr.match)
  }
  
  merged <- merge(tcr.match, tcr.epi, by = c("Clonotype", "ClonotypeId"), 
                  all.x = TRUE, all.y = FALSE)

  rownames(merged) <- merged$Sample
  return(merged)
}


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

getHighestExpressedGenes <- function( scRNA ) {
  highest.expressed = data.frame()
  for (i in levels(scRNA)) {
    sub <- subset(scRNA, idents=i)
    data <- GetAssayData(object = sub, slot = "data")
    top.genes <- sort(rowMeans(data), decreasing=TRUE)[1:10]
    if (length(highest.expressed) == 0) {
      highest.expressed <- data.frame("0"=names(top.genes))
    }
    else {
      highest.expressed <- cbind(highest.expressed, names(top.genes))
    }
  }
  
  colnames(highest.expressed) <- paste("Cluster", levels(scRNA), sep=" ")
  highest.expressed
}

geneNameToEnsembl <- function(genes) {
  features <- read.table(file = gzfile(file.path(dir_filtered_counts, "features.tsv.gz")))
  rownames(features) <- features$V2
  features[genes,"V1"]
}

writeDifferentialGenes <- function( sig.markers, out.file ) {
  sheets <- c()
  
  sig.markers$Ensembl.ID <- geneNameToEnsembl(sig.markers$gene)
  sig.markers$Gene <- str_replace(sig.markers$gene, "--.*", "")
  
  for (i in levels(sig.markers$cluster)) {
    sub <- filter(sig.markers, cluster==i)
    sub <- sub[order(sub$avg_log2FC, decreasing=TRUE),]
    
    # cut out "p_val", "cluster", and old "gene" columns
    sub <- sub[, c("Ensembl.ID", "Gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
    sheets[[str_replace(i, "/", "+")]] <- sub
  }
  
  write_xlsx(sheets, path=out.file)
}

writeClusterClonotypes <- function( scRNA, tcr.anno, clono.file, out.file ) {
  sheets <- c()
  
  clono <- read.csv(clono.file)
  rownames(clono) <- clono$clonotype_id
  
  all_counts <- table(tcr.anno$ClonotypeId)
  
  for (i in levels(scRNA)) {
    sub <- subset(scRNA, idents=i)
    cells <- colnames(sub)
    
    clust_tcr <- filter(tcr.anno, Sample %in% cells)
    cl_counts <- as.data.frame(table(clust_tcr$ClonotypeId, clust_tcr$Genotype), 
                               stringsAsFactors = FALSE)
    cl_counts <- subset(cl_counts, Freq > 0)
    cl_df <- data.frame(ClonotypeId = cl_counts$Var1, 
                        Clonotype = clono[cl_counts$Var1, "cdr3s_aa"],
                        Count = cl_counts$Freq,
                        Sample = cl_counts$Var2,
                        iNKT = (clono[cl_counts$Var1, "inkt_evidence"] != ""),
                        MAIT = (clono[cl_counts$Var1, "mait_evidence"] != ""))
    cl_df$Clonotype <- str_replace_all(cl_df$Clonotype, ";", "; ")
    
    cl_df$Percent.Cells.In.Clust <- cl_df$Count / sum(cl_df$Count) * 100
    
    total_counts <- all_counts[cl_df$ClonotypeId]
    cl_df$Percent.Cells.Out.Of.Clust <- (total_counts - cl_df$Count) / sum(total_counts) * 100
    
    cl_df <- cl_df[order(cl_df$Percent.Cells.In.Clust, decreasing=TRUE),]
    sheets[[str_replace(i, "/", "+")]] <- cl_df
  }
  
  write_xlsx(sheets, path=out.file)
}

writeGenotypeDifferentialGenes <- function( scRNA, genotypes, FDR, out.file ) {
  sheets <- c()
  
  for (G1 in genotypes) {
    markers <- FindMarkers(scRNA, ident.1 = G1, group.by = "genotype")
    sig.markers <- filter(markers, p_val_adj <= FDR & pct.1 >= 0.1)
    sig.markers$Ensembl.ID <- geneNameToEnsembl(rownames(sig.markers))
    sig.markers$Gene <- str_replace(rownames(sig.markers), "--.*", "")
    
    sub <- sig.markers[order(sig.markers$avg_log2FC, decreasing=TRUE),]
    
    # cut out "p_val"
    sub <- sub[, c("Ensembl.ID", "Gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]

    sheets[[paste(G1, "vs All")]] <- sub
    
    for (G2 in genotypes[genotypes != G1]) {
      markers <- FindMarkers(scRNA, ident.1 = G2, ident.2 = G1, group.by = "genotype")
      sig.markers <- filter(markers, p_val_adj <= FDR & pct.1 >= 0.1)
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

writeClusterVGenotypeDiffGenes <- function( scRNA, genotypes, FDR, out.file ) {
  sheets <- c()
  clust <- unique(Idents(scRNA))
  
  for (C in clust) {
    for (G in genotypes) {
      markers <- FindMarkers(scRNA, ident.1 = G, subset.ident = C, 
                             group.by = "genotype")
      
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


writeGenotypeClonotypes <- function( scRNA, tcr.anno, genotypes, clono.file, out.file ) {
  sheets <- c()
  
  clono <- read.csv(clono.file)
  rownames(clono) <- clono$clonotype_id
  
  all_counts <- table(tcr.anno$ClonotypeId)
  
  for (i in genotypes) {
    sub <- scRNA$genotype[scRNA$genotype == i]
    cells <- names(sub)
    
    clust_tcr <- filter(tcr.anno, Sample %in% cells)
    cl_counts <- as.data.frame(table(clust_tcr$ClonotypeId, clust_tcr$Genotype), 
                               stringsAsFactors = FALSE)
    cl_counts <- subset(cl_counts, Freq > 0)
    cl_df <- data.frame(ClonotypeId = cl_counts$Var1, 
                        Clonotype = clono[cl_counts$Var1, "cdr3s_aa"],
                        Count = cl_counts$Freq,
                        Sample = cl_counts$Var2,
                        iNKT = (clono[cl_counts$Var1, "inkt_evidence"] != ""),
                        MAIT = (clono[cl_counts$Var1, "mait_evidence"] != ""))
    cl_df$Clonotype <- str_replace_all(cl_df$Clonotype, ";", "; ")
    
    cl_df$Percent.Cells.In.Geno <- cl_df$Count / sum(cl_df$Count) * 100
    
    total_counts <- all_counts[cl_df$ClonotypeId]
    
    cl_df <- cl_df[order(cl_df$Percent.Cells.In.Geno, decreasing=TRUE),]
    
    sheets[[str_replace(i, "/", "+")]] <- cl_df
  }
  
  write_xlsx(sheets, path=out.file)
}


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

getCellTypeMarkers <- function(all.genes) {
  markers.pos <- list("CD8" = c("Cd8[a|b]"),
                     
                     "CD4" = c("^Cd4$"),
                     
                     "CD3" = c("^Cd3[d|e|g]$"),
                     
                     "Naive/CM" = c("Sell", "Ccr7", "Ptprc$", "Il7r", "Tcf7$", 
                                    "Cd27$", "Cd28", "Klf2", "S1pr1"),
                     #, "Lef1$", "Satb1", "Gpr183", "Ltb$", "S100a10"),
                     
                     "CM" = c("Sell", "Ccr7", "Cd44", "Cd27$", "Il7r", "Il2rb",
                              "Eomes", "Bcl6", "Tcf1$", "Stat3", "^Id3$", "^Fas$"),
                     
                     "Effector/EM" = c("Klrg1", "Cx3cr1", "S1pr5", "Tnf$", 
                                       "Ifng$", "Gzmb", "Prf1", "Cd44", "^Fas$"),
                     #"Stat4", "Prdm1", "Eomes
                     #"Id2", "Tbx21", "Bcl6"),
                     
                     "RM" = c("Itgae", "Itga1$", "Cxcr6", "Cd69", "Cxcr3", "Il7r",
                              "Zfp683", "Prdm1$"),
                     
                     "Exhausted" = c("Pdcd1$", "Havcr2", "Lag3", "Tigit", "Cd160", 
                                     "Ctla4"), #"Cxcl13", "Hspb1$", "Cxcr6", "Irf4",
                     #"Layn", "Gimap6", "Hsph1"),
                     
                     "Interferon" = c("Bst2", "Irf1", "Irf2$", "Irf7", "Stat1", 
                                      "Stat2", "^Mx1"), # Mx1 not present in data set
                     
                     "Proliferating" = c("Mki67"), #, "Kif"), #"Cdk[0-9]*$", 
                     #"^Cdc[0-9]*$"),
                     
                     "Activation" = c("Klrg1", "Il2ra", "Cd69", "Tnfrsf4", 
                                      "Tnfrsf9", "Cd44", "Cd40lg", "^Fas$", 
                                      "Ccl5", "Klrk1", "Cx3cr1", "Icos$"),
                     
                     "Cytotoxic" = c("Prf1", "Gzm", "Lamp1", "Ifng$", "Tnf$", 
                                     "Tnfsf10", "Nkg7"),
                     
                     "Anergic/Tolerant" = c("Lag3", "Pdcd1$", "Nfatc[12]$", 
                                            "Nfkb1", "Ikzf1", "Egr[12]"),
                     
                     "GammaDelta" = c("Trg", "Trdv", "Trdc", "Il17a", "Rorc", 
                                      "Cxcr5", "Il2rb", "Tcrg", "Sox13", "Blk",
                                      "Cd163l1", "5830411N06Rik"),
                     
                     #"AlphaBeta" = c("Trac", "Traj", "Trav", "^Trb"),
                     
                     "MAIT" = c("Trav1$", "Traj33", "Trbv19", "Trbv13", #"Trav19", "Trbv[68]"
                                "Klrb1c", "Zbtb16", "Rorc",
                                "Tbx21", "Ccr5", "Ccr6"), #"Slc4a10", 
                     
                     "NKT" = c("Trav11$", "Traj18", "Trbv13-2", "Trbv1$", 
                               "Trbv29", "Klrb1c", # Trav14-[123]", "Trbv[278]$"
                               "Rorc", "Zbtb16", "Gata3", "Tbx21", "Ifng$", 
                               "^Il4$", "Nkg7"),
                     #"Klrc2", "Zfp683", "^Xcl1$"),
                     
                     "Treg" = c("Il2ra", "Foxp3", "Ctla4", "Il7r", "^Il2$", "Tigit", "Nrp1"), 
                     #"Ccr10"), # Tregs can express Gata3, Rorc, and Tbx21 too
                     
                     "Th1" = c("Tbx21", "Ifng$", "Cxcr6", "Cxcr3", "Ccr5", "Il12rb2"),
                     
                     "Th2" = c("Gata3", "^Il4$", "Ccr4", "^Il5", "Il13"), # Il5 not present in data set
                     
                     "Th9" = c("Spi1", "Il9$", "Ccr3", "Ccr6"), # Il9 not present in data set
                     
                     "Th17" = c("Rorc", "Il17[a|f]", "Il21$", "Il22", "Il23r", "Ccr5", "Il10$", "Ccr6", "Ccr4", "Klrb1c"),
                     
                     "Th22" = c("Ahr", "Il22", "Ccr6", "Ccr10", "Ccr4"),
                     
                     "Tfh" = c("Bcl6$", "Cxcr5", "Tcf7$", "Il21$", "Icos$", "Cd40lg"),
                     
                     "Tr1" = c("Il10$"), # Subset of Tregs
                     
                     "Tfr" = c("Cxcr5", "Ctla4") # Subset of Tregs
  )

  for (item in names(markers.pos)) {
    markers.pos[item] <- list(c(unname(unlist(sapply(markers.pos[item][[1]], function(x,y) grep(x, y, ignore.case=TRUE, value=TRUE), all.genes)))))
  }
  
  markers.pos
}

convertMarkerListToDf <- function(markers) {
  markers <- stack(markers)
  colnames(markers) <- c("Gene", "Type")
  markers
}

# Fig 2 of Magen paper
getInterestingMarkers <- function() {
  markers <- list("Tissue Residence" = c("S1pr1", "Klf2", "Cd69", "Cxcr6"),
                  "Chemokine" = c("Cxcr3", "Cxcr5", "Cxcr4", "Ccr7", "Ccr2", 
                                  "Ccr5", "Ccl3", "Ccl4", "Ccl5", "Cxcl10", 
                                  "Ccrl2", "Ccl1", "Ccr8"),
                  "Cytokine" = c("^Il4$", "Il9r", "Il21$", "^Il2$", "Il7r", 
                                 "Il27ra", "Il6st", "Il1rl1", "Il2ra", "Il10$", 
                                 "Ifng$", "Il12rb2", "Il10ra", "Il18rap", 
                                 "Il2rb", "Il2rg", "Il18r1"),
                  "Costimulatory" = c("Tnfsf11", "Tnfsf8", "^Icos$", "Tnfrsf4",
                                      "Tnfrsf9", "Tnfrsf18"),
                  "Tx Factors" = c("^Id3$", "Lef1$", "Tcf7$", "Klf4", "Bcl6", 
                                   "Nr4a3", "Nr4a2", "Bhlhe40", "Gata3", 
                                   "Pdcd11", "Prdm1$", "Foxp3", "Ikzf2", 
                                   "^Maf$", "Tbx21", "Klf3", "Nr4a1", "Rora", 
                                   "^Sub1", "Stat4", "Irf9", "Irf7", "Stat1", 
                                   "Klf6", "Stat2", "Runx3", "^Id2$", "Runx2$", 
                                   "Irf8", "Stat3")
  )
  
  for (item in names(markers)) {
    markers[item] <- list(c(unname(unlist(sapply(markers[item][[1]], function(x,y) grep(x, y, ignore.case=TRUE, value=TRUE), all.genes)))))
  }
  
  markers
}

printCellTypeAssignments <- function(sig.markers.pos, sig.markers.neg) {
  for (clust in unique(sig.markers$cluster)) {
    genes.pos <- filter(sig.filtered.pos, cluster == clust)$gene
    genes.neg <- filter(sig.filtered.neg, cluster == clust)$gene
    
    print(paste("Cluster", clust, "POSITIVE:"))
    print(genes.pos)
    pos <- unstack(filter(markers.pos, Gene %in% genes.pos), form = Type ~ Gene)
    if (class(pos) == "list") {
      pos.df <- data.frame(Gene = names(pos), Types = sapply(pos, function(x) str_flatten(x, collapse=", ")))
      print(pos.df)
    }
    else {
      print(pos)
    }
    
    print("")
    
    print(paste("Cluster", clust, "NEGATIVE:"))
    print(genes.neg)
    neg <- unstack(filter(markers.pos, Gene %in% genes.neg), form = Type ~ Gene)
    if (class(neg) == "list") {
      neg.df <- data.frame(Gene = names(neg), Types = sapply(neg, function(x) str_flatten(x, collapse=", ")))
      print(neg.df)
    }
    else {
      print(neg)
    }
    
    print("")
  }
}

getGammaDeltas2 <- function( scRNA ) {
  DefaultAssay(scRNA) <- "RNA"
  scRNA <- DietSeurat(scRNA, assays = c("RNA"))
  
  all.genes <- rownames(scRNA)
  tras <- grep("^Tra[v|j]", all.genes, value=TRUE)
  tratrd <- grep("-dv", tras, value=TRUE) # Could be alpha or delta
  tras <- setdiff(tras, tratrd)
  
  receptors.df <- list("TRA" = tras, 
                       "TRA.C" = grep("^Trac", all.genes, value=TRUE),
                       "TRB" = grep("^Trbv", all.genes, value=TRUE),
                       "TRB.C" = grep("^Trbc", all.genes, value=TRUE),
                       #"TRA.TRD" = tratrd, 
                       "TRD" = grep("^Trdv", all.genes, value=TRUE),
                       "TRD.C" = grep("^Trdc", all.genes, value=TRUE),
                       "TRG" = c(grep("^Trgv", all.genes, value=TRUE), 
                                 grep("^Tcrg-V", all.genes, value=TRUE)),
                       "TRG.C" = grep("^Tcrg-C", all.genes, value=TRUE))
                       #"CD" = c("Cd4", "Cd8a", "Cd8b1"))
  
  receptors.df <- stack(receptors.df) %>% as.data.frame()
  colnames(receptors.df) <- c("Gene", "Category")
  rownames(receptors.df) <- receptors.df$Gene
  
  categories <- unique(receptors.df$Category)
  assay <- GetAssayData(scRNA, slot="counts")
  assay <- assay[receptors.df$Gene,]
  
  neg.cells <- colnames(scRNA)[which(!(colnames(scRNA) %in% tcr.anno$Sample))]
  tcr.ab <- subset(receptors.df, Category %in% c("TRA", "TRA.C", "TRB", "TRB.C"))
  tcr.gd <- subset(receptors.df, Category %in% c("TRD", "TRD.C", "TRG", "TRG.C"))
  assay.gd <- assay[tcr.gd$Gene, neg.cells]
  gd <- colnames(assay.gd)[which(colSums(assay.gd) > 1)]
  
  colmax.ab <- apply(assay[tcr.ab$Gene, neg.cells], 2, max)
  colmax.gd <- apply(assay[tcr.gd$Gene, neg.cells], 2, max)
  
  gd2 <- neg.cells[which(colmax.gd > colmax.ab)]
  
  scRNA2 <- scRNA[,neg.cells]
  
  for (cat in categories) {
    genes <- subset(receptors.df, Category == cat)
    if (length(genes$Gene) == 1) {
      colmax <- assay[genes$Gene,]
    }
    else {
      colmax <- apply(assay[genes$Gene,], 2, max)
    }
    scRNA <- AddMetaData(scRNA, colmax, col.name = paste0(cat, ".Max"))
    #scRNA2 <- AddModuleScore_UCell(scRNA2, features = genes$Gene, name = cat, 
    #                               maxRank = 10000, storeRanks = TRUE)
  }
  
  for (cat in categories) {
    #cols <- grep(paste0("^", cat, "[0-9]+"), colnames(scRNA2@meta.data), value=TRUE)
    cols <- grep(paste0(cat, "$"), colnames(scRNA2@meta.data), value=TRUE)
    meta <- scRNA2@meta.data[,cols]
    if (length(cols) == 1) {
      mx <- meta
    }
    else {
      mx <- apply(meta, 1, max)
    }
    scRNA2 <- AddMetaData(scRNA2, mx, col.name = paste0(cat, ".Max"))
  }
  
  abs <- paste0(c("TRA", "TRA.C", "TRB", "TRB.C"), ".Max")
  gds <- paste0(c("TRD", "TRD.C", "TRG", "TRG.C"), ".Max")
  
  scRNA2 <- AddMetaData(scRNA2, rowSums(scRNA2@meta.data[abs]), col.name = "AB.Score")
  scRNA2 <- AddMetaData(scRNA2, rowSums(scRNA2@meta.data[gds]), col.name = "GD.Score")
  
  scRNA2 <- AddMetaData(scRNA2, scRNA2$GD.Score > scRNA2$AB.Score, col.name = "Is.GD")
  
  gammas <- names(which(scRNA2$Is.GD))
}

getGammaDeltas <- function( scRNA, tcr.anno ) {
  all.genes <- rownames(scRNA)
  tras <- grep("^Tra[v|j]", all.genes, value=TRUE)
  tratrd <- grep("-dv", tras, value=TRUE) # Could be alpha or delta
  tras <- setdiff(tras, tratrd)
  
  receptors.df <- list("TRA" = tras, 
                       "TRA.C" = grep("^Trac", all.genes, value=TRUE),
                       "TRB" = grep("^Trbv", all.genes, value=TRUE),
                       "TRB.C" = grep("^Trbc", all.genes, value=TRUE),
                       #"TRA.TRD" = tratrd, 
                       "TRD" = grep("^Trdv", all.genes, value=TRUE),
                       "TRD.C" = grep("^Trdc", all.genes, value=TRUE),
                       "TRG" = c(grep("^Trgv", all.genes, value=TRUE), 
                                 grep("^Tcrg-V", all.genes, value=TRUE)),
                       "TRG.C" = grep("^Tcrg-C", all.genes, value=TRUE))
                       #"CD" = c("Cd4", "Cd8a", "Cd8b1"))
  
  receptors.df <- stack(receptors.df) %>% as.data.frame()
  colnames(receptors.df) <- c("Gene", "Category")
  rownames(receptors.df) <- receptors.df$Gene
  
  categories <- unique(receptors.df$Category)
  assay <- GetAssayData(scRNA, slot = "counts", assay = "RNA")
  assay <- assay[receptors.df$Gene,]
  
  neg.cells <- colnames(scRNA)[which(!(colnames(scRNA) %in% tcr.anno$Sample))]
  tcr.ab <- subset(receptors.df, Category %in% c("TRA", "TRA.C", "TRB", "TRB.C"))
  tcr.gd <- subset(receptors.df, Category %in% c("TRD", "TRD.C", "TRG", "TRG.C"))
  assay.gd <- assay[tcr.gd$Gene, neg.cells]
  gd <- colnames(assay.gd)[which(colSums(assay.gd) > 1)]
}


getSinglePositiveCells <- function( scRNA, tcr.anno ) {
  assay <- GetAssayData(scRNA, slot = "counts", assay = "RNA")
  all.genes <- rownames(assay)
  
  assay.cd4 <- assay["Cd4",]
  cd4.pos <- subset(names(assay.cd4), assay.cd4 >= 1)
  
  assay.cd8 <- assay[c("Cd8a", "Cd8b1"),]
  cd8.pos <- subset(colnames(assay.cd8), assay.cd8["Cd8a",] >= 1 | assay.cd8["Cd8b1",] >= 1)
  
  #assay.gd <- assay[grep("Trdc|Trdv4|Tcrg-V|Tcrg-C1|Sox13|Blk", all.genes, 
  #                       value=TRUE),]
  #assay.gd <- colSums(assay.gd)
  gd.pos <- getGammaDeltas(scRNA, tcr.anno) #subset(names(assay.gd), assay.gd >= 1)
  
  #assay.ab <- assay[grep("Trac|Traj|Trav|^Trb", all.genes, value=TRUE),]
  #assay.ab <- colSums(assay.ab)
  #ab.pos <- subset(names(assay.ab), assay.ab >= 1)
  
  cd4cd8 <- intersect(cd4.pos, cd8.pos)
  cd4gd <- c() #intersect(cd4.pos, gd.pos)
  cd8gd <- c() #intersect(cd8.pos, gd.pos)
  #abgd <- intersect(ab.pos, gd.pos)
  
  double.pos <- unique(c(cd4cd8, cd4gd, cd8gd))
  
  # Remove double positives
  cd8.pos <- cd8.pos[!(cd8.pos %in% double.pos)]
  cd4.pos <- cd4.pos[!(cd4.pos %in% double.pos)]

  neg <- setdiff(colnames(scRNA), c(cd8.pos, cd4.pos, gd.pos, double.pos))
  
  list("CD8" = cd8.pos, "CD4" = cd4.pos, "GD" = gd.pos, "Double" = double.pos, 
       "Negative" = neg)
}


getHumanMouseHomologs <- function(human.genes) {
  ens.human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  ens.mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
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


readSignatureGenes <- function(sig.gene.file) {
  signature.genes <- read.table(sig.gene.file, header=FALSE, sep="\t")
  colnames(signature.genes) <- c("Cluster", "H.gene")
  ms.genes <- getHumanMouseHomologs(signature.genes$H.gene)
  
  clusts <- unique(signature.genes$Cluster)
  ms.sig.genes <- list()
  
  for (id in clusts) {
    hu.genes <- subset(signature.genes, Cluster == id)
    for (gene in hu.genes$H.gene) {
      ms <- subset(ms.genes, HGNC.symbol == gene)$MGI.symbol
      ms.sig.genes[[id]] <- append(ms.sig.genes[[id]], ms)
    }
  }
  
  ms.sig.genes
}


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


goResultsToGem <- function( go.res.list, outfile ) {
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
    
    write.table(gem, file = outfile, sep = "\t", quote = F, row.names = F)
  }
}





