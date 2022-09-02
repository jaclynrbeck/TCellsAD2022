# This script was written to be run line-by-line because of exploratory plots
# to determine parameters. The values as entered in the script are those used
# for the paper.

# Author: Jaclyn Beck
# Final script used for paper as of April 14, 2022

library(Seurat)
library(dplyr)
library(stringr)
library(writexl)
library(ggplot2)
library(fgsea)
library(msigdbr)
library(gprofiler2)
source(file.path("functions", "Analysis_HelperFunctions.R"))
source(file.path("functions", "Analysis_PlottingFunctions.R"))
source(file.path("functions", "SeuratFromData_HelperFunctions.R"))
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

clusters <- list("0" = "CD8 Cytotoxic",
                 "1" = "Naive",
                 "2" = "CD4 Helper",
                 "3" = "CD8 CM + \u03b3\u03b4-1",
                 "4" = "\u03b3\u03b4-17", 
                 "5" = "CD8 Activated",
                 "6" = "Proliferating")
scRNA <- RenameIdents(scRNA, clusters)
scRNA$clusters <- Idents(scRNA)

saveRDS(scRNA, file_seurat_analyzed_allcells)

##### Find cluster markers #####

#scRNA <- readRDS(file_seurat_analyzed_allcells)
DefaultAssay(scRNA) <- "RNA"
scRNA <- NormalizeData(scRNA)

all.markers <- FindAllMarkers(scRNA, test.use = "MAST", 
                              logfc.threshold = 0.25, min.pct = 0.1,
                              latent.vars = c("nCount_RNA", "StressScoreVariable"))

saveRDS(all.markers, file = file_markers_combined_all)

#all.markers <- readRDS(file_markers_combined_all)

# Most expressed genes in each cluster
highest.expressed <- getHighestExpressedGenes(scRNA)
highest.expressed

# Get all diff genes in each cluster
FDR = 0.01
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

# Identify cell types as we did with flow cytometry

DefaultAssay(scRNA) <- "RNA"
scRNA <- NormalizeData(scRNA)
#scRNA <- AddModuleScore(scRNA, markers.core, name = "CellType")

#inds <- grep("CellType", colnames(scRNA@meta.data))
#colnames(scRNA@meta.data)[inds] <- paste0("CellType", names(markers.core))

#FeaturePlot(scRNA, features = colnames(scRNA@meta.data)[inds], order = TRUE, 
#            ncol = 6, cols = c("purple", "yellow"))

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

# Cluster identification -- This code is old and may not work anymore
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


##### GSEA - genotypes #####

Idents(scRNA) <- scRNA$genotype
DefaultAssay(scRNA) <- "RNA"

genotypes <- c("5XFAD", "PS19", "PS-5X")

for (G in genotypes) {
  markers <- FindMarkers(scRNA, ident.1 = G, ident.2 = "WT", 
                         logfc.threshold = 0.0, min.pct = 0, test.use = "MAST",
                         latent.vars = c("nCount_RNA", "StressScoreVariable"))
  saveRDS(markers, file = file.path(dir_allcells_go, paste0("markers_", G, "vWT.rds")))
}

# Only use genes expressed in at least 10 cells for GSEA. Gives 14840 genes
tmp <- GetAssayData(scRNA, slot = "counts")
tmp <- tmp > 0
ok <- rowSums(tmp)

good_genes <- names(ok[ok >= 10])
rm(tmp, ok)

genes.exclude = readRDS(file_excluded_genes)
good_genes = setdiff(good_genes, unlist(genes.exclude))

as.data.frame(msigdbr_collections())
reac <- msigdbr(species = "mouse", category = "C2", subcategory = "CP:REACTOME")
go_bp <- msigdbr(species = "mouse", category = "C5", subcategory = "GO:BP")
imm <- msigdbr(species = "mouse", category = "C7", subcategory = "IMMUNESIGDB")
hall <- msigdbr(species = "mouse", category = "H")

gene_set <- go_bp
gene_set_list <- split(x = gene_set$ensembl_gene, f = gene_set$gs_name)
set_name <- "GO:BP"
replace_string <- "GOBP" # "HALLMARK", "REACTOME", or "GOBP". ImmuneSigDB doesn't have anything to remove

# TODO I've just been re-running this with different lists, need to automate for
# final code distribution
sheets = list()

for (G in genotypes) {
  markers <- readRDS(file.path(dir_allcells_go, paste0("markers_", G, "vWT.rds")))
  
  ranks <- markers[good_genes, "avg_log2FC"]
  names(ranks) <- geneNameToEnsembl(good_genes)
  ranks <- ranks[abs(ranks) > 0.001] # remove near-0 values
  ranks <- ranks[order(ranks)]

  set.seed(112233)
  
  res <- fgsea(pathways = gene_set_list, stats = ranks, maxSize = 500)
  res <- subset(res, !is.na(padj))
  res <- subset(res, padj <= 0.05)
  
  collapsed <- collapsePathways(res, gene_set_list, ranks)

  res$leadingEdgeGeneSymbols <- lapply(res$leadingEdge, ensemblToGeneName)
  res$leadingEdgeGeneSymbols <- sapply(res$leadingEdgeGeneSymbols, str_c, collapse = ", ")
  res$leadingEdge <- sapply(res$leadingEdge, str_c, collapse = ", ")
  
  gs_source <- subset(gene_set, gs_name %in% res$pathway)
  gs_source <- unique(gs_source[,c("gs_name", "gs_exact_source")])
  
  res <- merge(res, gs_source, by.x = "pathway", by.y = "gs_name")
  
  res <- res[order(res$NES, decreasing = TRUE),]
  
  res$IsMainPathway <- res$pathway %in% collapsed$mainPathways
  
  res$pathway <- str_replace(res$pathway, replace_string, "") %>% 
                  str_replace_all("_", " ") %>% str_to_title()
  
  res.display <- res[res$IsMainPathway]
  res.display$pathway <- str_replace(res.display$pathway, replace_string, "") %>% 
                          str_replace_all("_", " ") %>% str_to_title()
  
  sheets[[paste0(G, " vs WT")]] <- res
  
  # TODO refine and save plots to disk
  plt <- ggplot(res.display, aes(reorder(pathway, NES), NES, fill = NES)) +
            geom_col() +
            coord_flip() +
            labs(x="Pathway", y="Normalized Enrichment Score",
                 title=paste0(set_name, " pathways NES from GSEA -- ", G, " vs WT")) + 
            theme_classic() + theme(axis.text.y = element_text(size = 5))
  print(plt)
}

write_xlsx(sheets, path = file.path(dir_allcells_go, "gsea_genotypes_gobp.xlsx"))

for (N in names(sheets)) {
  res <- sheets[[N]]
  #res <- subset(res, IsMainPathway == TRUE)
  
  # Have to add some dummy columns for cytoscape
  res2 <- data.frame("NAME" = res$gs_exact_source,
                     "GS<br> follow link to MSigDB" = res$pathway,
                     "GS DETAILS" = "",
                     "SIZE" = res$size, 
                     "ES" = res$ES,
                     "NES" = res$NES,
                     "NOM p-val" = res$pval,
                     "FDR q-val" = res$padj,
                     "FWER p-val" = res$padj,
                     "RANK AT MAX" = 0,
                     "LEADING EDGE" = res$leadingEdgeGeneSymbols,
                     check.names = FALSE
  )
  
  write.table(subset(res2, NES > 0), 
              file = file.path(dir_allcells_go, paste0("cytoscape_gsea_", N, "_up.txt")), 
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(subset(res2, NES < 0), 
              file = file.path(dir_allcells_go, paste0("cytoscape_gsea_", N, "_down.txt")), 
              sep = "\t", quote = FALSE, row.names = FALSE)
}


scRNA <- readRDS(file_seurat_analyzed_allcells)
scRNA.cd8 <- readRDS(file_seurat_analyzed_cd8)
scRNA.cd4 <- readRDS(file_seurat_analyzed_cd4)

scRNA$clusters <- "Unassigned"
scRNA$clusters[colnames(scRNA.cd8)] <- as.character(scRNA.cd8$clusters)
scRNA$clusters[colnames(scRNA.cd4)] <- as.character(scRNA.cd4$clusters)

Idents(scRNA) <- scRNA$clusters
DimPlot(scRNA, shuffle = TRUE, label = TRUE)
