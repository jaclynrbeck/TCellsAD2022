# Helper functions used by multiple processing steps. 

# Author: Jaclyn Beck
# Final script used for paper as of Sep 02, 2022



# Returns a Seurat object that only contains cells that are "positive" for a 
# specific gene (i.e. RNA count >= threshold). 
# pattern: a regular expression used for a grep search of gene names.
#          Example: "Cd3g" will only find the gene "Cd3g", while "Cd3[deg]"
#                   will find "Cd3d", "Cd3e" and "Cd3g".
# In cases where multiple genes are specified, the cell must express at least
# <threshold> of one of those genes to pass. 
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


# Converts gene names to Ensembl IDs. This function is necessary because
# gene names may have "--1" appended to them if they map to multiple 
# Ensembl IDs, and the "gene symbols" file accounts for this. 
geneNameToEnsembl <- function(genes) {
  features <- readRDS(file_gene_symbols)
  rownames(features) <- features$Gene.Symbol
  features[genes, "Ensembl.Id"]
}

# Converts Ensembl IDs to gene names. This function is necessary because
# gene names may have "--1" appended to them if they map to multiple 
# Ensembl IDs, and the "gene symbols" file accounts for this. 
ensemblToGeneName <- function(genes) {
  features <- readRDS(file_gene_symbols)
  rownames(features) <- features$Ensembl.Id
  features[genes, "Gene.Symbol"]
}
