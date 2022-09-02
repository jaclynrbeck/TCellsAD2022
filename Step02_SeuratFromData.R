# Reads 10X data from the output of cellranger aggr and turns it into a Seurat
# object that has been normalized and scaled. 

# Author: Jaclyn Beck
# Final script used for paper as of April 20, 2022

library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
source(file.path("functions", "SeuratFromData_HelperFunctions.R"))
source(file.path("functions", "Analysis_HelperFunctions.R"))
source("Filenames.R")


##### Create the Seurat object #####

scRNA <- makeSeurat(dir_filtered_counts)

##### Doublet identification #####

# Identify potential doublets from TCR data -- looking for multiple beta-chains
# only. Some T cells can express two alpha chains so we ignore that. 
tcr.anno <- read.csv(file_clonotypes_processed)
table(tcr.anno$Genotype)

# Genes in these columns are comma-separated if there's more than one. 
# Split on the comma so we can count how many genes there are
genes <- sapply(c("TRB.V", "TRB.D", "TRB.J", "TRB.C"),
                function(X) { str_split(tcr.anno[,X], pattern = ", ") })

# For each cell, get the length of the list of genes in each gene column.
# Returns true if all lengths are 1 or less. False otherwise, indicates a 
# doublet.
tcr.anno$Singlet <- sapply(1:nrow(genes), function(X) {
  lengths <- sapply(genes[X,], length)
  all(lengths <= 1)
})
(table(tcr.anno$Singlet) / nrow(tcr.anno)) * 100 # Doublet rate is ~4%

# Add this info to the Seurat object. Any cells without a recognized clonotype
# are marked as "Unknown" whether they're a singlet or doublet
scRNA$Singlet <- "Unknown"
scRNA$Singlet[tcr.anno$Sample] <- tcr.anno$Singlet
(table(scRNA$Singlet) / ncol(scRNA)) * 100 

##### Removal of CD3- cells and doublets #####

# Only take Cd3+ cells or cells with a recognized TCR
pos.cells <- filterOnGenePositive(scRNA, "Cd3[edg]$", 1)
scRNA <- subset(scRNA, cells = unique(c(pos.cells, tcr.anno$Sample)))
scRNA <- subset(scRNA, Singlet != FALSE)

scRNA <- removeLowExpressedGenes(scRNA, 10)

##### Exploratory plots ##### 

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

##### Remove low-quality cells #####

scRNA <- subset(scRNA, subset = nFeature_RNA >= 800 & 
                  percent.mt <= 5 &   # Clear cutoff at 4 or 5
                  percent.rb >= 8 &   # Clear cutoff at 8-10
                  percent.rb <= 55)   # Removes cells with low RNA + high rb

VlnPlot(scRNA, features = "nFeature_RNA", pt.size = 0) + 
  geom_hline(yintercept = 2*median(scRNA$nFeature_RNA))

VlnPlot(scRNA, features = "nCount_RNA", pt.size = 0.1) + 
  geom_hline(yintercept = 1000) +
  geom_hline(yintercept = 3*median(scRNA$nCount_RNA))

FeatureScatter(scRNA, feature1 = "complexity", feature2 = "nCount_RNA", shuffle=TRUE, pt.size = 0.1) +
  scale_y_log10(n.breaks=16) + scale_x_continuous(n.breaks=16) + theme_bw()

VlnPlot(scRNA, features = "complexity", pt.size = 0.1) + 
  geom_hline(yintercept = 4.75, col = "red")

# Remove potential doublets by thresholding on complexity
scRNA <- subset(scRNA, complexity <= 4.75 & nCount_RNA <= 20000)

##### Finalize un-normalized Seurat object #####

table(scRNA$orig.ident)
scRNA <- removeLowExpressedGenes(scRNA, 10)
scRNA <- DietSeurat(scRNA, scale.data = FALSE)
scRNA$genotype <- str_replace(scRNA$orig.ident, "-[12]", "")
scRNA$batch <- "Batch1"
scRNA$batch[scRNA$orig.ident == "WT-2" | scRNA$orig.ident == "5XFAD-2"] <- "Batch2"

# At this point scRNA contains all cells that passed QC
saveRDS(scRNA, file = file_seurat_unnorm)


##### Normalize and score #####

# Account for batch effect of cellular stress -- 
# After input of diff genes between batches into PANTHER, and looking at 
# Reactome pathways, the top results all relate to cellular stress. We use the 
# genes in the "cellular response to heat stress" to create a "stress score".
# This score is regressed out during feature scaling. 

scRNA <- readRDS(file_seurat_unnorm)
scRNA <- NormalizeData(scRNA)
scRNA <- FindVariableFeatures(scRNA, nfeatures = 4000)

reac <- read.table(file.path(dir_data, "External", "reactome_ResponseToHeatStress.tsv"), 
                   sep = "\t", skip = 4, header = TRUE)

reac_stress <- reac$Gene.Name[reac$Gene.Name %in% rownames(scRNA)]

assay <- GetAssayData(scRNA, slot = "counts")
assay <- assay[reac_stress,]
assay <- rowSums(assay > 0)

# Only look at genes expressed in > 5% of cells
reac_stress <- reac_stress[assay >= ncol(scRNA) * 0.05]

# Score on genes that are variable, to help discriminate between clusters
reac_v <- reac_stress[reac_stress %in% VariableFeatures(scRNA)]

scRNA <- AddModuleScore(scRNA, features = list(reac_stress, reac_v), name = "Stress")
nmeta <- ncol(scRNA@meta.data)
colnames(scRNA@meta.data)[(nmeta-1):nmeta] <- c("StressScore", "StressScoreVariable")

RidgePlot(scRNA, features = c("StressScore", "StressScoreVariable"), 
          group.by = "orig.ident", same.y.lims = TRUE)

##### Scale and integrate w/ regression #####

scRNA <- generateIntegratedData(scRNA, method = "SCT",
                                regress = c("nCount_RNA", "StressScoreVariable"))

saveRDS(scRNA, file = file_seurat_norm)


##### Create a list of genes to ignore in future analysis #####

DefaultAssay(scRNA) <- "RNA"

all.genes <- rownames(scRNA)

# For downstream analysis / visualization we can ignore mitochondrial, 
# ribosomal, and heat stress genes. These were allowed in for integration 
# because they might be biologically relevant to cell state, but for 
# visualization we want to focus more on genes related to immune response and
# less on variables that may be due to batch effect.
genes.exclude <- list("Ribosomal" = grep("Rp[s|l]", all.genes, value=TRUE),
                      "Mitochondrial" = c(grep("^mt", all.genes, value=TRUE),
                                          grep("Mrp[s|l]", all.genes, value = TRUE)),
                      "Stress" = all.genes[geneNameToEnsembl(all.genes) %in% reac$Gene.ID],
                      # Other stress-related and X-related genes that are sig. 
                      # different between batches. 
                      "Batch" = c("Dnajb1", "Hsph1", "Hspa1a", "Hspa8", 
                                  "Hsp90ab1", "Hspd1", "Dnaja1", "Chordc1", 
                                  "Dnajc15", "Dynll1", "Fkbp4", "Hsp90aa1", 
                                  "Vim", "Hist1h1c", "Gm11808", "H3f3b", 
                                  "Hist1h1d", "Hist1h4d", "Stip1", "H1f0", 
                                  "Hmgb2", "Fos", "Hspe1", "Uba52", "Hist1h2ap",
                                  "Cul3", "Aven", "Ywhaq", "Xist"))

saveRDS(genes.exclude, file_excluded_genes)

gc()

# Done
