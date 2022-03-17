# Output from CellRanger aggr doesn't carry genotype over into its barcodes. 
# It only increases the number at the end of the barcode after each sample it 
# adds to the aggregated data.
#
# This script re-writes the "barcodes.tsv.gz" file with the genotype embedded in
# the barcodes for easy reading into Seurat. It also does this for the 
# clonotypes file, and simplifies the clonotypes file a bit for easier data 
# manipulation.
#
# Lastly, gene names in the "features.tsv.gz" file aren't necessarily unique,
# which leads Seurat to ignore duplicate names. This is fixed in this script as
# well by adding "--1" to duplicated gene names and re-writing the file. 

# Author: Jaclyn Beck
# Final script for paper as of March 06, 2022

library(stringr)

source("Filenames.R")

file_barcodes <- file.path(dir_filtered_counts, "barcodes.tsv.gz")
file_features <- file.path(dir_filtered_counts, "features.tsv.gz")

# Copy the original barcodes and features files as a backup
file.copy(from = file_barcodes, 
          to = file.path(dir_filtered_counts, "barcodes.original.tsv.gz"))

file.copy(from = file_features, 
          to = file.path(dir_filtered_counts, "features.original.tsv.gz"))

# Read in barcodes and annotation
barcodes <- read.table(file = gzfile(file_barcodes))
barcodes <- barcodes$V1

geno.anno <- read.csv(file = file.path(dir_cellranger, "aggregation.csv"))
genotypes <- geno.anno$sample_id

# Replace the number at the end of each barcode with the genotype. Also change 
# the "-" separator to "_" so the genotype "PS-5X" is read correctly by Seurat
for (N in 1:length(genotypes)) {
  barcodes <- str_replace(barcodes, paste0("-", N), paste0("_", genotypes[N]))
}

# Write a new 'barcodes.tsv.gz' file for Read10X
write.table(barcodes, 
            file=gzfile(file_barcodes), sep=' ', row.names=FALSE, 
            col.names=FALSE, quote=FALSE)


## Rearrange clonotype data

anno <- read.csv(file_clonotypes_annotated)
clono <- read.csv(file_clonotypes_raw)

anno$clonotype <- paste(anno$chain, anno$cdr3, sep=": ")

# Put the data in an easier format. Combine the TRA and TRB entries for each 
# cell, and discard columns that aren't needed.
tcr.tmp <- data.frame(Barcode = anno$barcode,
                      Genotype = anno$origin,
                      TRA.V = " ",
                      TRA.J = " ",
                      TRA.C = " ",
                      TRB.V = " ",
                      TRB.D = " ",
                      TRB.J = " ",
                      TRB.C = " ",
                      Clonotype = " ", 
                      ClonotypeId = anno$raw_clonotype_id,
                      iNKT = FALSE,
                      MAIT = FALSE)
tcr.tmp <- unique(tcr.tmp)

for (R in 1:nrow(tcr.tmp)) {
  tcrs <- subset(anno, barcode == tcr.tmp$Barcode[R])
  tcr.tmp[R, "Clonotype"] <- str_flatten(sort(unique(tcrs$clonotype)), collapse=", ")
  
  # There might be more than one of each of these
  tra <- subset(tcrs, chain == "TRA")
  trb <- subset(tcrs, chain == "TRB")
  
  if (nrow(tra) > 0) {
    tcr.tmp[R, "TRA.V"] <- str_flatten(sort(unique(tra$v_gene)), collapse=", ")
    tcr.tmp[R, "TRA.J"] <- str_flatten(sort(unique(tra$j_gene)), collapse=", ")
    tcr.tmp[R, "TRA.C"] <- str_flatten(sort(unique(tra$c_gene)), collapse=", ")
  }
  
  if (nrow(trb) > 0) {
    tcr.tmp[R, "TRB.V"] <- str_flatten(sort(unique(trb$v_gene)), collapse=", ")
    tcr.tmp[R, "TRB.D"] <- str_flatten(sort(unique(trb$d_gene)), collapse=", ")
    tcr.tmp[R, "TRB.J"] <- str_flatten(sort(unique(trb$j_gene)), collapse=", ")
    tcr.tmp[R, "TRB.C"] <- str_flatten(sort(unique(trb$c_gene)), collapse=", ")
  }
}

# Make the barcode compatible with the new barcodes above
tcr.tmp$Sample <- str_replace(tcr.tmp$Barcode, "-[0-9]", 
                              paste0("_", tcr.tmp$Genotype))
tcr.tmp <- unique(tcr.tmp)

# Get iNKT and MAIT info
for (R in 1:nrow(clono)) {
  id <- clono$clonotype_id[R]
  inds <- which(tcr.tmp$ClonotypeId == id)
  tcr.tmp$iNKT[inds] <- (clono$inkt_evidence[R] != "")
  tcr.tmp$MAIT[inds] <- (clono$mait_evidence[R] != "")
}

write.csv(tcr.tmp, file = file_clonotypes_processed, row.names = FALSE)

# Fix the features file
features <- read.table(file = gzfile(file_features))
features$V2 <- make.unique(features$V2, sep = "--")

write.table(features, 
            file=gzfile(file_features), sep='\t', row.names=FALSE, 
            col.names=FALSE, quote=FALSE)

features <- features[,1:2]
colnames(features) <- c("Ensembl.Id", "Gene.Symbol")
saveRDS(features, file_gene_symbols)

# Clear data
rm(list=ls())

# Done. 
