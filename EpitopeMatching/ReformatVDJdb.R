# Reformats the raw export of VDJdb data into the format recognized by TCRMatch. 
# VDJdb export was generated using https://vdjdb.cdr3.net/search with the
# following parameters:
#   CDR3 tab
#       Species             - Human, Monkey, Mouse
#       Gene (chain)        - TRA, TRB
#       Variable segments   - (blank)
#       Joining segments    - (blank)
#       Sequence or pattern - (blank)
#       CDR3 length         - 5 to 30
#       By Levenstein distance - (blank)
#       Substitutions       - 0
#       Insertions          - 0
#       Deletions           - 0
#   Antigen tab - Everything left blank
#   MHC tab
#       Class     - MHCI, MHCII
#       Haplotype - (blank)
#   Meta tab
#       References - (blank)
#       Assay type - all boxes checked
#       Sequencing - all boxes checked
#       Minimal confidence score - 0
#       Spurious CDR3 - no boxes checked
#
# Search results were exported to a tsv file with paired gene export enabled, 
# downloaded on Jan 23, 2022.
# That tsv file has a lot of unnecessary extra fields/data, so this script pares
# the data down to only the fields required by TCRMatch: trimmed_seq, 
# original_seq, receptor_group, epitopes, source_organisms, source_antigens
#
# This script trims the starting "C" and ending "F/W" from each sequence, and 
# sequences that are < 5 AA in length were filtered out prior to downloading
# the data. VDJDB does not have a "source antigen" or "receptor group" field so  
# these columns are left blank in the final cleaned data. 
#
# TCRMatch source/reference: https://github.com/IEDB/TCRMatch
#
# Author: Jaclyn Beck
# Last revision: Jan 27, 2022

library(stringr)
library(dplyr)
source("Filenames.R")
source(file.path("EpitopeMatching", "TrimCDR3.R"))

orig.data <- read.table(file_vdjdb_export, header = TRUE, sep = "\t", comment.char = "")

# Remove some unneeded columns
needed <- c("complex.id", "Gene", "CDR3", "Epitope", "Epitope.species")
orig.data <- orig.data[,needed]

# Trim
orig.data$CDR3.Trimmed <- TrimCDR3(orig.data$CDR3)

# Split into alpha and beta TCR data, put in the order expected by TCRMatch
alpha.pre <- subset(orig.data, Gene == "TRA")
beta.pre <- subset(orig.data, Gene == "TRB")

alpha.data <- data.frame(trimmed_seq = alpha.pre$CDR3.Trimmed,
                         original_seq = alpha.pre$CDR3,
                         receptor_group = "",
                         epitopes = alpha.pre$Epitope, 
                         source_organisms = alpha.pre$Epitope.species, 
                         source_antigens = "")

alpha.data <- distinct(alpha.data)

beta.data  <- data.frame(trimmed_seq = beta.pre$CDR3.Trimmed,
                         original_seq = beta.pre$CDR3,
                         receptor_group = "",
                         epitopes = beta.pre$Epitope, 
                         source_organisms = beta.pre$Epitope.species, 
                         source_antigens = "")

beta.data <- distinct(beta.data)

write.table(alpha.data, file_vdjdb_clean_alpha, sep="\t", row.names = FALSE, 
            quote = FALSE)

write.table(beta.data, file_vdjdb_clean_beta, sep="\t", row.names = FALSE, 
            quote = FALSE)

rm(list=ls())
