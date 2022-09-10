# This script extracts CDR3 sequences from clonotype data output by CellRanger,
# and separates them into TRA and TRB data sets. CDR3 sequences are then written
# to a text file for use with TCRMatch. The initial "C" and ending "F/W" are 
# trimmed from each sequence.
#
# Author: Jaclyn Beck
# Last revision: Jan 27, 2021

library(stringr)
library(dplyr)
source("Filenames.R")
source(file.path("EpitopeMatching", "TrimCDR3.R"))

clono <- read.csv(file_clonotypes_annotated)
clono <- clono[,c("chain", "cdr3")]
clono <- distinct(clono)

clono$cdr3 <- TrimCDR3(clono$cdr3)

# Separate into alpha and beta sets to write to separate files
alpha <- subset(clono, chain == "TRA")
alpha <- alpha$cdr3

beta <- subset(clono, chain == "TRB")
beta <- beta$cdr3

write.table(as.data.frame(alpha), file_epitope_split_tra, 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(as.data.frame(beta), file_epitope_split_trb, 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

rm(list=ls())



