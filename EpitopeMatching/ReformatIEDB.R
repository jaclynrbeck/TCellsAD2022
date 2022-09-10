# Reformats the raw export of IEDB data into the format recognized by TCRMatch. 
# IEDB export was generated using the http://www.iedb.org/ search function with
# the following parameters:
#   Epitope         - Linear peptide
#   Epitope Source  - (blank)
#   Host            - Any
#   Assay           - T Cell only
#   Assay (Outcome) - Positive only
#   MHC Restriction - Any
#   Disease         - Any
#
# Search results were exported to a csv file, downloaded on Jan 23, 2022.
# That csv file has a lot of unnecessary extra fields/data, so this script pares
# the data down to only the fields required by TCRMatch: trimmed_seq, 
# original_seq, receptor_group, epitopes, source_organisms, source_antigens
#
# This script trims the starting "C" and ending "F/W" from each sequence, and 
# removes sequences that are < 5 AA in length.
#
# TCRMatch source/reference: https://github.com/IEDB/TCRMatch
#
# Author: Jaclyn Beck
# Last revision: Jan 27, 2022

library(stringr)
library(dplyr)
source("Filenames.R")
source(file.path("EpitopeMatching", "TrimCDR3.R"))

orig.data <- read.csv(file_iedb_export)

# Remove some initial unneeded columns
needed <- c("Group.Receptor.ID", "Description", "Antigen", "Organism", 
            "Chain.1.CDR3.Curated", "Chain.1.CDR3.Calculated", 
            "Chain.2.CDR3.Curated", "Chain.2.CDR3.Calculated")
orig.data <- orig.data[,needed]

# Remove things like "+ SCM(F5)" from epitopes (i.e. LLFGFPVYV + SCM(F5))
orig.data$Description <- str_replace(orig.data$Description, " \\+ .*", "")

# Find all the epitopes which have a recognized organism/antigen
nonblank.data <- subset(orig.data, Organism != "")
nonblank.data <- nonblank.data[,c("Description", "Organism", "Antigen")]
nonblank.data <- distinct(nonblank.data)
epitopes <- unique(nonblank.data$Description)

orig.data <- subset(orig.data, Description %in% epitopes)

# Fill in missing organism/antigen blanks for these epitopes. This is a very
# small number of samples (2 as of the writing of this script)
blank.org <- which(orig.data$Organism == "")

for (X in blank.org) {
  epi <- subset(nonblank.data, Description == orig.data$Description[X])
  
  # There are 1-2 cases where nonblank.data has the same epitope map to multiple
  # organisms/antigens (i.e. mouse vs rat). Here we just use the first one. 
  orig.data$Organism[X] <- epi$Organism[1]
  orig.data$Antigen[X] <- epi$Antigen[1]
}

orig.data <- distinct(orig.data)
rm(nonblank.data, blank.org, epitopes)

# Replace any blank curated CDR3s with calculated CDR3s. This only affects a
# small fraction of entries and allows removal of essentially duplicate rows
blank.alpha <- which(orig.data$Chain.1.CDR3.Curated == "")
blank.beta <- which(orig.data$Chain.2.CDR3.Curated == "")

orig.data$Chain.1.CDR3.Curated[blank.alpha] = orig.data$Chain.1.CDR3.Calculated[blank.alpha]
orig.data$Chain.2.CDR3.Curated[blank.beta] = orig.data$Chain.2.CDR3.Calculated[blank.beta]

orig.data <- distinct(orig.data)
rm(blank.alpha, blank.beta)

orig.data$CDR3.1.Trimmed <- TrimCDR3(orig.data$Chain.1.CDR3.Curated)
orig.data$CDR3.2.Trimmed <- TrimCDR3(orig.data$Chain.2.CDR3.Curated)

# Split into alpha and beta TCR data, put in the order expected by TCRMatch
alpha.data <- data.frame(trimmed_seq = orig.data$CDR3.1.Trimmed,
                         original_seq = orig.data$Chain.1.CDR3.Curated,
                         receptor_group = orig.data$Group.Receptor.ID,
                         epitopes = orig.data$Description, 
                         source_organisms = orig.data$Organism, 
                         source_antigens = orig.data$Antigen)

alpha.data <- subset(alpha.data, trimmed_seq != "")

len <- unlist(lapply(alpha.data$trimmed_seq, nchar))
alpha.data <- alpha.data[len >= 5,]

alpha.data <- distinct(alpha.data)

beta.data  <- data.frame(trimmed_seq = orig.data$CDR3.2.Trimmed,
                         original_seq = orig.data$Chain.2.CDR3.Curated,
                         receptor_group = orig.data$Group.Receptor.ID,
                         epitopes = orig.data$Description, 
                         source_organisms = orig.data$Organism, 
                         source_antigens = orig.data$Antigen)

beta.data <- subset(beta.data, trimmed_seq != "")

len <- unlist(lapply(beta.data$trimmed_seq, nchar))
beta.data <- beta.data[len >= 5,]

beta.data <- distinct(beta.data)

# Remove invalid amino acid characters from CDR3 chains
alphabet <- "ARNDCQEGHILKMFPSTWYV"
rma <- grep(paste0("[^", alphabet, "]"), alpha.data$original_seq)
rmb <- grep(paste0("[^", alphabet, "]"), beta.data$original_seq)

alpha.data <- alpha.data[-rma,]
beta.data <- beta.data[-rmb,]

write.table(alpha.data, file_iedb_clean_alpha, sep="\t", row.names = FALSE, 
            quote = FALSE)

write.table(beta.data, file_iedb_clean_beta, sep="\t", row.names = FALSE, 
            quote = FALSE)

rm(list=ls())
