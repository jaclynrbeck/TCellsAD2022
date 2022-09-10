# Reads the output of TCRMatch (TCRa and TCRb were scored separately). Currently
# all matches with a score >= 0.9 are returned. Then, for each clonotype present
# in the data set, it looks for the best score for the clonotype's alpha and 
# beta receptors using the following rules:
#   If both the TRA and TRB match to the same epitope, their scores are averaged
#     for that epitope
#   If only one of TRA or TRB match an epitope, the epitope's score is the score
#     of the receptor that matched it
#   Take the highest score of all epitopes, ignoring whether it was a paired or
#     single match
# The score and type of match (TRA, TRB, Both, or None) are recorded so they can
# be used to filter data later. 
#
# Note to self: The written files have to be TSVs because of commas in some of
# the fields. 
#
# Author: Jaclyn Beck
# Last revision: Sep 09, 2022

library(stringr)
library(dplyr)
library(reshape2)

source("Filenames.R")
source(file.path("EpitopeMatching", "TrimCDR3.R"))

ReadMatchOutput <- function(filename) {
  # Read in the output from TCRMatch
  scores <- read.table(filename, header = TRUE, sep = "\t", comment.char = "")
  
  # Get rid of unused data to pare down the size
  keep <- c("input_sequence", "epitope", "score", "organism")
  scores <- scores[,keep]
  scores <- distinct(scores)
  
  # Sometimes the same cdr3 will get matched to the same epitope with a slightly
  # different score (from merging 2 different databases), so use the highest 
  # score per cdr3/epitope combination
  scores <- aggregate(score ~ input_sequence + epitope + organism, scores, max)
  return(scores)
}

alpha.scores <- ReadMatchOutput(file_combined_alpha)
beta.scores <- ReadMatchOutput(file_combined_beta)

# Read clonotype data from CellRanger
clono <- read.csv(file_clonotypes_summary)

# Clonotypes in this file are concatenated like "TRB:XXXXX;TRA:XXXXXX".
# Separate them into TRA and TRB but keep the pairs associated. Cells with 
# multiple TRBs are assumed to be doublets and are skipped. Cells missing a TRA 
# or TRB are scored based on the one they do have. Cells with multiple TRAs
# will use the highest-scoring TRA.
uniq <- distinct(clono[,c("clonotype_id", "cdr3s_aa")])
unique.clonos <- data.frame(Clonotype = uniq$cdr3s_aa, 
                            ClonotypeId = uniq$clonotype_id, 
                            TRA1 = "", TRA2 = "", TRB = "", Score = 0)
rm(uniq)

combined <- str_split(unique.clonos$Clonotype, ";")

for (R in 1:nrow(unique.clonos)) {
  clono.list <- combined[[R]]
  
  tra <- grep("TRA:", clono.list, value = TRUE)
  trb <- grep("TRB:", clono.list, value = TRUE)
  
  # Cells with exactly one TRB or up to 2 TRAs get these fields filled in. 
  # Otherwise the TRA and TRB are left blank for the cell.
  if (length(trb) <= 1) {
    if (length(tra) <= 2 & length(tra) >= 1) {
      unique.clonos$TRA1[R] <- str_replace(tra[1], "TRA:", "")
      
      if (length(tra) == 2) {
        unique.clonos$TRA2[R] <- str_replace(tra[2], "TRA:", "")
      }
    }
    if (length(trb) == 1 & length(tra) <= 2) {
      unique.clonos$TRB[R] <- str_replace(trb, "TRB:", "")
    }
  }
}

rm(combined)

# Get rid of doublet cells, which would have empty entries for both TRA1 and TRB
unique.clonos <- subset(unique.clonos, TRA1 != "" | TRB != "")

# Trim sequences to match what's in the output file
unique.clonos$TRA1 <- TrimCDR3(unique.clonos$TRA1)
unique.clonos$TRA2 <- TrimCDR3(unique.clonos$TRA2)
unique.clonos$TRB <- TrimCDR3(unique.clonos$TRB)
unique.clonos$Chain.Match <- "None"
unique.clonos$Epitope <- ""
unique.clonos$Antigen <- ""

# Match clonotypes to scores in the TCRMatch output
for (R in 1:nrow(unique.clonos)) {
  tra <- c(unique.clonos$TRA1[R], unique.clonos$TRA2[R])
  trb <- unique.clonos$TRB[R]

  matches.a <- subset(alpha.scores, input_sequence %in% tra)
  matches.b <- subset(beta.scores, input_sequence == trb)
  
  merged <- merge(matches.a, matches.b, by = "epitope", all = TRUE)
  if (nrow(merged) == 0) {
    next
  }
  merged$match <- ""
  merged$organism <- ""
    
  # This will only catch cases where both TRA and TRB match an epitope
  merged$avg.score <- (merged$score.x + merged$score.y) / 2
  both <- !is.na(merged$avg.score)
  merged$match[both] <- "Both"
  merged$organism[both] <- merged$organism.x[both]
  
  # TRA or TRB matches only
  tra.only <- which(is.na(merged$input_sequence.y))
  trb.only <- which(is.na(merged$input_sequence.x))
  
  merged$avg.score[tra.only] <- merged$score.x[tra.only]
  merged$match[tra.only] <- "A"
  merged$organism[tra.only] <- merged$organism.x[tra.only]
  
  merged$avg.score[trb.only] <- merged$score.y[trb.only]
  merged$match[trb.only] <- "B"
  merged$organism[trb.only] <- merged$organism.y[trb.only]
  
  best <- subset(merged, avg.score == max(merged$avg.score))
  
  # There might be multiple matches, so these are concatenated
  unique.clonos$Chain.Match[R] <- paste(unique(best$match), collapse = " / ")
  unique.clonos$Score[R] <- best$avg.score[1]
  unique.clonos$Epitope[R] <- paste(best$epitope, collapse = " / ")
  unique.clonos$Antigen[R] <- paste(unique(best$organism), collapse = " / ")
}

# Fixes some duplicated antigen names
antigen <- str_split(unique.clonos$Antigen, pattern = " / ")
antigen <- lapply(antigen, FUN = unique) %>% sapply(FUN = paste, collapse = " / ")
unique.clonos$Antigen <- antigen

sum(unique.clonos$Score >= 0.97)
sum(unique.clonos$Score >= 0.95)

# Add some information to this data frame for the file
tcr.anno <- read.csv(file_clonotypes_processed)
table(tcr.anno$Genotype)

clono <- read.csv(file_clonotypes_summary)
clono <- clono[,-c(6,7)] # Get rid of inkt_evidence and mait_evidence

colnames(unique.clonos)[2] <- "clonotype_id"
unique.clonos <- unique.clonos[,-1]

clono <- subset(clono, clonotype_id %in% unique.clonos$clonotype_id)
tcr.match <- subset(tcr.anno, ClonotypeId %in% unique.clonos$clonotype_id)

genos <- tcr.match[,c("ClonotypeId", "Sample", "Genotype")]
colnames(genos) <- c("clonotype_id", "Sample", "Genotype")
genos <- distinct(genos)

clono <- merge(clono, genos, by = "clonotype_id")
clono <- merge(clono, unique.clonos, by = "clonotype_id", all.y = TRUE)

clono <- arrange(clono, desc(frequency))

# Keep all scores in this file. 
write.table(clono, file_matched_epitopes, sep="\t", row.names = FALSE, quote = FALSE)

# Make a file of probable matches, sorted by frequency
good <- subset(clono, Score >= 0.97)

write.table(good, 
            file.path(dir_epitope_output, "ProbableMatchedClonotypes.tsv"),
            sep="\t", row.names = FALSE, quote = FALSE)

# File of non-matches
bad <- subset(clono, Score < 0.97)

write.table(bad, 
            file.path(dir_epitope_output, "UnmatchedClonotypes.tsv"),
            sep="\t", row.names = FALSE, quote = FALSE)

rm(list=ls())
gc()

