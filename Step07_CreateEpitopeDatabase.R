# Reformats raw exports from IEDB and VDJdb into a format that TCRMatch can
# read. Same with the CellRanger clonotype output. 
# This file concatenates the alpha TCR data from IEDB and VDJdb as one file and
# the beta data from IEDB and VDJdb as another. Slight differences in the 
# naming of source organisms are corrected to make the databases consistent
# with each other. Duplicate entries between databases are removed. 
#
# After running this file, the TCRMatch program should be run from the
# command-line as follows:
# ./path/to/tcrmatch -s 0.90 -t 4 -i path/to/<dir_epitope_matching>/Parenchymal_TCRalpha.txt -d path/to/<dir_epitope_matching>/combined_epitope_db_alpha.tsv > path/to/<dir_epitope_output>/alpha_combined.tsv
# ./path/to/tcrmatch -s 0.90 -t 4 -i path/to/<dir_epitope_matching>/Parenchymal_TCRbeta.txt -d path/to/<dir_epitope_matching>/combined_epitope_db_beta.tsv > path/to/<dir_epitope_output>/beta_combined.tsv

#
# Author: Jaclyn Beck
# Last revision: June 23, 2022

library(stringr)
library(dplyr)

# This runs the code in each of these files to create the cleaned databases and
# trimmed sequences from our clonotypes. 
source(file.path("EpitopeMatching", "ReformatClonotypes.R"))
source(file.path("EpitopeMatching", "ReformatVDJdb.R"))
source(file.path("EpitopeMatching", "ReformatIEDB.R"))

source("Filenames.R")

db.files.a <- c(file_vdjdb_clean_alpha, file_iedb_clean_alpha)
db.files.b <- c(file_vdjdb_clean_beta, file_iedb_clean_beta)

db <- list()

# Combines the VDJdb and IEDB databases and corrects for slightly different
# naming between the two
for (db.files in list(db.files.a, db.files.b)) {

  df <- data.frame()

  for (F in db.files) {
    df <- rbind(df, read.table(F, header = TRUE, sep = "\t", comment.char = "",
                               quote = "\"", na.strings = ""))
  }
  
  # Receptor group and antigen are meaningless when databases are combined
  df$source_antigens <- ""
  df$receptor_group <- ""
  
  df <- distinct(df)
  
  # Rename some of the VDJdb organisms to match those in IEDB, and/or shorten
  # some names
  renames <- list("CMV" = "Human herpesvirus 5 (Human cytomegalovirus)",
                  "DENV1" = "Dengue virus 1 (Dengue virus type 1)",
                  "DENV2" = "Dengue virus 2 (Dengue virus type 2)",
                  "DENV3/4" = "Dengue virus 3 (Dengue virus serotype 3)",
                  "E.Coli" = "Escherichia coli",
                  "EBV" = "Human herpesvirus 4 (Epstein Barr virus)",
                  "GallusGallus" = "Gallus gallus (chicken)",
                  "HCV" = "Hepatitis C virus",
                  "HIV-1 M:B_HXB2R (Human immunodeficiency virus type 1 (HXB2 ISOLATE))" = "HIV-1",
                  "HomoSapiens" = "Homo sapiens (human)",
                  "HPV" = "Human papillomavirus type 16 (Human papilloma virus type 16)",
                  "HSV-2" = "Human herpesvirus 2",
                  "HTLV-1" = "Human T-cell leukemia virus type I",
                  "Human herpesvirus 5 strain AD169 (Human cytomegalovirus (strain AD169))" = "Human herpesvirus 5 (Human cytomegalovirus)",
                  "Human immunodeficiency virus 1 (human immunodeficiency virus 1 HIV-1)" = "HIV-1",
                  "Human respiratory syncytial virus A2 (Human respiratory syncytial virus (strain A2))" = "Human respiratory syncytial virus",
                  "Human T-cell leukemia virus type I (Human T cell leukemia virus type 1)" = "Human T-cell leukemia virus type I",
                  "Influenza A virus (A/bar-headed goose/Qinghai/3/2005(H5N1)) (Influenza A virus (A/bar headed goose/Qinghai/3/2005(H5N1)))" = "Influenza A virus (A/bar-headed goose/Qinghai/3/2005(H5N1))",
                  "Influenza A virus (A/chicken/Anhui/1/1998(H9N2)) (Influenza A virus (A/Chicken/Anhui/1/98(H9N2)))" = "Influenza A virus (A/chicken/Anhui/1/1998(H9N2))",
                  "Influenza A virus (A/Japan/305/1957(H2N2)) (Influenza A virus (A/Japan/305/57(H2N2)))" = "Influenza A virus (A/Japan/305/1957(H2N2))",
                  "Influenza A virus (A/Memphis/4/1973(H3N2)) (Influenza A virus (A/Memphis/4/73(H3N2)))" = "Influenza A virus (A/Memphis/4/1973(H3N2))",
                  "Influenza A virus (A/nt/60/1968(H3N2)) (Influenza A virus (A/NT/60/68/(H3N2)))" = "Influenza A virus (A/nt/60/1968(H3N2))",
                  "Influenza A virus (A/Puerto Rico/8/1934(H1N1)) (Influenza A virus (A/PR 8/34 (H1N1)))" = "Influenza A virus (A/Puerto Rico/8/1934(H1N1))",
                  "Influenza A virus (A/Puerto Rico/8/34/Mount Sinai(H1N1)) (Influenza A virus (A/Puerto Rico/8/1934(mt sinai)(H1N1)))" = "Influenza A virus (A/Puerto Rico/8/34/Mount Sinai(H1N1))",
                  "InfluenzaA" = "Influenza A virus",
                  "LCMV" = "Lymphocytic choriomeningitis mammarenavirus (Lymphocytic choriomeningitis virus)",
                  "M.tuberculosis" = "Mycobacterium tuberculosis",
                  "ManducaSexta" = "Manduca sexta (Carolina sphinx)",
                  "MCMV" = "Murid betaherpesvirus 1 (Mouse cytomegalovirus 1)",
                  "MCPyV" = "Merkel cell polyomavirus (MCPyV)",
                  "MusMusculus" = "Mus musculus (mouse)",
                  "PlasmodiumBerghei" = "Plasmodium berghei ANKA",
                  "PseudomonasAeruginosa" = "Pseudomonas aeruginosa",
                  "PseudomonasFluorescens" = "Pseudomonas fluorescens",
                  "RSV" = "Human respiratory syncytial virus",
                  "SaccharomycesCerevisiae" = "Saccharomyces cerevisiae (baker's yeast)",
                  "Salmonella enterica subsp. enterica serovar Typhi str. 404ty (Salmonella enterica subsp. enterica serovar Typhi 404ty)" = "Salmonella enterica subsp. enterica serovar Typhi str. 404ty", 
                  "SARS-CoV-2" = "SARS-CoV2",
                  "SelaginellaMoellendorffii" = "Selaginella moellendorffii",
                  "SIV" = "Simian immunodeficiency virus (Chimpanzee immunodeficiency virus)",
                  "StreptomycesKanamyceticus" = "Streptomyces kanamyceticus",
                  "synthetic" = "Synthetic",
                  "TriticumAestivum" = "Triticum aestivum (Canadian hard winter wheat)",
                  "VSV" = "Vesicular stomatitis virus (vesicular stomatitis virus VSV)",
                  "YFV" = "Yellow fever virus (Flavivirus febricis)")
  
  renames <- stack(renames)
  colnames(renames) <- c("New", "Old")
  for (R in 1:nrow(renames)) {
    rows <- which(df$source_organisms == renames$Old[R])
    df$source_organisms[rows] <- renames$New[R]
  }
  
  df <- distinct(df)
  
  # Fix cases where the databases have overlapping information but minor 
  # differences cause the rows to be different. 
  unique.epi <- distinct(df[,c("epitopes", "source_organisms")])
  epitopes <- unique(df$epitopes)
  
  for (R in 1:length(epitopes)) {
    rows <- which(unique.epi$epitopes == epitopes[R])
    
    # If multiple matches, concat the info together
    if (length(rows) > 1) {
      unique.epi$source_organisms[rows] <- paste(unique(unique.epi$source_organisms[rows]), collapse = " / ")
      print(paste(unique(unique.epi$source_organisms[rows]), collapse = " / "))
    } 
  }
  
  unique.epi <- distinct(unique.epi)
  
  for (R in 1:nrow(unique.epi)) {
    epi <- unique.epi$epitopes[R]
    rows <- which(df$epitopes == epi)
    df$source_organisms[rows] <- unique.epi$source_organisms[R]
  }
  
  df <- distinct(df)
  
  db <- append(db, list(df))
}

write.table(db[[1]], file_combined_db_alpha, sep="\t", row.names = FALSE, quote = FALSE)
write.table(db[[2]], file_combined_db_beta, sep="\t", row.names = FALSE, quote = FALSE)

rm(list=ls())



