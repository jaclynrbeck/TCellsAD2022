# A list of filenames used in multiple places in the pipeline, to avoid mistakes
# or using the wrong file. 

# Main data directory
dir_data <- "data"

# Raw data from CellRanger
dir_cellranger <- file.path(dir_data, "CellRanger")
dir_filtered_counts <- file.path(dir_cellranger, "count", "filtered_feature_bc_matrix")
dir_vdj <- file.path(dir_cellranger, "vdj_t")

# Processed and raw clonotypes files
file_clonotypes_raw <- file.path(dir_vdj, "clonotypes.csv")
file_clonotypes_annotated <- file.path(dir_vdj, "filtered_contig_annotations.csv")
file_clonotypes_processed <- file.path(dir_cellranger, "clonotypes_processed.csv")

# Epitope matching
dir_epitope_matching <- file.path("data", "EpitopeMatching")
dir_epitope_output <- file.path("data", "EpitopeMatching", "output")
file_epitope_split_tra <- file.path(dir_epitope_matching, "Parenchymal_TCRalpha.txt")
file_epitope_split_trb <- file.path(dir_epitope_matching, "Parenchymal_TCRbeta.txt")
file_vdjdb_export <- file.path(dir_epitope_matching, "raw", "VJDB_SearchTable-2022-01-23 20_11_27.798.tsv")
file_vdjdb_clean_alpha <- file.path(dir_epitope_matching, "VDJdb_clean_TCRalpha.tsv")
file_vdjdb_clean_beta <- file.path(dir_epitope_matching, "VDJdb_clean_TCRbeta.tsv")
file_iedb_export <- file.path(dir_epitope_matching, "raw", "IEDB_tcell_receptor_table_export_1642966412.csv")
file_iedb_clean_alpha <- file.path(dir_epitope_matching, "IEDB_clean_TCRalpha.tsv")
file_iedb_clean_beta <- file.path(dir_epitope_matching, "IEDB_clean_TCRbeta.tsv")
file_combined_db <- file.path(dir_epitope_matching, "combined_epitope_db.tsv")
file_combined_alpha <- file.path(dir_epitope_output, "alpha_combined.tsv")
file_combined_beta <- file.path(dir_epitope_output, "beta_combined.tsv")
file_matched_epitopes <- file.path(dir_data, "matched_epitopes.tsv")

# Un-normalized and normalized Seurat objects, no analysis included
dir_seurat <- file.path(dir_data, "Seurat")
file_seurat_unnorm <- file.path(dir_seurat, "seurat_umi800_unnorm_2022-02-03.rds")
file_seurat_norm <- file.path(dir_seurat, "seurat_umi800_integrated_SCT_reg_2022-02-03.rds")
file_seurat_norm_cd8 <- file.path(dir_seurat, "seurat_umi800_integrated_SCT_reg_CD8_2022-02-14.rds")

# Analyzed Seurat objects
file_seurat_analyzed_allcells <- file.path(dir_seurat, "seurat_analyzed_allcells_2022-03-03.rds")
file_seurat_analyzed_cd8 <- file.path(dir_seurat, "seurat_analyzed_CD8_2022-02-14.rds")

# All cells combined files: Cluster / Genotype markers
dir_allcells <- file.path(dir_data, "AllCells")
dir_allcells_go <- file.path(dir_allcells, "GO_Analysis")
file_markers_combined_all <- file.path(dir_allcells, "all_markers_combined_2022-02-25.rds")
file_markers_combined_clusters <- file.path(dir_allcells, 'DiffGenes_byCluster_combined_2022-02-25.xlsx')
file_markers_combined_genotypes <- file.path(dir_allcells, 'DiffGenes_byGenotype_combined_2022-02-25.xlsx')
file_markers_combined_clusters_vs_genotype <- file.path(dir_allcells, 'DiffGenes_byGenotypePerCluster_combined_2022-02-25.xlsx')
file_clonotypes_combined_clusters <- file.path(dir_allcells, 'Clonotypes_byCluster_combined_2022-02-25.xlsx')
file_clonotypes_combined_genotypes <- file.path(dir_allcells, 'Clonotypes_byGenotype_combined_2022-02-25.xlsx')


# CD8 files: Cluster / Genotype markers
dir_cd8 <- file.path(dir_data, "CD8")
file_markers_cd8_all <- file.path(dir_cd8, "all_markers_cd8_2022-02-25.rds")
file_markers_cd8_clusters <- file.path(dir_cd8, 'DiffGenes_byCluster_CD8_2022-02-25.xlsx')
file_markers_cd8_genotypes <- file.path(dir_cd8, 'DiffGenes_byGenotype_CD8_2022-02-25.xlsx')
file_markers_cd8_clusters_vs_genotype <- file.path(dir_cd8, 'DiffGenes_byGenotypePerCluster_cd8_2022-02-25.xlsx')
file_clonotypes_cd8_clusters <- file.path(dir_cd8, 'Clonotypes_byCluster_CD8_2022-02-25.xlsx')
file_clonotypes_cd8_genotypes <- file.path(dir_cd8, 'Clonotypes_byGenotype_CD8_2022-02-25.xlsx')

# Figures
dir_figures <- file.path("figures")
dir_figures_allcells <- file.path(dir_figures, "AllCells")
dir_figures_cd8 <- file.path(dir_figures, "CD8")



