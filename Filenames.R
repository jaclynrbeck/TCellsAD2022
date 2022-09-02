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
file_clonotypes_raw_edited <- file.path(dir_vdj, "clonotypes_edited.csv")
file_clonotypes_annotated <- file.path(dir_vdj, "filtered_contig_annotations.csv")
file_clonotypes_annotated_edited <- file.path(dir_vdj, "filtered_contig_annotations_edited.csv")
file_clonotypes_processed <- file.path(dir_cellranger, "clonotypes_processed.csv")

# Gene symbols / ensembl IDs files
file_gene_symbols <- file.path(dir_data, "gene_symbols.rds")
file_excluded_genes <- file.path(dir_data, "excluded_genes.rds")

# Epitope matching
dir_epitope_matching <- file.path(dir_data, "EpitopeMatching")
dir_epitope_output <- file.path(dir_epitope_matching, "output")
file_epitope_split_tra <- file.path(dir_epitope_matching, "Parenchymal_TCRalpha.txt")
file_epitope_split_trb <- file.path(dir_epitope_matching, "Parenchymal_TCRbeta.txt")
file_vdjdb_export <- file.path(dir_epitope_matching, "raw", "VJDB_SearchTable-2022-01-23 20_11_27.798.tsv")
file_vdjdb_clean_alpha <- file.path(dir_epitope_matching, "VDJdb_clean_TCRalpha.tsv")
file_vdjdb_clean_beta <- file.path(dir_epitope_matching, "VDJdb_clean_TCRbeta.tsv")
file_iedb_export <- file.path(dir_epitope_matching, "raw", "IEDB_tcell_receptor_table_export_1642966412.csv")
file_iedb_clean_alpha <- file.path(dir_epitope_matching, "IEDB_clean_TCRalpha.tsv")
file_iedb_clean_beta <- file.path(dir_epitope_matching, "IEDB_clean_TCRbeta.tsv")
file_combined_db_alpha <- file.path(dir_epitope_matching, "combined_epitope_db_alpha.tsv")
file_combined_db_beta <- file.path(dir_epitope_matching, "combined_epitope_db_beta.tsv")
file_combined_alpha <- file.path(dir_epitope_output, "alpha_combined.tsv")
file_combined_beta <- file.path(dir_epitope_output, "beta_combined.tsv")
file_matched_epitopes <- file.path(dir_data, "matched_epitopes.tsv")

# VDJTools
dir_vdjtools <- file.path(dir_data, "VDJTools")
dir_vdjtools_input <- file.path(dir_vdjtools, "input")
dir_vdjtools_output <- file.path(dir_vdjtools, "output")

# Un-normalized and normalized Seurat objects, no analysis included
dir_seurat <- file.path(dir_data, "Seurat")
file_seurat_unnorm <- file.path(dir_seurat, "seurat_unnorm_2022-04-12.rds")
file_seurat_norm <- file.path(dir_seurat, "seurat_integrated_SCT_2022-04-12.rds")
file_seurat_norm_cd8 <- file.path(dir_seurat, "seurat_integrated_SCT_CD8_2022-05-10.rds")
file_seurat_norm_cd4 <- file.path(dir_seurat, "seurat_integrated_SCT_CD4_2022-05-26.rds")
file_seurat_norm_gd <- file.path(dir_seurat, "seurat_SCT_GD_2022-06-28.rds")

# Analyzed Seurat objects
file_seurat_analyzed_allcells <- file.path(dir_seurat, "seurat_analyzed_allcells_2022-04-14.rds")
file_seurat_analyzed_cd8 <- file.path(dir_seurat, "seurat_analyzed_CD8_2022-05-10.rds")
file_seurat_analyzed_cd4 <- file.path(dir_seurat, "seurat_analyzed_CD4_2022-05-26.rds")
file_seurat_analyzed_gd <- file.path(dir_seurat, "seurat_analyzed_GD_2022-06-28.rds")

# All cells combined files: Cluster / Genotype markers
dir_allcells <- file.path(dir_data, "AllCells")
dir_allcells_go <- file.path(dir_allcells, "GO_Analysis")
file_markers_combined_all <- file.path(dir_allcells, "all_markers_combined_2022-04-14.rds")
file_markers_combined_clusters <- file.path(dir_allcells, 'DiffGenes_byCluster_combined_2022-04-14.xlsx')
file_markers_combined_genotypes <- file.path(dir_allcells, 'DiffGenes_byGenotype_combined_2022-04-14.xlsx')
file_markers_combined_clusters_vs_genotype <- file.path(dir_allcells, 'DiffGenes_byGenotypePerCluster_combined_2022-04-14.xlsx')
file_clonotypes_combined_clusters <- file.path(dir_allcells, 'Clonotypes_byCluster_combined_2022-04-14.xlsx')
file_clonotypes_combined_genotypes <- file.path(dir_allcells, 'Clonotypes_byGenotype_combined_2022-04-14.xlsx')


# CD8 files: Cluster / Genotype markers
dir_cd8 <- file.path(dir_data, "CD8")
dir_cd8_go <- file.path(dir_cd8, "GO_Analysis")
file_markers_cd8_all <- file.path(dir_cd8, "all_markers_cd8_2022-05-19.rds")
file_markers_cd8_clusters <- file.path(dir_cd8, 'DiffGenes_byCluster_CD8_2022-05-19.xlsx')
file_markers_cd8_genotypes <- file.path(dir_cd8, 'DiffGenes_byGenotype_CD8_2022-05-19.xlsx')
file_markers_cd8_clusters_vs_genotype <- file.path(dir_cd8, 'DiffGenes_byGenotypePerCluster_cd8_2022-05-19.xlsx')
file_clonotypes_cd8_clusters <- file.path(dir_cd8, 'Clonotypes_byCluster_CD8_2022-05-19.xlsx')
file_clonotypes_cd8_genotypes <- file.path(dir_cd8, 'Clonotypes_byGenotype_CD8_2022-05-19.xlsx')


# CD4 files: Cluster / Genotype markers
dir_cd4 <- file.path(dir_data, "CD4")
dir_cd4_go <- file.path(dir_cd4, "GO_Analysis")
file_markers_cd4_all <- file.path(dir_cd4, "all_markers_cd4_2022-05-26.rds")
file_markers_cd4_clusters <- file.path(dir_cd4, 'DiffGenes_byCluster_CD4_2022-05-26.xlsx')
file_markers_cd4_genotypes <- file.path(dir_cd4, 'DiffGenes_byGenotype_CD4_2022-05-26.xlsx')
file_markers_cd4_clusters_vs_genotype <- file.path(dir_cd4, 'DiffGenes_byGenotypePerCluster_cd4_2022-05-26.xlsx')
file_clonotypes_cd4_clusters <- file.path(dir_cd4, 'Clonotypes_byCluster_CD4_2022-05-26.xlsx')
file_clonotypes_cd4_genotypes <- file.path(dir_cd4, 'Clonotypes_byGenotype_CD4_2022-05-26.xlsx')

# Gamma delta files: Cluster / Genotype markers
dir_gd <- file.path(dir_data, "GammaDelta")
file_markers_gd_all <- file.path(dir_gd, "all_markers_gd_2022-06-28.rds")
file_markers_gd_clusters <- file.path(dir_gd, 'DiffGenes_byCluster_GD_2022-06-28.xlsx')
file_markers_gd_genotypes <- file.path(dir_gd, 'DiffGenes_byGenotype_GD_2022-06-28.xlsx')
file_markers_gd_clusters_vs_genotype <- file.path(dir_gd, 'DiffGenes_byGenotypePerCluster_GD_2022-06-28.xlsx')
file_clonotypes_gd_clusters <- file.path(dir_gd, 'Clonotypes_byCluster_GD_2022-06-28.xlsx')
file_clonotypes_gd_genotypes <- file.path(dir_gd, 'Clonotypes_byGenotype_GD_2022-06-28.xlsx')

# Figures
dir_figures <- file.path("figures")
dir_figures_allcells <- file.path(dir_figures, "AllCells")
dir_figures_cd8 <- file.path(dir_figures, "CD8")
dir_figures_cd4 <- file.path(dir_figures, "CD4")
dir_figures_gd <- file.path(dir_figures, "GammaDelta")

