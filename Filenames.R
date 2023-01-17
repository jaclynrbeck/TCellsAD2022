# A list of filenames used in multiple places in the pipeline, to avoid mistakes
# or using the wrong file. 

# Main data directory
dir_data <- "data"


# Raw data from CellRanger
dir_cellranger <- file.path(dir_data, "CellRanger")
dir_filtered_counts <- file.path(dir_cellranger, "count", "filtered_feature_bc_matrix")
dir_vdj <- file.path(dir_cellranger, "vdj_t")


# External data folder
dir_data_external <- file.path(dir_data, "External")


# Processed and raw clonotypes files
file_clonotypes_summary <- file.path(dir_vdj, "clonotypes.csv")
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
file_vdjdb_export <- file.path(dir_data_external, "VDJDBSearchTable-2023-01-17 03_51_57.099.tsv")
file_vdjdb_clean_alpha <- file.path(dir_epitope_matching, "VDJdb_clean_TCRalpha.tsv")
file_vdjdb_clean_beta <- file.path(dir_epitope_matching, "VDJdb_clean_TCRbeta.tsv")
file_iedb_export <- file.path(dir_data_external, "IEDB_tcell_receptor_table_export_1673928207.csv")
file_iedb_clean_alpha <- file.path(dir_epitope_matching, "IEDB_clean_TCRalpha.tsv")
file_iedb_clean_beta <- file.path(dir_epitope_matching, "IEDB_clean_TCRbeta.tsv")
file_combined_db_alpha <- file.path(dir_epitope_matching, "combined_epitope_db_alpha.tsv")
file_combined_db_beta <- file.path(dir_epitope_matching, "combined_epitope_db_beta.tsv")
file_combined_alpha <- file.path(dir_epitope_output, "alpha_combined.tsv")
file_combined_beta <- file.path(dir_epitope_output, "beta_combined.tsv")
file_matched_epitopes <- file.path(dir_epitope_output, "matched_epitopes.tsv")


# Un-normalized and normalized Seurat objects, no analysis included
dir_seurat <- file.path(dir_data, "Seurat")
file_seurat_unnorm <- file.path(dir_seurat, "seurat_unnorm_2022-09-02.rds")
file_seurat_norm <- file.path(dir_seurat, "seurat_integrated_2022-09-02.rds")
file_seurat_norm_cd8 <- file.path(dir_seurat, "seurat_integrated_CD8_2022-09-02.rds")
file_seurat_norm_cd4 <- file.path(dir_seurat, "seurat_SCT_CD4_2022-09-02.rds")
file_seurat_norm_gd <- file.path(dir_seurat, "seurat_SCT_GD_2022-09-02.rds")


# Analyzed Seurat objects
file_seurat_analyzed_allcells <- file.path(dir_seurat, "seurat_analyzed_allcells_2022-09-02.rds")
file_seurat_analyzed_cd8 <- file.path(dir_seurat, "seurat_analyzed_CD8_2022-09-08.rds")
file_seurat_analyzed_cd4 <- file.path(dir_seurat, "seurat_analyzed_CD4_2022-09-02.rds")
file_seurat_analyzed_gd <- file.path(dir_seurat, "seurat_analyzed_GD_2022-09-02.rds")


# All cells combined files: Cluster / Genotype markers
dir_allcells <- file.path(dir_data, "AllCells")
dir_allcells_go <- file.path(dir_allcells, "GO_Analysis")
dir_allcells_dpa <- file.path(dir_allcells, "DPA")
file_markers_combined_all <- file.path(dir_allcells, "all_markers_combined_2022-09-02.rds")
file_markers_combined_clusters <- file.path(dir_allcells, 'DiffGenes_byCluster_combined_2022-09-02.xlsx')
file_markers_combined_genotypes <- file.path(dir_allcells, 'DiffGenes_byGenotype_combined_2022-09-02.xlsx')
file_markers_combined_clusters_vs_genotype <- file.path(dir_allcells, 'DiffGenes_byGenotypePerCluster_combined_2022-09-02.xlsx')
file_clonotypes_combined_clusters <- file.path(dir_allcells, 'Clonotypes_byCluster_combined_2022-09-02.xlsx')
file_clonotypes_combined_genotypes <- file.path(dir_allcells, 'Clonotypes_byGenotype_combined_2022-09-02.xlsx')
file_ad_risk_combined <- file.path(dir_allcells, "AD_DiffGenes_byCluster_allCells_2022-09-02.xlsx")
file_dpa_allcells_glm_summary <- file.path(dir_allcells_dpa, "DPA_combined_GLM_summary_2022-09-09.csv")
file_dpa_allcells_pairwise <- file.path(dir_allcells_dpa, "DPA_combined_pairwise_comparisons_2022-09-09.xlsx")


# CD8 files: Cluster / Genotype markers
dir_cd8 <- file.path(dir_data, "CD8")
dir_cd8_go <- file.path(dir_cd8, "GO_Analysis")
dir_cd8_dpa <- file.path(dir_cd8, "DPA")
file_markers_cd8_all <- file.path(dir_cd8, "all_markers_cd8_2022-09-08.rds")
file_markers_cd8_clusters <- file.path(dir_cd8, 'DiffGenes_byCluster_CD8_2022-09-08.xlsx')
file_markers_cd8_genotypes <- file.path(dir_cd8, 'DiffGenes_byGenotype_CD8_2022-09-08.xlsx')
file_markers_cd8_clusters_vs_genotype <- file.path(dir_cd8, 'DiffGenes_byGenotypePerCluster_cd8_2022-09-08.xlsx')
file_clonotypes_cd8_clusters <- file.path(dir_cd8, 'Clonotypes_byCluster_CD8_2022-09-08.xlsx')
file_clonotypes_cd8_genotypes <- file.path(dir_cd8, 'Clonotypes_byGenotype_CD8_2022-09-08.xlsx')
file_dpa_cd8_glm_summary <- file.path(dir_cd8_dpa, "DPA_CD8_GLM_summary_2022-09-09.csv")
file_dpa_cd8_pairwise <- file.path(dir_cd8_dpa, "DPA_CD8_pairwise_comparisons_2022-09-09.xlsx")


# CD4 files: Cluster / Genotype markers
dir_cd4 <- file.path(dir_data, "CD4")
file_markers_cd4_all <- file.path(dir_cd4, "all_markers_cd4_2022-09-08.rds")
file_markers_cd4_clusters <- file.path(dir_cd4, 'DiffGenes_byCluster_CD4_2022-09-08.xlsx')
file_markers_cd4_genotypes <- file.path(dir_cd4, 'DiffGenes_byGenotype_CD4_2022-09-08.xlsx')
file_markers_cd4_clusters_vs_genotype <- file.path(dir_cd4, 'DiffGenes_byGenotypePerCluster_cd4_2022-09-08.xlsx')
file_clonotypes_cd4_clusters <- file.path(dir_cd4, 'Clonotypes_byCluster_CD4_2022-09-08.xlsx')
file_clonotypes_cd4_genotypes <- file.path(dir_cd4, 'Clonotypes_byGenotype_CD4_2022-09-08.xlsx')


# Gamma delta files: Cluster / Genotype markers
dir_gd <- file.path(dir_data, "GammaDelta")
file_markers_gd_all <- file.path(dir_gd, "all_markers_gd_2022-09-08.rds")
file_markers_gd_clusters <- file.path(dir_gd, 'DiffGenes_byCluster_GD_2022-09-08.xlsx')
file_markers_gd_genotypes <- file.path(dir_gd, 'DiffGenes_byGenotype_GD_2022-09-08.xlsx')
file_markers_gd_clusters_vs_genotype <- file.path(dir_gd, 'DiffGenes_byGenotypePerCluster_GD_2022-09-08.xlsx')
file_clonotypes_gd_clusters <- file.path(dir_gd, 'Clonotypes_byCluster_GD_2022-09-08.xlsx')
file_clonotypes_gd_genotypes <- file.path(dir_gd, 'Clonotypes_byGenotype_GD_2022-09-08.xlsx')


# Figures
dir_figures <- file.path("figures")
dir_figures_allcells <- file.path(dir_figures, "AllCells")
dir_figures_cd8 <- file.path(dir_figures, "CD8")
dir_figures_cd4 <- file.path(dir_figures, "CD4")
dir_figures_gd <- file.path(dir_figures, "GammaDelta")
dir_figures_paper <- file.path(dir_figures, "Paper")


# Ensure all directories we plan to read/write to exist. 
dirs <- grep("dir_", ls(), value = TRUE)
for (D in dirs) {
  dir_eval <- eval(parse(text = D))
  if (!file.exists(dir_eval)) {
    dir.create(dir_eval)
    print(paste0("Created directory ", dir_eval))
  }
}

rm(D)

