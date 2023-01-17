# TCellsAD2022

This repository contains the code used to process single cell sequencing data as described in (Sanchez, J, Beck, J, et al 2022) [manuscript in preparation].

Description of folders:

- `EpitopeMatching` - Helper functions for Steps 07 and 08 in the pipeline  
- `FigurePlotting` - Each R file in this folder will generate the corresponding 
figure from the paper. Some figures require post-editing in photoshop but most 
are usable as-is.
- `cellranger_scripts` - Shell scripts used to run CellRanger on UCI's High 
Performance Computing (HPC) cluster, using Slurm for job management. 
- `data` - Contains count and VDJ matrices output by CellRanger, the AD risk
gene list, IEDB and VDJdb exports, and the Reactome heat stress pathway genes. 
This folder is also where data from the pipeline is output by default. 
- `functions` - Helper functions for various steps in the pipeline.

---

## (Optional) Run CellRanger

This step has already been performed and the output is in this repository, but
in case CellRanger ever needs to be re-run on the raw data again, follow these
steps:

1. Download the fastq files from NCBI: [GSE212680](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?&acc=GSE212680). 
2. Download the mouse reference data for gene expression ([refdata-gex-mm10-2020-A.tar.gz](https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz)) and immune profiling ([refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0.tar.gz](https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0.tar.gz)) from CellRanger. 
3. Every csv file in `cellranger_scripts/configs_parenchymal` needs to be edited
so that file paths point to your file structure. 
4. Shell scripts will need to be tweaked if not using Slurm or if using a 
different way of configuring Slurm than UCI's HPC. 
5. Run the shell scripts in order from Step01 to Step03.

---

## Run Pipeline

Most of the filenames in `Filenames.R` should be set to good default values, and
file paths are relative to the main directory. Edit filenames as necessary.

Scripts are intended to be run in order from Step01 through Step09. See
comments at the header of individual files for more information about each 
script. 

Output from the scripts will be put in `data` under various sub-folders.

If run as-is, this pipeline will output exactly the data and figures used for
the paper, including supplementary tables. 



