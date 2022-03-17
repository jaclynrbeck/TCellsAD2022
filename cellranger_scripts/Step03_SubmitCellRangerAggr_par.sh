#!/bin/bash

CSV=/pub/jaclynb1/Parenchymal_SingleCell/scripts

sbatch ${CSV}/slurm_scripts/Step03_slurm_runCellRangerAggr.sh ParenchymalAggr ${CSV}/cellranger_aggr_par.csv
