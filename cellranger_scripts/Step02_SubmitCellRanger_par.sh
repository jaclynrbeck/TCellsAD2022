#!/bin/bash

CSV=/pub/jaclynb1/Parenchymal_SingleCell/scripts

samples=("WT-1" "WT-2" "5XFAD-1" "5XFAD-2" "PS19" "PS-5X")

for s in "${samples[@]}"
do
  sbatch ${CSV}/slurm_scripts/Step02_slurm_runCellRanger.sh ${s}_par ${CSV}/configs_parenchymal/cellranger_config_par_${s}.csv
done
