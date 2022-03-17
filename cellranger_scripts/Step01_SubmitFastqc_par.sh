SCRIPTS=/pub/jaclynb1/Parenchymal_SingleCell/scripts/slurm_scripts
FPATH=/pub/jaclynb1/Parenchymal_SingleCell/samples21042216
OUTPATH=/pub/jaclynb1/Parenchymal_SingleCell/Fastqc
FASTQC_DONE=${OUTPATH}/fastqc.done

for sample in ${FPATH}/*/*_R*; do
  sbatch ${SCRIPTS}/Step01_slurm_runFastqc.sh $OUTPATH $sample
  echo -e $sample >> $FASTQC_DONE
done
