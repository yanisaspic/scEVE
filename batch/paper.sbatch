#!/bin/bash
#SBATCH --cores-per-socket=8
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --array=0-13%5
#SBATCH --time=24:00:00

module load r/4.3.1

echo "#############################" 
echo "User:" $USER
echo "Date:" `date`
echo "Host:" `hostname`
echo "Directory:" `pwd`
echo "SLURM_JOBID:" $SLURM_JOBID
echo "SLURM_SUBMIT_DIR:" $SLURM_SUBMIT_DIR
echo "SLURM_JOB_NODELIST:" $SLURM_JOB_NODELIST
echo "#############################" 

REAL_DATASETS=("Baron_HumPan_1" \
"Baron_HumPan_2" \
"Baron_HumPan_3" \
"Baron_HumPan_4" \
"Baron_MouPan_1" \
"Baron_MouPan_2" \
"Darmanis_HumGBM" \
"JerbyArnon_HumMLM" \
"Lambrechts_HumNSCLC" \
"Li_HumCRC_a" \
"Li_HumCRC_b" \
"Peng_HumPDAC" \
"Tasic_MouBra" \
"VanGalen_HumAML")

Rscript ./benchmark.R ${DATASETS[$SLURM_ARRAY_TASK_ID]}
