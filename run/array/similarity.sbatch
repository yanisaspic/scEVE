#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --array=0-149%10
#SBATCH --time=24:00:00

module load r/4.3.1

IDS=({1..150})
id=${IDS[$SLURM_ARRAY_TASK_ID]}

echo "#############################" 
echo "User:" $USER
echo "Date:" `date`
echo "Host:" `hostname`
echo "Directory:" `pwd`
echo "SLURM_JOBID:" $SLURM_JOBID
echo "SLURM_SUBMIT_DIR:" $SLURM_SUBMIT_DIR
echo "SLURM_JOB_NODELIST:" $SLURM_JOB_NODELIST
echo "DATASET ID:" $id
echo "#############################" 

Rscript ./run/get_similarity.R $id
