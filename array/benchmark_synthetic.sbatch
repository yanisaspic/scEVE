#!/bin/bash
#SBATCH --cores-per-socket=3
#SBATCH --cpus-per-task=3
#SBATCH --mem=64G
#SBATCH --array=0-299%10
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

IDS=({1..300})

Rscript ./benchmark_synthetic.R ${IDS[$SLURM_ARRAY_TASK_ID]}
