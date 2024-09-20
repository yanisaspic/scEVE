#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=24:00:00

module load r/4.3.1

dataset=$1

echo "#############################" 
echo "User:" $USER
echo "Date:" `date`
echo "Host:" `hostname`
echo "Directory:" `pwd`
echo "SLURM_JOBID:" $SLURM_JOBID
echo "SLURM_SUBMIT_DIR:" $SLURM_SUBMIT_DIR
echo "SLURM_JOB_NODELIST:" $SLURM_JOB_NODELIST
echo "DATASET:" $dataset
echo "#############################" 

Rscript ./run/get_real.R $dataset

mv ./$dataset ./results/figures/$dataset
mv ./$dataset.xlsx ./results/records/$dataset.xlsx
