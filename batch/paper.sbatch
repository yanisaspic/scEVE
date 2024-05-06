#!/bin/bash
#SBATCH --cpus-per-task 20
#SBATCH --mem 128GB

./setup.sh
module load r/4.3.1
Rscript ./paper.R
