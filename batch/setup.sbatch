#!/bin/bash
#SBATCH --cpus-per-task 2
#SBATCH --mem 10GB
module load r/4.3.1
module load python/3.9
./setup.sh