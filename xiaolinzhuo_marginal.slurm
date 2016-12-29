#!/bin/bash
#SBATCH -p stats 
#SBATCH -n 1
#SBATCH -t 01:00:00 
#SBATCH --mem-per-cpu=1000M
#SBATCH -o "%A.o"
#SBATCH -e "%A.e" 

## LOAD SOFTWARE ENV ## 
module load R
input=xiaolinzhuo_marginal.R

## EXECUTE CODE ##
R CMD BATCH --no-save --no-restore --no-timing $input $input.$SLURM_ARRAY_TASK_ID.out
