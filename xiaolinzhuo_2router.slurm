#!/bin/bash
#SBATCH -p stats 
#SBATCH -n 10
#SBATCH -t 01:30:00 
#SBATCH --mem-per-cpu=800M 
#SBATCH -o "dump/%A-%a.o"
#SBATCH -e "dump/%A-%a.e" 
#SBATCH --array=0-287:10

## LOAD SOFTWARE ENV ## 
module load R
input=xiaolinzhuo_2router.R

## EXECUTE CODE ##
srun time R CMD BATCH --no-save --no-restore --no-timing $input "dump/$input.$SLURM_ARRAY_TASK_ID.out" 
