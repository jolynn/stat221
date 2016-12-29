#!/bin/bash
#SBATCH -p stats 
#SBATCH -n 25
#SBATCH -t 01:00:00 
#SBATCH --mem-per-cpu=1G
#SBATCH -o "dump/%A-%a.o"
#SBATCH -e "dump/%A-%a.e" 
#SBATCH --array=2-287:1


## LOAD SOFTWARE ENV ## 
module load R
input=xiaolinzhuo_1router_all_times_informative.R

## EXECUTE CODE ##
srun time R CMD BATCH --no-save --no-restore --no-timing $input "dump/$input.$SLURM_ARRAY_TASK_ID.out" 