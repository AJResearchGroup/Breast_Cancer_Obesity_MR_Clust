#!/bin/bash -l
#SBATCH -A sens2017538
#SBATCH -n 16
#SBATCH -t 4:00:00
#SBATCH -J pgs_lm
#SBATCH --mail-type=ALL
#SBATCH --mail-user=marina.marchal@igp.uu.se


module purge # clear any inherited modules
module load R_packages/4.3.1

#Run Rscript in a clean R instance, output log file
Rscript --vanilla --verbose /proj/sens2017538/nobackup/marina/codes/pgs_lm.R > slurm-${SLURM_JOBID}.Rout 2>&1

#Append log file to scripts log file
cat slurm-${SLURM_JOBID}.Rout >> slurm-${SLURM_JOBID}.out

#Remove R.out log
rm slurm-${SLURM_JOBID}.Rout