#!/bin/bash -l
#SBATCH -A sens2017538
#SBATCH -n 16
#SBATCH -t 8:00:00
#SBATCH -J bmi_prs
#SBATCH --mail-type=ALL
#SBATCH --mail-user=marina.marchal@igp.uu.se


module purge # clear any inherited modules
module load R_packages/4.1.1
module load bioinfo-tools
module load plink/1.90b4.9




#Run Rscript in a clean R instance, output log file
Rscript --vanilla --verbose /proj/sens2017538/nobackup/marina/codes/bmi_pgs.R > slurm-${SLURM_JOBID}.Rout 2>&1

#Append log file to scripts log file
cat slurm-${SLURM_JOBID}.Rout >> slurm-${SLURM_JOBID}.out

#Remove R.out log
rm slurm-${SLURM_JOBID}.Rout
