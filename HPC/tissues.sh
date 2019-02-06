#!/bin/bash
#SBATCH --time=1440
#SBATCH --mem=30000
#SBATCH --ntasks=4
#SBATCH	--cpus-per-task=4
#SBATCH --nodes=2
#SBATCH --output=output_wasp6PI_%j.txt
#SBATCH --error=error_out_wasp6pi_%j.txt
#SBATCH --job-name=WASP6
#SBATCH --partition=BIOINF_Std
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fernando.buenogutierrez@wur.nl
#SBATCH --array=1-1887


#script for something my sc
module load R/3.4.0

Rscript tissues_1.R ${SLURM_ARRAY_TASK_ID} &
Rscript tissues_2.R ${SLURM_ARRAY_TASK_ID} &
Rscript tissues_3.R ${SLURM_ARRAY_TASK_ID} &
Rscript tissues_4.R ${SLURM_ARRAY_TASK_ID} 
