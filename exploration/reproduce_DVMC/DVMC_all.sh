#!/bin/bash
#SBATCH -A p31931
#SBATCH -p normal
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=96G
#SBATCH --job-name=dvmc_all
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/p30791/methylation/out/DVMC_all.out

module purge all 
module load R/4.3.0

cd ~/breast-methylation/exploration/reproduce_DVMC/

Rscript --vanilla DVMC_all.R

