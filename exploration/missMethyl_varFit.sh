#!/bin/bash
#SBATCH -A p31931
#SBATCH -p normal
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=72G
#SBATCH --job-name=diffvar
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/p30791/methylation/out/missMethyl_varFit.out

module purge all 
module load R/4.3.0

cd ~/breast-methylation/exploration/

Rscript --vanilla missMethyl_varFit.R

