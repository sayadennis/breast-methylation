#!/bin/bash
#SBATCH -A p30791
#SBATCH -p normal
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=128G
#SBATCH --job-name=modeling
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/p30791/methylation/out/05_differential_methylation.out

module purge all 
module load R/4.3.0

cd ~/breast-methylation/pipeline/

mkdir -p /projects/p30791/methylation/differential_methylation/
Rscript --vanilla 05_differential_methylation.R

