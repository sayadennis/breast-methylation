#!/bin/bash
#SBATCH -A p30791
#SBATCH -p short
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=16G
#SBATCH --job-name=pca
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/p30791/methylation/out/pca_betas_by_sample.out

module purge all 
module load python-miniconda3/4.12.0

source activate methylation

cd ~/breast-methylation/exploratory/

python pca_betas_by_sample.py

