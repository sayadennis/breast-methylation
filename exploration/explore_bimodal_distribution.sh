#!/bin/bash
#SBATCH -A p30791
#SBATCH -p normal
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=48G
#SBATCH --job-name=bimodal
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/p30791/methylation/out/explore_bimodal_distribution.out

module purge all 
module load python-miniconda3/4.12.0
source activate methylation

cd ~/breast-methylation/exploration/

python explore_bimodal_distribution.py

