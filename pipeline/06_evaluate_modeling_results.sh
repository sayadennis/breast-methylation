#!/bin/bash
#SBATCH -A p30791
#SBATCH -p short
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=4G
#SBATCH --job-name=dmleval
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/p30791/methylation/out/05_evaluate_modeling_results.out

module purge all 
module load python-miniconda3/4.12.0
source activate methylation

cd ~/breast-methylation/pipeline/

mkdir -p /projects/p30791/methylation/sesame_out/differential_methylation/
python 05_evaluate_modeling_results.py
