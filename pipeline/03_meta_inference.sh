#!/bin/bash
#SBATCH -A p30791
#SBATCH -p normal
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=96G
#SBATCH --job-name=infermeta
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/p30791/methylation/out/03_infer_meta.out

module purge all 
module load R/4.3.0

cd ~/breast-methylation/pipeline/

Rscript --vanilla 03_meta_inference.R

module load python-miniconda3/4.12.0
source activate methylation

python 03_meta_inference_eval.py

