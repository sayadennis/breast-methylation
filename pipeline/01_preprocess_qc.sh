#!/bin/bash
#SBATCH -A p31931
#SBATCH -p gengpu
#SBATCH --gres=gpu:a100:1
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=140G
#SBATCH --job-name=preproc
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/p30791/methylation/out/01_preprocess_qc.out

module purge all 
module load R/4.3.0

cd ~/breast-methylation/pipeline/

Rscript --vanilla 01_preprocess_qc.R "/projects/p30791/methylation/raw_data/IDAT_all" "/projects/p30791/methylation/raw_data/meta.csv" "/projects/p30791/methylation/sesame_data"

