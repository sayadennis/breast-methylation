#!/bin/bash
#SBATCH -A p30791
#SBATCH -p normal
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=200G
#SBATCH --job-name=preproc
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/p30791/methylation/out/01_preprocess_qc.out

module purge all 
module load R/4.3.0

cd ~/breast-methylation/pipeline/

Rscript --vanilla 01_preprocess_qc.R "/projects/p30791/methylation/data/IDAT_all" "/projects/p30791/methylation/sesame_out"

