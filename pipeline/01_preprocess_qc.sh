#!/bin/bash
#SBATCH -A p31931
#SBATCH -p normal
#SBATCH -t 8:00:00
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=140G
#SBATCH --job-name=preproc
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/p30791/methylation/out/01_preprocess_qc.out

module purge all 
module load R/4.3.0

cd ~/breast-methylation/pipeline/

Rscript --vanilla 01_preprocess_qc.R "/projects/p30791/methylation/data/IDAT_all" "/projects/p30791/methylation/data/meta.csv" "/projects/p30791/methylation/sesame_out"

