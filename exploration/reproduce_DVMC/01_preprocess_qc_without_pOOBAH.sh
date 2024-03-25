#!/bin/bash
#SBATCH -A p31931
#SBATCH -p normal
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=140G
#SBATCH --job-name=preproc_nopoobah
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/p30791/methylation/out/01_preprocess_qc_without_pOOBAH.out

module purge all 
module load R/4.3.0

cd ~/breast-methylation/exploration/reproduce_DVMC/

Rscript --vanilla 01_preprocess_qc_without_pOOBAH.R "/projects/p30791/methylation/data/IDAT_all" "/projects/p30791/methylation/data/meta.csv" "/projects/p30791/methylation/sesame_out"

