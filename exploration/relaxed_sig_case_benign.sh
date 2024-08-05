#!/bin/bash
#SBATCH -A p31931
#SBATCH -p normal
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=24G
#SBATCH --job-name=relsig_kycg
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/p30791/methylation/out/relaxed_sig_case_benign.out

module purge all 
module load R/4.3.0

cd ~/breast-methylation/exploration/

Rscript --vanilla relaxed_sig_case_benign.R

