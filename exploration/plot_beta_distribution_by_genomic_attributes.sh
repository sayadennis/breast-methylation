#!/bin/bash
#SBATCH -A p31931
#SBATCH -p normal
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=128G
#SBATCH --job-name=cgi_chrom
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/p30791/methylation/out/plot_beta_distribution_by_genomic_attributes.out

module purge all 
module load R/4.3.0

cd ~/breast-methylation/exploration/

Rscript --vanilla plot_beta_distribution_by_genomic_attributes.R

