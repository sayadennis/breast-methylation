#!/bin/bash
#SBATCH -A p31931
#SBATCH -p normal
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=48G
#SBATCH --job-name=cnv
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/p30791/methylation/out/copy_number.out

module purge all 
module load R/4.3.0

cd ~/breast-methylation/exploration/

Rscript --vanilla copy_number.R "/projects/p30791/methylation/sesame_out" "/projects/p30791/methylation/plots/copy_number" "/projects/p30791/methylation/data/meta.csv"

