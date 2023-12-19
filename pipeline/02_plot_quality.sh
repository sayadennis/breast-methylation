#!/bin/bash
#SBATCH -A p30791
#SBATCH -p short
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=72G
#SBATCH --job-name=plotqual
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/p30791/methylation/out/02_plot_quality.out

module purge all
module load R/4.3.0

cd ~/breast-methylation/pipeline/

Rscript --vanilla 02_plot_quality.R "/projects/p30791/methylation/sesame_out" "/projects/p30791/methylation/plots" "/projects/p30791/methylation/data/meta.csv"

