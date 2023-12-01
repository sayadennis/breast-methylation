#!/bin/bash
#SBATCH -A p30791
#SBATCH -p normal
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=48G
#SBATCH --job-name=plotqual
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/p30791/methylation/out/02_plot_quality.out

module purge all
module load R/4.3.0

cd ~/breast-methylation/pipeline/

Rscript --vanilla 02_plot_quality.R

