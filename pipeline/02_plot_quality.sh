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

cd ~/breast-methylation/pipeline/

module purge all
module load R/4.3.0

Rscript --vanilla 02_plot_quality.R "/projects/p30791/methylation/sesame_data" "/projects/p30791/methylation/plots" "/projects/p30791/methylation/raw_data/meta.csv"

module load python-miniconda3/4.12.0
source activate methylation

python 02_plot_quality.py "/projects/p30791/methylation/sesame_data" "/projects/p30791/methylation/raw_data/meta.csv" "/projects/p30791/methylation/plots"
