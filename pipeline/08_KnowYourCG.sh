#!/bin/bash
#SBATCH -A p31931
#SBATCH -p normal
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=24G
#SBATCH --job-name=kycg
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/p30791/methylation/out/KnowYourCG.out

module purge all 
module load R/4.3.0

cd ~/breast-methylation/pipeline/

Rscript --vanilla 08_KnowYourCG.R "probesets.json" "/projects/p30791/methylation/differential_methylation/KYCG" "/projects/p30791/methylation/plots/differential_methylation"

