#!/bin/bash
#SBATCH -A p30791
#SBATCH -p normal
#SBATCH -t 2:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=24G
#SBATCH --job-name=gometh
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/p30791/methylation/out/missMethyl_gometh.out

module purge all 
module load R/4.3.0
module load hdf5/1.8.12-serial
module load openssl/1.1.1
module load libxml2/2.9.10-gcc-4.8.5

# cd /projects/p30791/methylation/
cd ~/breast-methylation/exploration/

Rscript --vanilla missMethyl_gometh.R

