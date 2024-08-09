#!/bin/bash
#SBATCH -A p30791
#SBATCH -p normal
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=96G
#SBATCH --job-name=diffvar
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/p30791/methylation/out/06_differential_variability.out

module purge all 
# module load singularity/3.8.1
module load R/4.3.0
module load hdf5/1.8.12-serial
module load openssl/1.1.1
module load libxml2/2.9.10-gcc-4.8.5
module load blas-lapack/3.12.0-gcc-11.2.0

# cd /projects/p30791/methylation/
cd ~/breast-methylation/pipeline/

mkdir -p /projects/p30791/methylation/differential_variability/
# singularity exec /projects/p30791/r_missmethyl_4.3.1.sif \
#     Rscript --vanilla ~/breast-methylation/pipeline/06_differential_variability.R \
#     -B /projects:/projects 

Rscript --vanilla 06_differential_variability.R

