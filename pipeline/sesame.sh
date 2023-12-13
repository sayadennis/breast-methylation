#!/bin/bash
#SBATCH -A p30791
#SBATCH -p normal
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=1G
#SBATCH --job-name=nfsesame
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/p30791/methylation/out/sesame_nextflow.out

module purge all
module load nextflow/23.04.3

cd ~/breast-methylation/pipeline/

nextflow run -c ./sesame.nextflow.config -cache false -work-dir "/projects/p30791/methylation/nf_workdir" sesame.nf
