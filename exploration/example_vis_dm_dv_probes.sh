#!/bin/bash
#SBATCH -A p31931
#SBATCH -p normal
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=72G
#SBATCH --job-name=dummyvis
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/p30791/methylation/out/example_vis_dm_dv_probes.out

module purge all 
module load python-miniconda3/4.12.0
source activate methylation

cd ~/breast-methylation/exploration/

python example_vis_dm_dv_probes.py

