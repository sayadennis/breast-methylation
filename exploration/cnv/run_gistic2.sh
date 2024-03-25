#!/bin/bash
#SBATCH -A p31931
#SBATCH -p normal
#SBATCH -t 24:00:00
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="methgistic"
#SBATCH --output=/projects/p30791/methylation/out/run_gistic2.out

module purge all
module load python-miniconda3/4.12.0
source activate methylation

#####################
#### Run GISTIC2 ####
#####################

din="/projects/p30791/methylation/sesame_out/copy_number"
refgenefile="./refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat"

cd /projects/p30791/GISTIC2/

for tissue_category in "Normal" "CUB" "OQ" "AN" "TU"
do
    dout="/projects/p30791/methylation/sesame_out/copy_number/gistic2_out_$tissue_category"
    mkdir -p $dout
    ./gistic2 \
        -b $dout \
        -seg $din/segs_${tissue_category}.tsv \
        -refgene $refgenefile \
        -genegistic 0 \
        -smallmem 0 \
        -brlen 0.98 \
        -conf 0.90 \
        -armpeel 0 \
        -savegene 1 \
        -gcm extreme
done

#########################
#### Process outputs ####
#########################

# python ~/bbcar/repo/01_processing/input/cnv/09_gistic_to_model_features.py

