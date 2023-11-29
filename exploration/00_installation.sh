#!/bin/bash

module purge all
module load R/4.3.0

Rscript --vanilla 00_installation.R
