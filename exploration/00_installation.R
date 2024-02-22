if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("sesame")
install.packages("pals")
BiocManager::install("DNAcopy")
install.packages("ggrepel")
install.packages("gprofiler2")
BiocManager::install("GenoGAM")

library(devtools)
devtools::install_github("davidsjoberg/ggsankey")
