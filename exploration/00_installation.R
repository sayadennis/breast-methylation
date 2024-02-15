if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("sesame")
install.packages("pals")
BiocManager::install("DNAcopy")
install.packages("ggrepel")
