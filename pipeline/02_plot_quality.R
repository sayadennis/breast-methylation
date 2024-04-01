library(stringr)
library(sesame)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Please provide 3 arguments: ",
    "1) input directory, ",
    "2) output directory to write plots to, and ",
    "3) metadata file path.",
    call. = FALSE
  )
} else {
  din <- args[1] # e.g. "/projects/p30791/methylation/sesame_data"
  dout <- args[2] # e.g. "/projects/p30791/methylation/plots"
  meta_fn <- args[3] # e.g. "/projects/p30791/methylation/raw_data/meta.csv"
}

## Load data
sdfs <- readRDS(paste0(din, "/sdf_raw.RDS"))
qcs <- readRDS(paste0(din, "/qcs_raw.RDS"))
qc_df <- read.csv(paste0(din, "/qc_raw.csv"))
meta <- read.csv(meta_fn)

# Note samples with high distortion
rm_sampleids <- readLines(paste0(din, "/exclude_IDATs.txt"))
print(paste0(
  "Samples to be removed from downstream analysis: ",
  paste0(rm_sampleids, collapse = ", ")
))
print("Meta data of these samples:")
print(meta[meta$IDAT %in% rm_sampleids, ])

## Plot QC scatter plots
plot_sampleids <- c(rm_sampleids, meta$IDAT[1])
for (idat_id in plot_sampleids) {
  sdf <- sdfs[idat_id][[1]]

  png(
    filename = paste0(dout, "/qc_plotRedGrnQQ_", idat_id, ".png"),
    type = "cairo"
  )
  sesameQC_plotRedGrnQQ(sdf, main = idat_id)
  dev.off()

  png(
    filename = paste0(dout, "/qc_sesameQC_plotIntensVsBetas_", idat_id, ".png"),
    type = "cairo"
  )
  sesameQC_plotIntensVsBetas(sdf)
  dev.off()
}

png(filename = paste0(dout, "/qc_sesameQC_plotHeatSNPs.png"), type = "cairo")
sesameQC_plotHeatSNPs(sdfs) # plots SNP probes to detect sample swaps
dev.off()
