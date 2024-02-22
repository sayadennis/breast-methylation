library(stringr)
library(sesame)
library(ggplot2)

source("~/breast-methylation/exploration/cnv_custom.R")

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Please provide 3 arguments: ",
    "1) input directory, ",
    "2) output directory to write plots to, and ",
    "3) metadata file path.",
    call. = FALSE
  )
} else {
  din <- args[1] # e.g. "/projects/p30791/methylation/sesame_out"
  dout <- args[2] # e.g. "/projects/p30791/methylation/plots/copy_number"
  meta_fn <- args[3] # e.g. "/projects/p30791/methylation/data/meta.csv"
}

if (!file.exists(dout)) {
  dir.create(dout)
  cat("Directory created:", dout, "\n")
}

meta <- read.csv(meta_fn)

## Load data
sdfs <- readRDS(paste0(din, "/sdf_processed.RDS"))

for (i in seq_along(sdfs)) {
  segs <- cnSegmentation(sdfs[[i]])
  idat_id <- names(sdfs)[[i]]
  tissue_category <- meta[meta$IDAT == idat_id, "Sample.Region"]
  cat("#### ", idat_id, "(", tissue_category, ") ####\n")
  # TODO: Explore removing Y chromosome from this plot
  gg_plot <- visualizeSegments(segs, to.plot = paste0("chr", c(1:22))) +
    ggtitle(paste0("CN Segments: ", idat_id, "(", tissue_category, ")")) +
    theme(legend.position = "none") +
    ylim(-1.5, 1.1)
  ggsave(
    filename = paste0(
      dout, "/cnv_segments_", tissue_category, "_", idat_id, ".png"
    ),
    plot = gg_plot
  )
}
