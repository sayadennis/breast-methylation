library(stringr)
library(sesame)
library(ggplot2)
library(gridExtra)

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

## Load data and metadata
sdfs <- readRDS(paste0(din, "/sdf_processed.RDS"))
meta <- read.csv(meta_fn)

meta <- meta[meta$IDAT %in% names(sdfs), ]
sdfs <- sdfs[meta$IDAT]

############################################################
#### Create CN Segment tables and plots for each sample ####
############################################################

tissue_categories <- c("Normal", "CUB", "OQ", "AN", "TU")
seg.cols <- c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean", "seg.sd", "seg.median")
seg.dfs <- list(
  Normal = data.frame(matrix(ncol = length(seg.cols), nrow = 0)),
  CUB = data.frame(matrix(ncol = length(seg.cols), nrow = 0)),
  OQ = data.frame(matrix(ncol = length(seg.cols), nrow = 0)),
  AN = data.frame(matrix(ncol = length(seg.cols), nrow = 0)),
  TU = data.frame(matrix(ncol = length(seg.cols), nrow = 0))
)
for (tissue_category in tissue_categories) {
  names(seg.dfs[[tissue_category]]) <- seg.cols
}

for (i in seq_along(sdfs)) {
  segs <- cnSegmentation(sdfs[[i]])
  idat_id <- names(sdfs)[[i]]
  tissue_category <- meta[meta$IDAT == idat_id, "Sample.Region"]
  cat("#### ", idat_id, "(", tissue_category, ") ####\n")
  segs$seg.signals$ID <- idat_id
  seg.dfs[[tissue_category]] <- rbind(seg.dfs[[tissue_category]], segs$seg.signals)
  patient_id <- meta[meta$IDAT == idat_id, "ID"]
  gg_plot <- visualizeSegments(segs, to.plot = paste0("chr", c(1:22))) +
    ggtitle(paste0("CN Segments: ", idat_id, "(", tissue_category, ")")) +
    theme(legend.position = "none") +
    ylim(-1, 1)
  ggsave(
    filename = paste0(
      dout, "/cnv_segments_", tissue_category, "_", idat_id, ".png"
    ),
    plot = gg_plot
  )
}

if (!file.exists(paste0(din, "/copy_number"))) {
  dir.create(paste0(din, "/copy_number"))
}

for (tissue_category in names(seg.dfs)) {
  write.table(
    seg.dfs[[tissue_category]][, c(1:6)], # only the first 6 columns
    file = paste0(din, "/copy_number/segs_", tissue_category, ".csv"),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}

#########################################
#### Combine Case CNVs into one plot ####
#########################################

tissue_categories <- c("CUB", "OQ", "AN", "TU")

num_patients <- 4

aspect_ratio <- 4 # ratio for each panel
plot_iter <- 1 # plot iteration number for layout_matrix argument
nonnamed_plotlist <- c() # non-named list to store plots for grid.arrange

widths <- rep(aspect_ratio, 4)
heights <- rep(1, num_patients)
layout_matrix <- matrix(nrow = num_patients, ncol = 4)
case_ids <- unique(meta[meta$Case.Control == "case", "ID"])

# Write the case ID list so that we can reference it later
writeLines(case_ids, paste0(dout, "/case_id_cn_plot_order.txt"))

for (i in seq_along(case_ids[1:num_patients])) {
  # Get patient ID and all their metadata
  patient_id <- case_ids[i]
  patient_meta <- meta[meta$ID == patient_id, ]
  # Iterate over tissue categories
  for (j in seq_along(tissue_categories)) {
    tissue_category <- tissue_categories[j]
    # If this patient has this tissue category
    if (tissue_category %in% patient_meta$Sample.Region) {
      # Get IDAT ID of this sample
      idat_id <- patient_meta[patient_meta$Sample.Region == tissue_category, "IDAT"]
      # Get CN segments and plot
      segs <- cnSegmentation(sdfs[[idat_id]])
      p <- visualizeSegments(segs, to.plot = paste0("chr", c(1:22))) +
        ggtitle(NULL) +
        theme(
          legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
        ) +
        ylim(-1, 1)
      # Add plot and its number to grob and layout_matrix
      nonnamed_plotlist <- append(nonnamed_plotlist, list(p)) # append
      layout_matrix[i, j] <- plot_iter
      # Move to the next number
      plot_iter <- plot_iter + 1
    }
  }
}

arranged_plot <- grid.arrange(
  grobs = nonnamed_plotlist,
  widths = widths,
  heights = heights,
  layout_matrix = layout_matrix
)

ggsave(
  filename = paste0(dout, "/cnv_segments_combined_", num_patients, ".png"),
  plot = arranged_plot,
  width = sum(widths) * 2,
  height = (sum(heights) + 2) * 2,
  units = "cm",
  limitsize = FALSE,
)
