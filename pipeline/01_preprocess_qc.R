library(sesame)
library(ggplot2)

## As sesame and sesameData are under active development, this documentation is
## specific to the following version of R, sesame, sesameData and ExperimentHub:
sesame_checkVersion()

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(paste("Expected two commandArgs: Please provide the IDAT directory path",
    "and output directory path.",
    call. = FALSE
  ))
} else {
  idat_dir <- args[1] # e.g. "/projects/p30791/methylation/data/IDAT_all"
  meta_fn <- args[2] # e.g. "/projects/p30791/methylation/data/meta.csv"
  out_dir <- args[3] # e.g. "/projects/p30791/methylation/sesame_out"
}

if (!dir.exists(out_dir)) {
  dir.create(out_dir)
  print(paste0("Output directory created: ", out_dir))
} else {
  print(paste0("Output directory already exists: ", out_dir, ". Proceeding..."))
}

meta <- read.csv(meta_fn)

#### The openSesame Pipeline ####

print("Running openSesame() to obtain raw beta values...")
betas_raw <- openSesame(
  idat_dir,
  func = getBetas,
)
print("Successfully read raw betas! Now reading SDFs...")
sdf_raw <- openSesame(
  idat_dir,
  func = NULL,
)
print("Successfully read raw SDFs!")

print(paste0(
  "Dimensions of raw beta matrix: ",
  dim(betas_raw)[1], " rows x ", dim(betas_raw)[2], " columns."
))
print(paste0(
  "Selecting IDATs that are found in the metadata file: ", meta_fn, "..."
))
betas_raw <- betas_raw[, colnames(betas_raw) %in% meta$IDAT]
print(paste0(
  "Dimensions of raw beta matrix following sample filter: ",
  dim(betas_raw)[1], " rows x ", dim(betas_raw)[2], " columns."
))

write.csv(betas_raw, file = paste0(out_dir, "/betas_raw.csv"), row.names = TRUE)
saveRDS(sdf_raw, file = paste0(out_dir, "/sdf_raw.RDS"))
print(paste0("Wrote raw beta and SDF data to ", out_dir))

#### Data Preprocessing ####

print("Running openSesame() to obtain processed beta values...")
betas_proc <- openSesame(
  idat_dir,
  prep = "QCDPB",
  func = getBetas,
)
print("Successfully read processed betas! Now reading SDFs...")
sdf_proc <- openSesame(
  idat_dir,
  prep = "QCDPB",
  func = NULL,
)
print("Successfully read processed SDFs!")

print(paste0(
  "Dimensions of processed beta matrix: ",
  dim(betas_proc)[1], " rows x ", dim(betas_proc)[2], " columns."
))
betas_proc <- betas_proc[, colnames(betas_proc) %in% meta$IDAT]
print(paste0(
  "Dimensions of processed beta matrix following sample filter: ",
  dim(betas_proc)[1], " rows x ", dim(betas_proc)[2], " columns."
))

#### Calculate quality metrics ####

print("Calculating raw and processed quality metrics...")
qcs_raw <- openSesame(
  idat_dir,
  prep = "",
  func = sesameQC_calcStats,
)
qcs_proc <- openSesame(
  idat_dir,
  prep = "QCDPB",
  func = sesameQC_calcStats,
)

qcs_raw_df <- do.call(rbind, lapply(qcs_raw, as.data.frame))
qcs_proc_df <- do.call(rbind, lapply(qcs_proc, as.data.frame))

print("Successfully calculated quality metrics! Writing...")

# Save QC metrics
saveRDS(qcs_raw, file = paste0(out_dir, "/qcs_raw.RDS"))
saveRDS(qcs_proc, file = paste0(out_dir, "/qcs_processed.RDS"))

write.csv(
  qcs_raw_df,
  file = paste0(out_dir, "/qc_raw.csv"),
  row.names = TRUE
)
write.csv(
  qcs_proc_df,
  file = paste0(out_dir, "/qc_processed.csv"),
  row.names = TRUE
)

# Remove samples with poor quality metrics (Dye bias >1.5 or <0.5)
print("Finding samples with high distortion...")
rm_criteria <- qcs_raw_df$frac_dt < 0.9
rm_sampleids <- rownames(qcs_raw_df[rm_criteria, ])
print("Samples that will be removed:")
print(rm_sampleids)
writeLines(rm_sampleids, paste0(out_dir, "/exclude_IDATs.txt"))
print("Wrote samples to exclude to TXT! Removing from beta and SDF objects...")

sdf_proc <- sdf_proc[!(names(sdf_proc) %in% rm_sampleids)]
betas_proc <- betas_proc[, !(colnames(betas_proc) %in% rm_sampleids)]

print("Removed samples from betas and SDFs! Writing...")

# Save SDF object as Rdata
write.csv(
  betas_proc,
  file = paste0(out_dir, "/betas_processed.csv"),
  row.names = TRUE
)
saveRDS(sdf_proc, file = paste0(out_dir, "/sdf_processed.RDS"))

print(paste0("Wrote processed beta and SDF data to ", out_dir))

## Plot metagene statistic

KYCG_plotMeta_custom <- function(betas, platform = "EPIC", alpha = 0.2) {
  if (!is.matrix(betas)) {
    betas <- cbind(sample = betas)
  }
  stopifnot(!is.null(platform))

  dbs <- KYCG_getDBs(sprintf("%s.metagene", platform))
  df <- dbStats(betas, dbs, long = TRUE)
  dflabel <- data.frame(
    ord = as.integer(names(dbs)),
    reg = vapply(dbs, function(x) attr(x, "label"), character(1))
  )

  ggplot(df) +
    annotate("rect",
      xmin = -1, xmax = 10, ymin = -Inf,
      ymax = Inf, fill = "grey80", alpha = .5, color = NA
    ) +
    geom_line(aes_string("db", "value", group = "query", alpha = alpha)) +
    scale_x_continuous(breaks = dflabel$ord, labels = dflabel$reg) +
    ylab("Mean DNA Methylation Level") +
    xlab("") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}


plotmeta <- KYCG_plotMeta_custom(
  betas_proc,
  platform = "EPIC",
  alpha = 0.1
) + ggtitle("Metagene Statistics of All Samples")
file_path <- paste0(out_dir, "/plotmeta.png")
ggsave(file_path, plot = plotmeta)
