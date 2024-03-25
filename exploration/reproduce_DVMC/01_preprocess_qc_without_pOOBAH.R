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

#### Data Preprocessing ####

print("Running openSesame() to obtain processed beta values WITHOUT pOOBAH...")
betas_proc <- openSesame(
  idat_dir,
  prep = "QCDB",
  func = getBetas,
)
print("Successfully read processed betas!")

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

print("Calculating processed quality metrics...")
qcs_proc <- openSesame(
  idat_dir,
  prep = "QCDPB",
  func = sesameQC_calcStats,
)

qcs_proc_df <- do.call(rbind, lapply(qcs_proc, as.data.frame))

print("Successfully calculated quality metrics! Writing...")

write.csv(
  qcs_proc_df,
  file = paste0(out_dir, "/qc_processed_without_pOOBAH.csv"),
  row.names = TRUE
)

# Remove samples with poor quality metrics
print("Finding samples with low rate of successful detection...")
rm_criteria <- qcs_proc_df$frac_dt < 0.9
rm_sampleids <- rownames(qcs_proc_df[rm_criteria, ])
print("Samples that will be removed:")
print(rm_sampleids)
writeLines(rm_sampleids, paste0(out_dir, "/exclude_IDATs_without_pOOBAH.txt"))
print("Wrote samples to exclude to TXT! Removing from betas object...")

betas_proc <- betas_proc[, !(colnames(betas_proc) %in% rm_sampleids)]

print("Removed samples from betas! Writing...")

# Save betas object
write.csv(
  betas_proc,
  file = paste0(out_dir, "/betas_processed_without_pOOBAH.csv"),
  row.names = TRUE
)

print(paste0("Wrote processed betas to ", out_dir))
