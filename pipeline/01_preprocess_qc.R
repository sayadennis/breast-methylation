library(sesame)

## As sesame and sesameData are under active development, this documentation is
## specific to the following version of R, sesame, sesameData and ExperimentHub:
sesame_checkVersion()

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Expected two commandArgs: Please provide the IDAT directory path and output directory path.", call.=FALSE)
} else {
  idat_dir = args[1]  # e.g. "/projects/p30791/methylation/data/IDAT_all"
  out_dir = args[2]  # e.g. "/projects/p30791/methylation/sesame_out"
}

if (!dir.exists(out_dir)) {
  dir.create(out_dir)
  print(paste0("Output directory created: ", out_dir))
} else {
  print(paste0("Output directory already exists: ", out_dir, " -- Proceeding..."))
}

#### The openSesame Pipeline ####

# idat_dir = "/projects/p30791/methylation/copied_from_b1122/data/IDAT_Files/IDAT_only" # system.file("extdata/", package = "sesameData")
print("Running openSesame() to obtain raw beta values...")
betas_raw = openSesame(idat_dir, func=getBetas, BPPARAM = BiocParallel::MulticoreParam(8))
print("Successfully read raw betas! Now reading SDFs...")
sdf_raw = openSesame(idat_dir, func=NULL, BPPARAM = BiocParallel::MulticoreParam(8))
print("Successfully read raw SDFs!")

print(paste0("Dimensions of raw beta matrix: ", dim(betas_raw)[1], "rows x ", dim(betas_raw)[2], "columns."))

write.csv(betas_raw, file = paste0(out_dir, "/betas_raw.csv"), row.names = TRUE)
saveRDS(sdf_raw, file=paste0(out_dir, "/sdf_raw.RDS"))
print(paste0("Wrote raw beta and SDF data to ", out_dir))

#### Data Preprocessing ####

print("Running openSesame() to obtain processed beta values...")
betas_proc = openSesame(idat_dir, prep="QCDPB", func=getBetas, BPPARAM = BiocParallel::MulticoreParam(8))
print("Successfully read processed betas! Now reading SDFs...")
sdf_proc = openSesame(idat_dir, prep="QCDPB", func=NULL, BPPARAM = BiocParallel::MulticoreParam(8))
print("Successfully read processed SDFs!")

#### Calculate quality metrics ####

print("Calculating raw and processed quality metrics...")
qcs_raw = openSesame(idat_dir, prep="", func=sesameQC_calcStats, BPPARAM = BiocParallel::MulticoreParam(8))
qcs_proc = openSesame(idat_dir, prep="QCDPB", func=sesameQC_calcStats, BPPARAM = BiocParallel::MulticoreParam(8))

qcs_raw_df = do.call(rbind, lapply(qcs_raw, as.data.frame))
qcs_proc_df = do.call(rbind, lapply(qcs_proc, as.data.frame))

print("Successfully calculated quality metrics! Writing...")

# Save QC metrics
saveRDS(qcs_raw, file=paste0(out_dir, "/qcs_raw.RDS"))
saveRDS(qcs_proc, file=paste0(out_dir, "/qcs_processed.RDS"))

write.csv(qcs_raw_df, file = paste0(out_dir, "/qc_raw.csv"), row.names = TRUE)
write.csv(qcs_proc_df, file = paste0(out_dir, "/qc_processed.csv"), row.names = TRUE)

# Remove samples with extremely poor raw quality metrics
print("Finding samples with high distortion...")
distorted_sampleids = rownames(qcs_raw_df[abs(qcs_raw_df$RGdistort-1)>0.5,])  # Dye bias >1.5 or <0.5
print("Samples that will be removed:")
print(distorted_sampleids)
writeLines(distorted_sampleids, paste0(out_dir, "/exclude_IDATs.txt"))
print("Wrote samples to exclude to TXT! Removing from beta and SDF objects...")

sdf_proc = sdf_proc[!(names(sdf_proc) %in% distorted_sampleids)]
betas_proc = betas_proc[,!(colnames(betas_proc) %in% distorted_sampleids)]

print("Removed samples from betas and SDFs! Writing...")

# Save SDF object as Rdata
write.csv(betas_proc, file = paste0(out_dir, "/betas_processed.csv"), row.names = TRUE)
saveRDS(sdf_proc, file=paste0(out_dir, "/sdf_processed.RDS"))

print(paste0("Wrote processed beta and SDF data to ", out_dir))
