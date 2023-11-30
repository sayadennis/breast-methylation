library(sesame)

## As sesame and sesameData are under active development, this documentation is
## specific to the following version of R, sesame, sesameData and ExperimentHub:
sesame_checkVersion()

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Please provide the path IDAT directory and output directory.", call.=FALSE)
} else {
  idat_dir = args[1]  # e.g. "/projects/p30791/methylation/data/IDAT_all"
  out_dir = args[2]  # e.g. "/projects/p30791/methylation/sesame_out"
}

if (!dir.exists(out_dir)) {
  dir.create(out_dir)
  print(paste0("Output directory created: ", out_dir))
} else {
  print(paste0("Output directory already exists: ", out_dir, "-- Proceeding..."))
}

#### The openSesame Pipeline ####

# idat_dir = "/projects/p30791/methylation/copied_from_b1122/data/IDAT_Files/IDAT_only" # system.file("extdata/", package = "sesameData")
betas = openSesame(idat_dir, BPPARAM = BiocParallel::MulticoreParam(8))

print(paste0("Dimensions of beta matrix: ", dim(betas)[1], "rows x ", dim(betas)[2], "columns."))
print("Snippet of beta matrix:")
print(betas[1:5,1:5])

write.csv(betas, file = paste0(out_dir, "/betas.csv"), row.names = TRUE)

#### Data Preprocessing ####

sdf_prepped = openSesame(idat_dir, prep="QCDPB", func=NULL, BPPARAM = BiocParallel::MulticoreParam(8))

print("Top rows of the first item of SDF prepped") 
print(head(sdf_prepped[[1]]))
print(paste0("Length of SDF-prepped: ", length(sdf_prepped)))

# Save SDF object as Rdata
saveRDS(sdf_prepped, file=paste0(out_dir, "/sdf_prepped.RDS"))

#### Calculate quality metrics ####

qcs = openSesame(idat_dir, prep="", func=sesameQC_calcStats, BPPARAM = BiocParallel::MulticoreParam(8))

print(paste0("Length of QC object: ", length(qcs)))

# Save QC list object as Rdata
saveRDS(qcs, file=paste0(out_dir, "/qcs.RDS"))

qcs_df = do.call(rbind, lapply(qcs, as.data.frame))
print(paste0("Dimensions of QC dataframe: ", dim(qcs_df)[1], "rows x ", dim(qcs_df)[2], "columns."))
print("Snippet of QC dataframe:")
print(head(qcs_df))

write.csv(qcs_df, file = paste0(out_dir, "/qc_metrics.csv"), row.names = TRUE)

# Rank quality score metrics against public datasets

sesameQC_rankStats(qcs[[1]], platform="EPIC")

