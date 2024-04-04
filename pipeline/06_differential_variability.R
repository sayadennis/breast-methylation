library(sesame)
library(minfi)
library(missMethyl)

dout <- "/projects/p30791/methylation/differential_variability/missMethylDiffVar"

if (!file.exists(dout)) {
  dir.create(dout)
}

## Load data
betas <- read.table(
  "/projects/p30791/methylation/sesame_data/betas_processed.csv",
  row.names = 1, sep = ",", header = TRUE
) # nrows=2000 for testing
meta <- read.csv("/projects/p30791/methylation/raw_data/meta.csv")

betas <- as.matrix(betas)
meta <- meta[paste0("X", meta$IDAT) %in% colnames(betas), ]

refs <- c("CFN", "CUB", "OQ", "CFN", "CFN", "CUB", "OQ", "AN")
comps <- c("TU", "TU", "TU", "AN", "CUB", "OQ", "AN", "TU")

for (i in seq_along(refs)) {
  ref <- refs[i]
  comp <- comps[i]
  print(paste0("Running for comparison with reference ", ref, " and comparison ", comp, "..."))

  ## Subset the betas and metadata
  meta_sub <- meta[meta$Sample.Region %in% c(ref, comp), ]
  meta_sub["(Intercept)"] <- 1 # needed for varFit
  meta_sub$Sample.Region.bin <- as.numeric(meta_sub$Sample.Region == comp)
  betas_sub <- betas[, paste0("X", meta_sub$IDAT)]
  print("Finished aligning metadata and betas for this comparison!")

  ## Exclude probes that are missing levels on sample region etc.
  cpg_ok <- checkLevels(betas_sub, meta_sub$Sample.Region)
  print(paste0(sum(cpg_ok), " probes have sufficient levels for sample region."))
  betas_sub <- betas_sub[cpg_ok, ]

  ## Exclude probes that include any missing values (too strict but testing quick & dirty for now)
  keep <- apply(betas_sub, 1, function(row) all(!is.na(row)))
  betas_sub <- betas_sub[keep, ]
  print("Finished narrowing down CpGs of interest!")

  fitvar <- varFit(
    betas_sub,
    design = meta_sub[c("(Intercept)", "Sample.Region.bin")],
    coef = c(1, 2) # Intercept and binary sample region
  )
  print("Finished running varFit! Selecting top DVs...")
  topDV <- topVar(
    fitvar,
    coef = 2,
    number = dim(betas_sub)[1] # number of top probes to include
  )
  print("Writing top DV results to CSV...")
  write.csv(
    topDV,
    paste0(dout, "/topDV_", ref, "_vs_", comp, ".csv"),
    quote = FALSE, row.names = TRUE
  )
  print("Done! Writing top DV probe IDs to TXT...")
  writeLines(
    rownames(topDV[(topDV["Adj.P.Value"] < 0.05) & (topDV["LogVarRatio"] > 0.2), ]),
    paste0(dout, "/hyperDV_missMethyl_", ref, "_vs_", comp, ".txt")
  )
  writeLines(
    rownames(topDV[(topDV["Adj.P.Value"] < 0.05) & (topDV["LogVarRatio"] < -0.2), ]),
    paste0(dout, "/hypoDV_missMethyl_", ref, "_vs_", comp, ".txt")
  )
}
