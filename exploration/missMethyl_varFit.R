library(sesame)
library(minfi)
library(missMethyl)

dout <- "/projects/p30791/methylation/missMethylDiffVar"

if (!file.exists(dout)) {
  dir.create(dout)
}

## Load data
betas <- read.table(
  "/projects/p30791/methylation/sesame_out/betas_processed.csv",
  row.names = 1, sep = ",", header = TRUE
) # nrows=2000 for testing
meta <- read.csv("/projects/p30791/methylation/data/meta.csv")

betas <- as.matrix(betas)
meta <- meta[paste0("X", meta$IDAT) %in% colnames(betas), ]

refs <- c("Normal", "Normal", "CUB", "OQ", "AN")
comps <- c("AN", "CUB", "OQ", "AN", "TU")

for (i in seq_along(refs)) {
  ref <- refs[i]
  comp <- comps[i]

  ## Subset the betas and metadata
  meta_sub <- meta[meta$Sample.Region %in% c(ref, comp), ]
  meta_sub["(Intercept)"] <- 1 # needed for varFit
  meta_sub$Sample.Region.bin <- as.numeric(meta_sub$Sample.Region == comp)
  betas_sub <- betas[, paste0("X", meta_sub$IDAT)]

  ## Exclude probes that are missing levels on sample region etc.
  cpg_ok <- checkLevels(betas_sub, meta_sub$Sample.Region)
  print(paste0(sum(cpg_ok), " probes have sufficient levels for sample region."))
  betas_sub <- betas_sub[cpg_ok, ]

  ## Exclude probes that include any missing values (too strict but testing quick & dirty for now)
  keep <- apply(betas_sub, 1, function(row) all(!is.na(row)))
  betas_sub <- betas_sub[keep, ]

  fitvar <- varFit(
    betas_sub,
    design = meta_sub[c("(Intercept)", "Sample.Region.bin")],
    coef = c(1, 2) # Intercept and binary sample region
  )
  topDV <- topVar(
    fitvar,
    coef = 2,
    number = dim(betas_sub)[1] # number of top probes to include
  )
  write.csv(
    topDV,
    paste0(dout, "/topDV_", ref, "_vs_", comp, ".csv"),
    quote = FALSE, row.names = TRUE
  )
}
