source("iEVORA.R")

din <- "/projects/p30791/methylation/sesame_out"
dout <- "/projects/p30791/methylation/sesame_out/differential_methylation"

## Load data
betas <- read.table(
  paste0(din, "/betas_processed.csv"),
  row.names = 1, sep = ",", header = TRUE
) # nrows=2000 for testing
meta <- read.csv("/projects/p30791/methylation/data/meta.csv")

betas <- as.matrix(betas)
meta <- meta[paste0("X", meta$IDAT) %in% colnames(betas), ]

## Keep only CFN and AN
meta_sub <- meta[meta$Sample.Region %in% c("Normal", "AN"), ]
betas_sub <- betas[, paste0("X", meta_sub$IDAT)]

## Run DVMC
data.m <- betas_sub
pheno.v <- as.numeric(as.vector(meta_sub$Sample.Region == "AN"))

# Remove probes that are all NaN in either tissue category
notnan <- !is.na(data.m)
remove <- (apply(notnan[, pheno.v == 1], 1, sum) < 2) | (apply(notnan[, pheno.v == 0], 1, sum) < 2)
data.m <- data.m[!remove, ]

topDVMC.m <- iEVORA(data.m, pheno.v)

write.csv(
  topDVMC.m,
  paste0(dout, "/topDVMC.csv"),
  quote = FALSE,
  row.names = TRUE
)
