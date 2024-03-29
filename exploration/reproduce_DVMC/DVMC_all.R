source("iEVORA.R")

din <- "/projects/p30791/methylation/sesame_out"
dout <- "/projects/p30791/methylation/sesame_out/differential_variability"
probemeta <- "/projects/p30791/methylation/illumina_array_meta"

## Load data
betas <- read.table(
  paste0(din, "/betas_processed.csv"),
  row.names = 1, sep = ",", header = TRUE
) # nrows=2000 for testing
meta <- read.csv("/projects/p30791/methylation/data/meta.csv")

betas <- as.matrix(betas)
meta <- meta[paste0("X", meta$IDAT) %in% colnames(betas), ]

#################################################
#### Loop through comparisons and run iEVORA ####
#################################################

refs <- c("Normal", "Normal", "CUB", "OQ", "AN")
comps <- c("AN", "CUB", "OQ", "AN", "TU")

for (i in seq_along(refs)) {
  ref <- refs[i]
  comp <- comps[i]

  ## Create subset data with categories of interest
  meta_sub <- meta[meta$Sample.Region %in% c(ref, comp), ]
  betas_sub <- betas[, paste0("X", meta_sub$IDAT)]

  ## Create arguments for DVMC
  data.m <- betas_sub
  pheno.v <- as.numeric(as.vector(meta_sub$Sample.Region == comp))

  # Remove probes that are all NaN in either tissue category
  notnan <- !is.na(data.m)
  remove <- (apply(notnan[, pheno.v == 1], 1, sum) < 2) |
    (apply(notnan[, pheno.v == 0], 1, sum) < 2)
  data.m <- data.m[!remove, ]

  ## Run iEVORA
  DVMC.m <- iEVORA(data.m, pheno.v)

  # Write all results to CSV
  write.csv(
    DVMC.m,
    paste0(dout, "/DVMC_", ref, "_vs_", comp, ".csv"),
    quote = FALSE,
    row.names = TRUE
  )

  # Write DV-signidicant probe IDs to TXT
  writeLines(
    rownames(DVMC.m[(DVMC.m[, "q(BT)"] < 0.05) & (DVMC.m[, "P(TT)"] < 0.05), ]),
    paste0(dout, "/DV_iEVORA_", ref, "_vs_", comp, ".txt")
  )
}

########################################################
#### Recreate Normal vs. AN results with just HM450 ####
########################################################

hm450 <- read.csv(
  paste0(probemeta, "/humanmethylation450_15017482_v1-2.csv"),
  skip = 7
)
epic <- read.csv(
  paste0(probemeta, "/EPIC-8v2-0_A1.csv"),
  skip = 7
)

# Select only CpG probes to avoid irregularities
hm450 <- unique(hm450[grep("^cg", hm450$Name), c("Name", "Forward_Sequence")])
epic <- unique(epic[grep("^cg", epic$Name), c("Name", "Forward_Sequence")])

overlap <- intersect(hm450$Name, epic$Name)

epic <- epic[(epic$Name != "cg18801637"), ] # weird duplicate probe name

## Are the probe designs consistent?
rownames(hm450) <- hm450$Name
rownames(epic) <- epic$Name

hm450_seqs <- hm450[overlap, "Forward_Sequence"]
epic_seqs <- epic[overlap, "Forward_Sequence"]

inconsistent_ct <- sum(!(hm450_seqs == epic_seqs))
inconsistent_pct <- 100 * inconsistent_ct / length(overlap)
print(paste(
  "####", inconsistent_ct, "probes out of", length(overlap),
  "(", round(inconsistent_pct, 2), "% ) have different probe sequences."
))

ref <- "Normal"
comp <- "AN"

## Create subset data with categories of interest
meta_sub <- meta[meta$Sample.Region %in% c(ref, comp), ]
betas_sub <- betas[, paste0("X", meta_sub$IDAT)]

## Create arguments for DVMC
data.m <- betas_sub
pheno.v <- as.numeric(as.vector(meta_sub$Sample.Region == comp))

# Remove probes that are all NaN in either tissue category
notnan <- !is.na(data.m)
remove <- (apply(notnan[, pheno.v == 1], 1, sum) < 2) |
  (apply(notnan[, pheno.v == 0], 1, sum) < 2)
data.m <- data.m[!remove, ]

## Select only HM450 probes
data.m.hm450 <- data.m[rownames(data.m) %in% overlap, ]

## Run iEVORA
DVMC.m <- iEVORA(data.m.hm450, pheno.v)

write.csv(
  DVMC.m,
  paste0(dout, "/DVMC_HM450.csv"),
  quote = FALSE,
  row.names = TRUE
)
