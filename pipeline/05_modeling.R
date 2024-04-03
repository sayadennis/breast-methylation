library(sesame)
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(ggplot2)

plotdir <- "/projects/p30791/methylation/plots"
raw_dir <- "/projects/p30791/methylation/raw_data"
sesame_dir <- "/projects/p30791/methylation/sesame_data"
dout <- "/projects/p30791/methylation/differential_methylation"

############################
#### Load and prep data ####
############################

## Load data
betas <- read.table(
  paste0(sesame_dir, "/betas_processed.csv"),
  row.names = 1, sep = ",", header = TRUE
) # nrows=2000 for testing
meta <- read.csv(paste0(raw_dir, "/meta.csv"))

betas <- as.matrix(betas)
meta <- meta[paste0("X", meta$IDAT) %in% colnames(betas), ]

##############################################################
#### Differential methylation analysis by tissue category ####
##############################################################

refs <- c("CFN", "CUB", "OQ", "AN", "CFN", "CFN", "CUB", "OQ", "AN")
comps <- c("TU", "TU", "TU", "TU", "AN", "CUB", "OQ", "AN", "TU")

for (i in seq_along(refs)) {
  ref <- refs[i]
  comp <- comps[i]

  ## Subset the betas and metadata
  meta_sub <- meta[meta$Sample.Region %in% c(ref, comp), ]
  meta_sub$Sample.Region <- relevel(factor(meta_sub$Sample.Region), ref)
  betas_sub <- betas[, paste0("X", meta_sub$IDAT)]

  ## Exclude probes that are missing levels on sample region etc.
  cpg_ok <- checkLevels(betas_sub, meta_sub$Sample.Region)
  print(paste0(sum(cpg_ok), " probes have sufficient levels for sample region."))
  betas_sub <- betas_sub[cpg_ok, ]

  ## Differential methylation analysis
  smry <- DML(betas_sub, ~ Sample.Region + Age, meta = meta_sub) # cant control for BMI
  saveRDS(
    smry,
    paste0(
      dout, "/model_summary_ref", ref, "_comp", comp, ".RDS"
    )
  )
  test_result <- summaryExtractTest(smry)
  write.csv(
    test_result,
    paste0(
      dout, "/DML_results_ref", ref, "_comp", comp, ".csv"
    ),
    quote = FALSE, row.names = FALSE
  )

  ## DMR (Differentially methylated regions)
  for (contrast in dmContrasts(smry)) {
    # merge probes to regions
    merged <- DMR(betas_sub, smry, contrast, platform = "EPIC")
    write.csv(
      merged,
      paste0(
        dout, "/DMR_results_", contrast, "_ref", ref, ".csv"
      ),
      quote = FALSE, row.names = FALSE
    )
  }
}

# nolint start
# df <- data.frame(
#   Age = meta$Age,
#   BetaValue = betas[test_result$Probe_ID[nrow(test_result)], ]
# )
#
# ggplot(df, aes(Age, BetaValue)) +
#   geom_smooth(method = "lm") +
#   geom_point()
# ggsave(paste0(plotdir, "/age_vs_betas.png"))
# nolint end

###################################################################
#### Differential methylation AN vs. TU separated by ER status ####
###################################################################

refs <- c("AN")
comps <- c("TU")

for (i in seq_along(refs)) {
  ref <- refs[i]
  comp <- comps[i]

  # Subset the betas and metadata
  meta_sub <- meta[meta$Sample.Region %in% c(ref, comp), ]
  betas_sub <- betas[, paste0("X", meta_sub$IDAT)]

  for (er_status in unique(meta_sub$ER)) {
    # Subset metadata and betas by ER status
    meta_sub_er <- meta_sub[meta_sub$ER == er_status, ]
    betas_sub_er <- betas_sub[, paste0("X", meta_sub_er$IDAT)]
    # Exclude probes that are missing levels on sample region etc.
    cpg_ok <- checkLevels(betas_sub_er, meta_sub_er$Sample.Region)
    print(paste0(sum(cpg_ok), " probes have sufficient levels for sample region."))
    betas_sub_cpgfilter <- betas_sub_er[cpg_ok, ]
    # Run DML
    smry <- DML(betas_sub_cpgfilter, ~ Sample.Region + Age, meta = meta_sub_er)
    test_result <- summaryExtractTest(smry)
    # Save results
    write.csv(
      test_result,
      paste0(
        dout, "/DML_results_ER", er_status, "_only_", ref, "_vs_", comp, ".csv"
      ),
      quote = FALSE, row.names = FALSE
    )
  }
}

#############################################################
#### Differential methylation analysis by tumor metadata ####
#############################################################

comps <- c("CUB", "OQ", "AN", "TU")

meta_cases <- meta[meta$Case.Control == "case", ]
betas_cases <- betas[, paste0("X", meta_cases$IDAT)]

# Create useful metadata columns
meta_cases$HER2 <- relevel(
  factor(
    ifelse(
      meta_cases$HER2 %in% c(2, 3),
      "Positive", "Negative"
    )
  ), "Negative"
)
meta_cases$ER <- relevel(
  factor(
    ifelse(
      meta_cases$ER == "+",
      "Positive", "Negative"
    )
  ), "Negative"
)

for (tissue_category in comps) {
  print(paste0("Running analysis for ", tissue_category, "..."))
  # Subset the betas and metadata
  meta_sub <- meta_cases[meta_cases$Sample.Region == tissue_category, ]
  betas_sub <- betas_cases[, paste0("X", meta_sub$IDAT)]

  for (contrast_feature in c("ER", "HER2")) {
    print(paste0("Running ", contrast_feature, "..."))
    # Exclude probes that are missing levels on sample region etc.
    cpg_ok <- checkLevels(betas_sub, meta_sub[[contrast_feature]])
    print(paste0(sum(cpg_ok), " probes have sufficient levels for sample region."))
    betas_sub_cpgfilter <- betas_sub[cpg_ok, ]
    # Create formula string with contrast feature to plug into DML
    formula_string <- paste("~", contrast_feature, "+ Age")
    formula_expression <- as.formula(formula_string)
    # Run DML
    smry <- DML(betas_sub_cpgfilter, formula_expression, meta = meta_sub)
    test_result <- summaryExtractTest(smry)
    # Save results
    write.csv(
      test_result,
      paste0(
        dout, "/DML_results_of_", tissue_category, "_by_", contrast_feature, ".csv"
      ),
      quote = FALSE, row.names = FALSE
    )
  }
}
